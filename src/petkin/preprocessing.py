"""Native PET/MRI preprocessing, implemented with NumPy, SciPy and SimpleITK.

No PETPrep executable or source code is used. Transform files map output-grid
LPS millimetres to source-image LPS millimetres (ITK's resampling convention).
"""
import csv
import html
import json
from pathlib import Path

import numpy as np
import SimpleITK as sitk
from scipy.ndimage import gaussian_filter
from scipy.optimize import nnls

from .bids import validate_pet_metadata
from .imaging import (split_dynamic, preprocess_dynamic, apply_composed_transforms,
                      _duration_mean)
from .reporting import file_manifest


def _same_grid(a, b):
    return a.GetDimension() == b.GetDimension() and a.GetSize() == b.GetSize() and all(
        np.allclose(getattr(a, k)(), getattr(b, k)(), atol=1e-6, rtol=0)
        for k in ('GetSpacing', 'GetOrigin', 'GetDirection'))


def _check_image(image, name):
    if image.GetDimension() != 3 or image.GetNumberOfComponentsPerPixel() != 1:
        raise ValueError(name + ' must be a scalar 3D image')
    if not np.isfinite(sitk.GetArrayViewFromImage(image)).all():
        raise ValueError(name + ' contains nonfinite voxels')


def _mask_on_grid(mask, image):
    _check_image(mask, 'Brain mask')
    if not _same_grid(mask, image):
        raise ValueError('Brain mask must share the MRI physical grid')
    a = sitk.GetArrayFromImage(mask)
    if not np.all((a == 0) | (a == 1)) or np.count_nonzero(a) < 30:
        raise ValueError('Brain mask must be binary and contain at least 30 voxels')
    return sitk.Cast(mask, sitk.sitkUInt8)


def n4_correct(mri, mask, iterations=(30, 20, 10), shrink=2):
    """Estimate a smooth MRI bias field inside a supplied anatomical brain mask."""
    _check_image(mri, 'MRI'); mask = _mask_on_grid(mask, mri)
    values = sitk.GetArrayViewFromImage(mri)[sitk.GetArrayViewFromImage(mask) > 0]
    if np.any(values < 0) or np.mean(values > 0) < .95 or np.std(values) <= 1e-8:
        raise ValueError('MRI brain mask must contain nonnegative, nonconstant intensities with at least 95% positive voxels')
    if not iterations or any(int(i) != i or i < 1 for i in iterations):
        raise ValueError('Positive integer N4 iteration counts required')
    if not isinstance(shrink, int) or shrink < 1:
        raise ValueError('N4 shrink factor must be a positive integer')
    image = sitk.Cast(mri, sitk.sitkFloat32)
    # N4 estimates its log field only from positive voxels. Isolated zeros at
    # mask boundaries remain zero in the corrected image and do not alter masks.
    fit_mask = sitk.Cast((mask > 0) & (image > 0), sitk.sitkUInt8)
    n4 = sitk.N4BiasFieldCorrectionImageFilter(); n4.SetNumberOfThreads(1)
    n4.SetMaximumNumberOfIterations(list(iterations))
    n4.Execute(sitk.Shrink(image, [shrink]*3), sitk.Shrink(fit_mask, [shrink]*3))
    bias = sitk.Exp(n4.GetLogBiasFieldAsImage(image))
    corrected = image / bias
    return corrected, bias


def segment_t1(mri, mask, iterations=100):
    """Three Gaussian intensity classes, ordered as candidate CSF, GM, WM.

    This is an unsupervised T1 intensity model, not an anatomical parcellation.
    Pathology, contrast, incomplete masks and intensity inversion require review.
    """
    if not isinstance(iterations, int) or iterations < 1: raise ValueError('Positive segmentation iteration count required')
    _check_image(mri, 'MRI'); mask = _mask_on_grid(mask, mri)
    inside = sitk.GetArrayFromImage(mask).astype(bool)
    data = sitk.GetArrayFromImage(mri).astype(float); y = data[inside]
    scale = float(np.std(y))
    if scale < 1e-8: raise ValueError('Tissue segmentation requires intensity variation')
    center = float(np.mean(y)); x = (y-center)/scale
    # Equal-frequency quantiles can initialize two classes in the same dominant
    # tissue. Spread the initial centers across the robust intensity range.
    means = np.linspace(*np.quantile(x, [.01, .99]), 3); variances = np.full(3, .25)
    weights = np.full(3, 1/3)
    previous = -np.inf
    for _ in range(iterations):
        logp = np.log(weights) - .5*np.log(2*np.pi*variances) - .5*(x[:, None]-means)**2/variances
        maximum = logp.max(axis=1, keepdims=True)
        p = np.exp(logp-maximum); normalizer = p.sum(axis=1, keepdims=True)
        likelihood = float(np.sum(maximum + np.log(normalizer)))
        p /= normalizer
        count = p.sum(axis=0).clip(min=1e-12)
        weights = count / len(x); means = (p*x[:, None]).sum(axis=0)/count
        variances = np.maximum((p*(x[:, None]-means)**2).sum(axis=0)/count, 1e-4)
        if abs(likelihood-previous) <= 1e-7*len(x): break
        previous = likelihood
    # Recompute posteriors with the final parameters before assigning classes.
    order = np.argsort(means); means = means[order]; variances = variances[order]; weights = weights[order]
    logp = np.log(weights) - .5*np.log(variances) - .5*(x[:, None]-means)**2/variances
    p = np.exp(logp-logp.max(axis=1, keepdims=True)); p /= p.sum(axis=1, keepdims=True)
    labels = np.zeros(data.shape, np.uint8); labels[inside] = 1+p.argmax(axis=1)
    label_image = sitk.GetImageFromArray(labels); label_image.CopyInformation(mri)
    probabilities = []
    for k in range(3):
        a = np.zeros(data.shape, np.float32); a[inside] = p[:, k]
        image = sitk.GetImageFromArray(a); image.CopyInformation(mri); probabilities.append(image)
    qc = dict(class_names=['candidate_CSF', 'candidate_GM', 'candidate_WM'],
              intensity_means=(means*scale+center).tolist(), fractions=weights.tolist(),
              review_required=True, limitation='Intensity classes are not expert tissue labels')
    return label_image, probabilities, qc


def normalize_anatomy(fixed, moving, *, deformable=True, iterations=60):
    """Estimate fixed-to-moving affine plus optional smooth B-spline mapping.

    A positive sampled Jacobian is a numerical screen, not anatomical accuracy.
    No inverse transform is approximated or silently substituted.
    """
    if not isinstance(iterations, int) or iterations < 1: raise ValueError('Positive registration iteration count required')
    for image in (fixed, moving): _check_image(image, 'Anatomical image')
    if min(np.std(sitk.GetArrayViewFromImage(i)) for i in (fixed, moving)) < 1e-8:
        raise ValueError('Normalization requires nonconstant images')
    a, b = [sitk.Cast(i, sitk.sitkFloat32) for i in (fixed, moving)]
    initial = sitk.CenteredTransformInitializer(a, b, sitk.AffineTransform(3),
                                               sitk.CenteredTransformInitializerFilter.GEOMETRY)
    reg = sitk.ImageRegistrationMethod(); reg.SetNumberOfThreads(1)
    reg.SetMetricAsMattesMutualInformation(40); reg.SetMetricSamplingStrategy(reg.NONE)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetShrinkFactorsPerLevel([4, 2, 1]); reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    reg.SetOptimizerAsRegularStepGradientDescent(1., .001, iterations)
    reg.SetOptimizerScalesFromPhysicalShift(); reg.SetInitialTransform(initial, inPlace=False)
    affine = reg.Execute(a, b)
    qc = dict(affine_metric=float(reg.GetMetricValue()), affine_stop=reg.GetOptimizerStopConditionDescription())
    transform = sitk.CompositeTransform([affine])
    if deformable:
        mesh = [max(1, int(round((n-1)*s/40))) for n, s in zip(a.GetSize(), a.GetSpacing())]
        bspline = sitk.BSplineTransformInitializer(a, mesh, order=3)
        # Start a new optimizer so affine parameter scaling is not inherited by
        # LBFGSB, whose B-spline parameters already represent millimetres.
        reg = sitk.ImageRegistrationMethod(); reg.SetNumberOfThreads(1)
        reg.SetMetricAsMattesMutualInformation(40); reg.SetMetricSamplingStrategy(reg.NONE)
        reg.SetInterpolator(sitk.sitkLinear)
        reg.SetShrinkFactorsPerLevel([2, 1]); reg.SetSmoothingSigmasPerLevel([1, 0])
        reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
        reg.SetMovingInitialTransform(affine)
        reg.SetInitialTransform(bspline, inPlace=False)
        reg.SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-5,
                                numberOfIterations=iterations, maximumNumberOfCorrections=5,
                                maximumNumberOfFunctionEvaluations=iterations*10)
        warp = reg.Execute(a, b); transform.AddTransform(warp)
        stop = reg.GetOptimizerStopConditionDescription()
        qc.update(deformable_metric=float(reg.GetMetricValue()), deformable_stop=stop,
                  nonlinear_stop_requires_review='tolerance' not in stop.lower() and 'convergence' not in stop.lower())
    field = sitk.TransformToDisplacementField(transform, sitk.sitkVectorFloat64,
                                             a.GetSize(), a.GetOrigin(), a.GetSpacing(), a.GetDirection())
    jacobian = sitk.GetArrayFromImage(sitk.DisplacementFieldJacobianDeterminant(field))
    qc.update(minimum_jacobian=float(jacobian.min()), maximum_jacobian=float(jacobian.max()),
              convention='Fixed/output to moving/source physical LPS mm', review_required=True)
    if not np.isfinite(jacobian).all() or np.any(jacobian <= 0):
        raise RuntimeError('Anatomical normalization has nonpositive/nonfinite sampled Jacobians')
    transform.FlattenTransform()
    return transform, qc


def gtm_correct(dynamic, labels, fwhm_mm, *, maximum_condition=1e6):
    """Regional geometric transfer matrix correction with nonnegative fitting.

    Assumes a spatially invariant Gaussian PSF and homogeneous activity within
    each supplied label. Label zero is an explicit, fitted FOV-background region.
    Outputs regional estimates, not a fabricated voxelwise corrected PET image.
    """
    frames = split_dynamic(dynamic); reference = frames[0]
    if not _same_grid(labels, reference): raise ValueError('GTM labels must match the PET grid')
    lab = sitk.GetArrayFromImage(labels)
    if not np.isfinite(lab).all() or np.any(lab < 0) or not np.all(lab == np.round(lab)):
        raise ValueError('GTM labels must be nonnegative integers')
    ids = np.unique(lab).astype(int)
    if not np.any(ids > 0): raise ValueError('GTM requires positive anatomical labels')
    if len(ids) > 512: raise ValueError('GTM supports at most 512 regions including background')
    fwhm = np.asarray(fwhm_mm, float)
    if fwhm.shape != (3,) or not np.isfinite(fwhm).all() or np.any(fwhm <= 0):
        raise ValueError('GTM needs three positive finite PSF FWHM values in x,y,z mm')
    sigma = (fwhm / np.sqrt(8*np.log(2)) / np.array(reference.GetSpacing()))[::-1]
    masks = [lab == i for i in ids]; matrix = np.empty((len(ids), len(ids)))
    for j, mask in enumerate(masks):
        blurred = gaussian_filter(mask.astype(float), sigma, mode='constant', cval=0, truncate=4.)
        matrix[:, j] = [blurred[m].mean() for m in masks]
    condition = float(np.linalg.cond(matrix))
    if not np.isfinite(condition) or condition > maximum_condition:
        raise ValueError('GTM matrix is ill-conditioned; revise regions or PSF')
    observed = np.array([[sitk.GetArrayViewFromImage(f)[m].mean() for f in frames] for m in masks])
    if not np.isfinite(observed).all(): raise ValueError('Nonfinite GTM observations')
    fitted = []; residuals = []
    for column in observed.T:
        solution, residual = nnls(matrix, column); fitted.append(solution); residuals.append(float(residual))
    return dict(region_ids=ids.tolist(), corrected=np.array(fitted).T, observed=observed,
                matrix=matrix, condition_number=condition, residual_norms=residuals,
                psf_fwhm_xyz_mm=fwhm.tolist(),
                assumptions=['Uniform activity per label including label-zero FOV background',
                             'Spatially invariant Gaussian PSF; zero activity outside the image FOV',
                             'Nonnegative regional estimates; uncertainty is not propagated'])


def preprocess_native(pet_path, mri_path, timing_path, output_dir, *, brain_mask=None,
                      labels=None, template=None, template_mask=None, atlas_labels=None,
                      normalize=False, affine_only=False, reference_index=None,
                      preserve_blank=False, pvc_fwhm=None, n4_iterations=(30, 20, 10)):
    """Run native anatomical preprocessing, motion correction and TAC extraction.

    Supply an MRI-grid brain mask or an anatomical template and its brain mask.
    Automatic template-mask propagation remains a reviewed research estimate.
    Original PET is resampled once for each final anatomical/standard output.
    """
    if labels and atlas_labels: raise ValueError('Supply MRI labels or template atlas labels, not both')
    if not brain_mask and not (template and template_mask):
        raise ValueError('Supply --brain-mask or both --template and --template-mask')
    if (normalize or atlas_labels) and not template: raise ValueError('Template required for normalization or atlas labels')
    if template_mask and not template: raise ValueError('Template image required with template mask')
    if pvc_fwhm is not None and not (labels or atlas_labels): raise ValueError('PVC requires anatomical region labels')
    if pvc_fwhm is not None:
        psf = np.asarray(pvc_fwhm, float)
        if psf.shape != (3,) or not np.isfinite(psf).all() or np.any(psf <= 0):
            raise ValueError('PVC requires three positive finite PSF FWHM values in x,y,z mm')
    inputs = [p for p in (pet_path, mri_path, timing_path, brain_mask, labels, template, template_mask, atlas_labels) if p]
    out = Path(output_dir).resolve()
    if any(out == Path(p).resolve() or out in Path(p).resolve().parents for p in inputs):
        raise ValueError('Output directory must not contain source inputs')
    if out.exists() and (not out.is_dir() or any(out.iterdir())): raise ValueError('Output directory must be empty')
    pet = sitk.ReadImage(str(pet_path)); mri = sitk.ReadImage(str(mri_path)); _check_image(mri, 'MRI')
    if pet.GetNumberOfComponentsPerPixel() != 1: raise ValueError('PET must be a scalar 4D image')
    metadata = json.loads(Path(timing_path).read_text()); checked = validate_pet_metadata(metadata)
    parts = split_dynamic(pet)
    if len(parts) != len(checked['frames'].start): raise ValueError('PET frame count does not match timing')
    if any(not np.isfinite(sitk.GetArrayViewFromImage(f)).all() for f in parts): raise ValueError('PET contains nonfinite values')
    atlas = sitk.ReadImage(str(template)) if template else None
    if atlas is not None: _check_image(atlas, 'Template')
    mask = sitk.ReadImage(str(brain_mask)) if brain_mask else None
    region_labels = sitk.ReadImage(str(labels)) if labels else None
    if region_labels is not None: _validate_labels(region_labels, mri)
    template_lab = sitk.ReadImage(str(atlas_labels)) if atlas_labels else None
    if template_lab is not None: _validate_labels(template_lab, atlas)
    template_brain_mask = _mask_on_grid(sitk.ReadImage(str(template_mask)), atlas) if template_mask else None
    if mask is not None: mask = _mask_on_grid(mask, mri)
    out.mkdir(parents=True, exist_ok=True)
    from . import __version__
    import scipy
    qc = dict(status='running', fully_validated=False, provenance=file_manifest(inputs),
              software=dict(petkin=__version__, numpy=np.__version__, scipy=scipy.__version__, SimpleITK=sitk.Version_VersionString()),
              settings=dict(normalize=normalize, affine_only=affine_only, reference_index=reference_index,
                            preserve_blank=preserve_blank, pvc_fwhm=None if pvc_fwhm is None else psf.tolist(), n4_iterations=list(n4_iterations)),
              backend='Native NumPy/SciPy/SimpleITK; no PETPrep execution',
              limitations=['Rigid interframe correction; no intraframe or attenuation-map correction',
                           'Brain masks, tissue classes, parcellations and alignments require review',
                           'No independent measured-human accuracy claim'])
    def save_status(): (out/'preprocessing.json').write_text(json.dumps(qc, indent=2, allow_nan=False))
    save_status()
    try:
        if mask is None:
            tx, registration_qc = normalize_anatomy(mri, atlas, deformable=not affine_only)
            mask = sitk.Resample(template_brain_mask, mri, tx, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
            mask = _mask_on_grid(mask, mri)
            sitk.WriteTransform(tx, str(out/'mri_to_template_mask_pull.h5'))
            qc['brain_extraction'] = dict(method='Template brain-mask propagation', registration=registration_qc, review_required=True)
        else: qc['brain_extraction'] = dict(method='Supplied anatomical mask', review_required=True)
        values = sitk.GetArrayViewFromImage(mri)[sitk.GetArrayViewFromImage(mask) > 0]
        qc['n4_zero_voxels_excluded_from_estimation'] = int(np.count_nonzero(values == 0))
        corrected, bias = n4_correct(mri, mask, iterations=n4_iterations)
        brain = sitk.Mask(corrected, mask)
        for name, image in [('brain_mask', mask), ('T1w_bias_corrected', corrected), ('T1w_brain', brain), ('T1w_bias_field', bias)]:
            sitk.WriteImage(image, str(out/(name+'.nii.gz')))
        segmentation, probabilities, qc['tissue_segmentation'] = segment_t1(corrected, mask)
        sitk.WriteImage(segmentation, str(out/'tissue_classes.nii.gz'))
        for name, image in zip(('candidate_CSF', 'candidate_GM', 'candidate_WM'), probabilities):
            sitk.WriteImage(image, str(out/(name+'_probability.nii.gz')))
        qc['brain_volume_mL'] = float(sitk.GetArrayViewFromImage(mask).sum()*np.prod(mask.GetSpacing())/1000)
        result = preprocess_dynamic(pet, brain, checked['frames'].duration, reference_index,
                                    'identity' if preserve_blank else 'raise')
        qc['motion'] = result['motion_qc']; qc['pet_mri_registration'] = {k: result[k] for k in ('registration_metric', 'registration_stop', 'registration_candidates')}
        qc['final_interpolations_per_frame_per_output'] = 1
        sitk.WriteImage(result['image'], str(out/'pet_in_mri.nii.gz'))
        sitk.WriteImage(result['motion_corrected'], str(out/'motion_corrected.nii.gz'))
        for i, tx in enumerate(result['composed_transforms']):
            tx.FlattenTransform(); sitk.WriteTransform(tx, str(out/f'mri_to_original_frame_{i:03d}.h5'))
        for i, tx in enumerate(result['frame_transforms']): sitk.WriteTransform(tx, str(out/f'reference_to_frame_{i:03d}.tfm'))
        sitk.WriteTransform(result['mri_to_reference'], str(out/'mri_to_pet_reference.tfm'))
        with (out/'motion.csv').open('w', newline='') as f:
            fields = ['frame', 'mean_landmark_displacement_mm', 'maximum_landmark_displacement_mm', 'relative_to_previous_mm']
            writer = csv.DictWriter(f, fields, extrasaction='ignore'); writer.writeheader(); writer.writerows(qc['motion'])
        # Preserve the exact original metadata bytes; no inferred clock/decay correction.
        (out/'source_pet.json').write_bytes(Path(timing_path).read_bytes())
        if atlas_labels:
            tx, registration_qc = normalize_anatomy(brain, sitk.Mask(atlas, template_brain_mask) if template_brain_mask is not None else atlas,
                                                   deformable=not affine_only)
            region_labels = sitk.Resample(template_lab, mri, tx, sitk.sitkNearestNeighbor, 0, template_lab.GetPixelID())
            qc['atlas_propagation'] = registration_qc
            missing = sorted(set(np.unique(sitk.GetArrayFromImage(template_lab))) - set(np.unique(sitk.GetArrayFromImage(region_labels))))
            qc['atlas_labels_missing_after_resampling'] = [int(i) for i in missing if i > 0]
            sitk.WriteTransform(tx, str(out/'mri_to_atlas_pull.h5'))
        if region_labels is not None:
            _validate_labels(region_labels, mri)
            sitk.WriteImage(region_labels, str(out/'regions_in_mri.nii.gz'))
            from .workflow import extract_regions
            manifest = extract_regions(out/'pet_in_mri.nii.gz', out/'regions_in_mri.nii.gz', out/'source_pet.json', out/'regions')
            qc['regions'] = len(manifest['regions'])
            lab = sitk.GetArrayFromImage(region_labels); in_brain = sitk.GetArrayFromImage(mask) > 0
            qc['label_brain_overlap_fraction'] = float(np.mean(in_brain[lab > 0]))
            qc['regional_anatomical_qc'] = [dict(id=int(i), brain_overlap_fraction=float(np.mean(in_brain[lab == i])))
                                           for i in np.unique(lab) if i > 0]
            if qc['label_brain_overlap_fraction'] < .9:
                qc.setdefault('warnings', []).append('Less than 90% of labelled voxels lie inside the supplied/propagated brain mask; review the regional geometry')
            if pvc_fwhm is not None:
                pvc = gtm_correct(result['image'], region_labels, pvc_fwhm)
                folder = out/'regions_gtm'; folder.mkdir(); corrected_regions = []
                for i, key in enumerate(pvc['region_ids']):
                    if key == 0: continue
                    name = f'region_{key}.csv'
                    np.savetxt(folder/name, np.c_[checked['frames'].start, checked['frames'].end, pvc['corrected'][i]],
                               delimiter=',', header='start,end,activity', comments='')
                    corrected_regions.append(dict(id=key, name=f'region_{key}', path=name))
                gtm_manifest = dict(manifest, regions=corrected_regions, role='GTM-corrected regional estimates',
                                    preprocessing_record='../preprocessing.json')
                (folder/'regions.json').write_text(json.dumps(gtm_manifest, indent=2))
                np.savetxt(out/'gtm_matrix.csv', pvc['matrix'], delimiter=',')
                qc['gtm'] = {k: v for k, v in pvc.items() if k not in ('corrected', 'observed', 'matrix')}
        if normalize:
            target = sitk.Mask(atlas, template_brain_mask) if template_brain_mask is not None else atlas
            tx, qc['normalization'] = normalize_anatomy(target, brain, deformable=not affine_only)
            sitk.WriteTransform(tx, str(out/'template_to_mri_pull.h5'))
            pull = sitk.CompositeTransform([result['mri_to_reference'], tx])
            standard, _ = apply_composed_transforms(pet, atlas, result['frame_transforms'], pull)
            sitk.WriteImage(standard, str(out/'pet_in_template.nii.gz'))
            sitk.WriteImage(sitk.Resample(corrected, atlas, tx, sitk.sitkLinear, 0, sitk.sitkFloat32), str(out/'T1w_in_template.nii.gz'))
            if qc['normalization'].get('nonlinear_stop_requires_review'):
                qc.setdefault('warnings', []).append('Nonlinear normalization did not report a tolerance/convergence stop; inspect its optimizer status and overlays')
            _write_template_review(out, target)
        _write_review(out, brain, _duration_mean(result['image'], checked['frames'].duration), segmentation, qc, region_labels)
        qc['status'] = 'completed'; qc['review_required'] = True; save_status()
        return qc
    except Exception as exc:
        qc.update(status='failed', error=str(exc)); save_status(); raise


def _validate_labels(labels, reference):
    _check_image(labels, 'Region labels')
    if not _same_grid(labels, reference): raise ValueError('Region labels must match their anatomical image grid')
    a = sitk.GetArrayViewFromImage(labels)
    if np.any(a < 0) or not np.all(a == np.round(a)) or not np.any(a > 0):
        raise ValueError('Region labels must be nonnegative integers with at least one positive ID')


def _write_review(out, brain, pet_mean, tissues, qc, regions=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a = sitk.GetArrayFromImage(brain); p = sitk.GetArrayFromImage(pet_mean); t = sitk.GetArrayFromImage(tissues)
    rows = 3 if regions is not None else 2
    lab = sitk.GetArrayFromImage(regions) if regions is not None else None
    fig, axes = plt.subplots(rows, 3, figsize=(12, rows*3.5))
    for col, fraction in enumerate((.35, .5, .65)):
        z = int((len(a)-1)*fraction)
        axes[0, col].imshow(a[z], cmap='gray', origin='lower')
        axes[0, col].imshow(np.ma.masked_where(p[z] <= 0, p[z]), cmap='hot', alpha=.45, origin='lower')
        axes[1, col].imshow(a[z], cmap='gray', origin='lower')
        axes[1, col].imshow(np.ma.masked_where(t[z] == 0, t[z]), cmap='viridis', alpha=.5, vmin=1, vmax=3, origin='lower')
        if lab is not None:
            axes[2, col].imshow(a[z], cmap='gray', origin='lower')
            axes[2, col].imshow(np.ma.masked_where(lab[z] == 0, lab[z]), cmap='tab20', alpha=.55, origin='lower', interpolation='nearest')
        for ax in axes[:, col]: ax.set_axis_off(); ax.set_title(f'Array axial slice {z}')
    fig.suptitle('PET/MRI overlay; candidate tissue classes; supplied regions (if present). Inspect full volumes.')
    fig.tight_layout(); fig.savefig(out/'alignment_and_tissues.png', dpi=120); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot([x['frame'] for x in qc['motion']], [x['mean_landmark_displacement_mm'] for x in qc['motion']], marker='.')
    ax.set(xlabel='Frame', ylabel='Mean landmark displacement (mm)', title='Interframe motion estimate; not Power FD')
    fig.tight_layout(); fig.savefig(out/'motion.png', dpi=120); plt.close(fig)
    standard = '<p><a href="template_alignment.png">Template/MRI normalization overlay</a></p>' if (out/'template_alignment.png').exists() else ''
    (out/'index.html').write_text('<!doctype html><html><meta charset="utf-8"><title>PetKin preprocessing</title>'
        '<body><h1>Native PetKin preprocessing</h1><p>Processing completed. Anatomical review is required; '
        'this report does not certify accuracy.</p><img style="max-width:100%" src="alignment_and_tissues.png">'
        '<img style="max-width:100%" src="motion.png"><p><a href="preprocessing.json">Processing record</a> | '
        '<a href="motion.csv">Motion table</a></p>'+standard+'<ul>' + ''.join('<li>'+html.escape(s)+'</li>' for s in qc['limitations']+qc.get('warnings', [])) + '</ul></body></html>')


def _write_template_review(out, template):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fixed = sitk.GetArrayFromImage(template)
    moving = sitk.GetArrayFromImage(sitk.ReadImage(str(out/'T1w_in_template.nii.gz')))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, fraction in zip(axes, (.35, .5, .65)):
        z = int((len(fixed)-1)*fraction)
        ax.imshow(fixed[z], cmap='gray', origin='lower')
        ax.imshow(np.ma.masked_where(moving[z] <= 0, moving[z]), cmap='hot', origin='lower', alpha=.4)
        ax.set_axis_off(); ax.set_title(f'Array axial slice {z}')
    fig.suptitle('Template (gray) and normalized MRI (color); visual review required')
    fig.tight_layout(); fig.savefig(out/'template_alignment.png', dpi=120); plt.close(fig)
