import json
import numpy as np
import pytest
import SimpleITK as sitk
from scipy.ndimage import gaussian_filter

from petkin.preprocessing import (segment_t1, n4_correct, normalize_anatomy,
                                  gtm_correct, preprocess_native)


def anatomy():
    z, y, x = np.mgrid[:32, :32, :32]
    r = ((x-15)/12)**2 + ((y-16)/11)**2 + ((z-15)/10)**2
    classes = np.where(r < .2, 3, np.where(r < .6, 2, np.where(r < 1, 1, 0)))
    data = np.choose(classes, [0., 30., 65., 100.])
    data += np.random.default_rng(4).normal(0, .8, data.shape)*(classes > 0)
    def image(a):
        im = sitk.GetImageFromArray(a); im.SetSpacing([2., 2., 2.]); im.SetOrigin([-30., -32., -30.]); return im
    return image(data.astype(np.float32)), image((classes > 0).astype(np.uint8)), image(classes.astype(np.uint8))


def test_tissue_segmentation_recovers_distinct_classes():
    mri, mask, truth = anatomy()
    labels, probabilities, qc = segment_t1(mri, mask)
    inside = sitk.GetArrayViewFromImage(mask) > 0
    assert np.mean(sitk.GetArrayFromImage(labels)[inside] == sitk.GetArrayFromImage(truth)[inside]) > .99
    total = sum(sitk.GetArrayFromImage(p) for p in probabilities)
    assert np.allclose(total[inside], 1, atol=1e-6)
    assert np.all(total[~inside] == 0)
    assert qc['review_required']


def test_n4_improves_smooth_bias_and_preserves_source():
    mri, mask, truth = anatomy()
    data = sitk.GetArrayFromImage(mri); original = data.copy()
    x = np.indices(data.shape)[2]
    bias = np.exp(.6*(x-16)/32)
    biased = sitk.GetImageFromArray((data*bias).astype(np.float32)); biased.CopyInformation(mri)
    corrected, field = n4_correct(biased, mask, iterations=(20, 10))
    within = sitk.GetArrayFromImage(truth) == 2
    before = sitk.GetArrayFromImage(biased)[within]; after = sitk.GetArrayFromImage(corrected)[within]
    assert after.std()/after.mean() < before.std()/before.mean()
    assert np.array_equal(sitk.GetArrayFromImage(mri), original)
    assert np.all(sitk.GetArrayFromImage(field) > 0)


def test_gtm_recovers_known_blurred_region_activities():
    _, _, labels = anatomy(); lab = sitk.GetArrayFromImage(labels)
    fwhm = np.array([5., 6., 7.]); sigma = (fwhm/np.sqrt(8*np.log(2))/2)[::-1]
    truth = np.array([[2., 3.], [8., 10.], [20., 25.], [40., 35.]])
    images = []
    for t in range(2):
        activity = truth[:, t][lab]
        blurred = gaussian_filter(activity, sigma, mode='constant', truncate=4.)
        image = sitk.GetImageFromArray(blurred.astype(np.float32)); image.CopyInformation(labels); images.append(image)
    result = gtm_correct(sitk.JoinSeries(images), labels, fwhm)
    assert result['region_ids'] == [0, 1, 2, 3]
    assert np.max(np.abs(result['corrected']-truth)) < 1e-4
    assert np.linalg.norm(result['corrected']-truth) < np.linalg.norm(result['observed']-truth)*.01
    with pytest.raises(ValueError, match='FWHM'): gtm_correct(sitk.JoinSeries(images), labels, [-1, 2, 3])


def test_affine_normalization_physical_translation():
    mri, _, _ = anatomy()
    translated = sitk.Image(mri); translated.SetOrigin(np.array(mri.GetOrigin())+[4., -2., 3.])
    tx, qc = normalize_anatomy(mri, translated, deformable=False, iterations=25)
    p = [0., 0., 0.]
    assert np.linalg.norm(np.array(tx.TransformPoint(p))-[4., -2., 3.]) < .5
    assert qc['minimum_jacobian'] > 0


def test_deformable_normalization_identity():
    mri, _, _ = anatomy()
    tx, qc = normalize_anatomy(mri, mri, iterations=5)
    assert np.linalg.norm(tx.TransformPoint([0., 0., 0.])) < 1.
    assert qc['minimum_jacobian'] > 0


def sources(tmp_path):
    mri, mask, labels = anatomy()
    # PET has a different contrast, with an identifiable asymmetric focus.
    pet = sitk.Sqrt(mri)*2
    a = sitk.GetArrayFromImage(pet); a[14:18, 17:21, 11:15] += 10
    pet = sitk.GetImageFromArray(a); pet.CopyInformation(mri)
    dynamic = sitk.JoinSeries([pet, pet*1.1])
    for name, im in [('pet', dynamic), ('mri', mri), ('mask', mask), ('labels', labels)]:
        sitk.WriteImage(im, str(tmp_path/(name+'.nii.gz')))
    metadata = dict(FrameTimesStart=[0, 30], FrameDuration=[30, 30], Units='Bq/mL',
                    TimeZero='12:00:00', InjectionStart=0, ImageDecayCorrected=True, ImageDecayCorrectionTime=0)
    (tmp_path/'pet.json').write_text(json.dumps(metadata))
    return tmp_path/'pet.nii.gz', tmp_path/'mri.nii.gz', tmp_path/'pet.json'


def test_full_native_pipeline_and_original_metadata(tmp_path):
    paths = sources(tmp_path); original = paths[2].read_bytes(); out = tmp_path/'out'
    result = preprocess_native(*paths, out, brain_mask=tmp_path/'mask.nii.gz',
                               labels=tmp_path/'labels.nii.gz', pvc_fwhm=np.array([5, 5, 5]), n4_iterations=(5,))
    assert result['status'] == 'completed' and not result['fully_validated']
    assert result['regions'] == 3 and result['final_interpolations_per_frame_per_output'] == 1
    assert (out/'source_pet.json').read_bytes() == original == paths[2].read_bytes()
    manifest = json.loads((out/'regions_gtm/regions.json').read_text())
    assert len(manifest['regions']) == 3
    assert (out/'index.html').is_file() and (out/'motion.csv').is_file()
    with pytest.raises(ValueError, match='empty'): preprocess_native(*paths, out, brain_mask=tmp_path/'mask.nii.gz')


def test_invalid_options_fail_without_creating_output(tmp_path):
    paths = sources(tmp_path); out = tmp_path/'out'
    with pytest.raises(ValueError, match='brain-mask'): preprocess_native(*paths, out)
    with pytest.raises(ValueError, match='PVC requires'):
        preprocess_native(*paths, out, brain_mask=tmp_path/'mask.nii.gz', pvc_fwhm=[5, 5, 5])
    assert not out.exists()


def test_failed_processing_is_recorded(tmp_path, monkeypatch):
    paths = sources(tmp_path); out = tmp_path/'out'
    def fail(*args, **kwargs): raise RuntimeError('unmeasurable motion')
    monkeypatch.setattr('petkin.preprocessing.preprocess_dynamic', fail)
    with pytest.raises(RuntimeError, match='unmeasurable'):
        preprocess_native(*paths, out, brain_mask=tmp_path/'mask.nii.gz', n4_iterations=(2,))
    record = json.loads((out/'preprocessing.json').read_text())
    assert record['status'] == 'failed' and 'unmeasurable' in record['error']


def test_template_mask_atlas_and_direct_standard_output(tmp_path):
    paths = sources(tmp_path); out = tmp_path/'standard'
    result = preprocess_native(*paths, out, template=tmp_path/'mri.nii.gz',
                               template_mask=tmp_path/'mask.nii.gz', atlas_labels=tmp_path/'labels.nii.gz',
                               normalize=True, affine_only=True, n4_iterations=(3,))
    assert result['status'] == 'completed' and result['regions'] == 3
    assert result['brain_extraction']['method'] == 'Template brain-mask propagation'
    assert (out/'template_to_mri_pull.h5').is_file()
    standard = sitk.ReadImage(str(out/'pet_in_template.nii.gz'))
    source = sitk.ReadImage(str(paths[0]))
    assert standard.GetSize() == source.GetSize()
    assert result['normalization']['minimum_jacobian'] > 0


def test_n4_excludes_isolated_zeros_without_changing_mask():
    mri, mask, _ = anatomy(); a = sitk.GetArrayFromImage(mri); a[15, 16, 15] = 0
    image = sitk.GetImageFromArray(a); image.CopyInformation(mri)
    mask_before = sitk.GetArrayFromImage(mask)
    corrected, _ = n4_correct(image, mask, iterations=(3,))
    assert sitk.GetArrayFromImage(corrected)[15, 16, 15] == 0
    assert np.array_equal(mask_before, sitk.GetArrayFromImage(mask))


def test_gtm_condition_rejection():
    mri, _, labels = anatomy()
    with pytest.raises(ValueError, match='ill-conditioned'):
        gtm_correct(sitk.JoinSeries([mri]), labels, [5, 5, 5], maximum_condition=1)
