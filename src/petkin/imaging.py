"""Physical-space rigid registration. SimpleITK uses LPS coordinates in mm."""
from dataclasses import dataclass
import numpy as np
import SimpleITK as sitk


@dataclass
class Registration:
    image: sitk.Image
    transform: sitk.Transform  # fixed/reference physical points -> moving points
    metric: float
    stop: str
    candidates: list | None = None


def _register_rigid_single(fixed, moving, multimodal=True, iterations=200, initialization="geometry", fixed_mask=None, moving_mask=None):
    for image in (fixed,moving):
        if image.GetDimension() != 3: raise ValueError("Registration requires 3D images")
        a = sitk.GetArrayViewFromImage(image)
        if not np.all(np.isfinite(a)) or np.std(a) < 1e-8: raise ValueError("Nonfinite or constant image")
    fixed, moving = sitk.Cast(fixed,sitk.sitkFloat32), sitk.Cast(moving,sitk.sitkFloat32)
    init = sitk.CenteredTransformInitializer(fixed,moving,sitk.Euler3DTransform(),
                                           sitk.CenteredTransformInitializerFilter.GEOMETRY if initialization=="geometry" else sitk.CenteredTransformInitializerFilter.MOMENTS)
    reg = sitk.ImageRegistrationMethod()
    reg.SetNumberOfThreads(1)
    for mask,image,setter in [(fixed_mask,fixed,reg.SetMetricFixedMask),(moving_mask,moving,reg.SetMetricMovingMask)]:
        if mask is not None:
            if mask.GetDimension()!=3 or mask.GetSize()!=image.GetSize() or any(
                not np.allclose(getattr(mask,k)(),getattr(image,k)()) for k in ("GetOrigin","GetSpacing","GetDirection")):
                raise ValueError("Registration mask must match its image geometry")
            if not np.any(sitk.GetArrayViewFromImage(mask)>0):raise ValueError("Empty registration mask")
            setter(sitk.Cast(mask>0,sitk.sitkUInt8))
    if multimodal: reg.SetMetricAsMattesMutualInformation(40)
    else: reg.SetMetricAsCorrelation()
    reg.SetMetricSamplingStrategy(reg.NONE)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetShrinkFactorsPerLevel([4,2,1]); reg.SetSmoothingSigmasPerLevel([2,1,0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    reg.SetOptimizerAsRegularStepGradientDescent(learningRate=1.,minStep=.001,
                                               numberOfIterations=iterations,gradientMagnitudeTolerance=1e-7)
    reg.SetOptimizerScalesFromPhysicalShift(); reg.SetInitialTransform(init,inPlace=False)
    transform = reg.Execute(fixed,moving)
    image = sitk.Resample(moving,fixed,transform,sitk.sitkLinear,0.,sitk.sitkFloat32)
    return Registration(image,transform,float(reg.GetMetricValue()),reg.GetOptimizerStopConditionDescription())


def register_rigid(fixed,moving,multimodal=True,iterations=200,*,initializations=None,fixed_mask=None,moving_mask=None):
    """Deterministic rigid multi-initialization; metric choice is not accuracy QC."""
    choices=initializations or (("geometry","moments") if multimodal else ("geometry",))
    if not choices or any(c not in ("geometry","moments") for c in choices):raise ValueError("Unknown registration initialization")
    candidates=[];results=[]
    for choice in choices:
        try:
            result=_register_rigid_single(fixed,moving,multimodal,iterations,choice,fixed_mask,moving_mask)
            if not np.isfinite(result.metric):raise RuntimeError("Nonfinite registration metric")
            candidates.append(dict(initialization=choice,metric=result.metric,stop=result.stop,status="completed"))
            results.append(result)
        except RuntimeError as error:
            candidates.append(dict(initialization=choice,status="failed",reason=str(error)))
    if not results:raise RuntimeError("All rigid registration initializations failed: "+str(candidates))
    best=min(results,key=lambda r:r.metric);best.candidates=candidates
    return best


def split_dynamic(dynamic):
    if dynamic.GetDimension() != 4: raise ValueError("Expected scalar 4D PET (x,y,z,time)")
    # Reject direction matrices that couple time to spatial coordinates.
    d = np.array(dynamic.GetDirection()).reshape(4,4)
    if not np.allclose(d[:3,3],0) or not np.allclose(d[3,:3],0):
        raise ValueError("Time-space direction coupling is unsupported")
    size = list(dynamic.GetSize()); size[3] = 0
    return [sitk.Extract(dynamic,size,[0,0,0,i]) for i in range(dynamic.GetSize()[3])]


def motion_correct(dynamic, reference_index=None, blank_policy="raise"):
    """Register each frame independently to a selected PET frame.

    Default: frame with highest total activity. Early low-count frames can fail;
    inspect returned metrics/transforms and overlay images before interpretation.
    """
    if blank_policy not in ("raise", "identity"): raise ValueError("Invalid blank policy")
    frames = split_dynamic(dynamic)
    if reference_index is None:
        reference_index = int(np.argmax([np.sum(sitk.GetArrayViewFromImage(f)) for f in frames]))
    if not 0 <= reference_index < len(frames): raise ValueError("Invalid reference frame")
    fixed = frames[reference_index]; corrected=[]; transforms=[]; qc=[]
    if np.std(sitk.GetArrayViewFromImage(fixed)) < 1e-8: raise ValueError("Blank reference frame")
    for i,frame in enumerate(frames):
        values = sitk.GetArrayViewFromImage(frame)
        if not np.all(np.isfinite(values)): raise ValueError("Nonfinite PET frame")
        if np.std(values) < 1e-8:
            if blank_policy == "raise": raise ValueError(f"Blank frame {i}: motion is unmeasurable")
            corrected.append(frame); transforms.append(sitk.Transform(3,sitk.sitkIdentity))
            qc.append(dict(frame=i,status="blank_preserved",motion_estimated=False))
        elif i == reference_index:
            corrected.append(frame); transforms.append(sitk.Transform(3,sitk.sitkIdentity))
            qc.append(dict(frame=i,reference=True))
        else:
            r = register_rigid(fixed,frame,multimodal=False)
            corrected.append(r.image); transforms.append(r.transform)
            qc.append(dict(frame=i,metric=r.metric,stop=r.stop))
    out = sitk.JoinSeries(corrected,dynamic.GetOrigin()[3],dynamic.GetSpacing()[3])
    return out,transforms,qc


def pet_to_mri(corrected_dynamic, mri, durations):
    """Register duration-weighted mean PET to MRI; resample every PET frame."""
    frames = split_dynamic(corrected_dynamic)
    w = np.asarray(durations,float)
    if w.shape != (len(frames),) or not np.all(np.isfinite(w)) or np.any(w<=0):
        raise ValueError("Positive frame durations required")
    mean = sitk.Cast(frames[0],sitk.sitkFloat32)*float(w[0]/w.sum())
    for f,v in zip(frames[1:],w[1:]): mean += sitk.Cast(f,sitk.sitkFloat32)*float(v/w.sum())
    result = register_rigid(mri,mean,True)
    resampled = [sitk.Resample(f,mri,result.transform,sitk.sitkLinear,0.,sitk.sitkFloat32) for f in frames]
    out = sitk.JoinSeries(resampled,corrected_dynamic.GetOrigin()[3],corrected_dynamic.GetSpacing()[3])
    return out,result


def roi_tacs(dynamic, labels):
    frames = split_dynamic(dynamic)
    ref = frames[0]
    if labels.GetDimension()!=3 or labels.GetSize()!=ref.GetSize() or any(
        not np.allclose(getattr(labels,k)(),getattr(ref,k)()) for k in ("GetSpacing","GetOrigin","GetDirection")):
        raise ValueError("Labels must be on the PET spatial grid; resample with nearest neighbor first")
    lab = sitk.GetArrayFromImage(labels)
    if not np.all(np.isfinite(lab)) or np.any(lab<0) or not np.all(lab==np.round(lab)):
        raise ValueError("Labels must be nonnegative integers")
    return {int(k):np.array([sitk.GetArrayViewFromImage(f)[lab==k].mean() for f in frames])
            for k in np.unique(lab) if k>0}


def _duration_mean(dynamic,durations):
    frames=split_dynamic(dynamic);w=np.asarray(durations,float)
    if w.shape!=(len(frames),) or not np.all(np.isfinite(w)) or np.any(w<=0):
        raise ValueError('One finite positive duration per frame required')
    result=sitk.Cast(frames[0],sitk.sitkFloat32)*float(w[0]/w.sum())
    for frame,weight in zip(frames[1:],w[1:]):result+=sitk.Cast(frame,sitk.sitkFloat32)*float(weight/w.sum())
    return result


def apply_composed_transforms(dynamic,mri,reference_to_frames,mri_to_reference):
    """Final PET output: one interpolation from each original frame onto MRI.

    Transform order is T_frame(T_mri(point)). SimpleITK evaluates the last added
    transform first. External rigid transforms (e.g. PETPrep outputs) may be used
    only after explicitly converting them to this LPS fixed-to-moving convention.
    """
    frames=split_dynamic(dynamic)
    if mri.GetDimension()!=3 or len(reference_to_frames)!=len(frames):raise ValueError('Transform count or MRI dimension mismatch')
    if mri_to_reference.GetDimension()!=3 or any(t.GetDimension()!=3 for t in reference_to_frames):
        raise ValueError('All transforms must be three-dimensional')
    resampled=[];composed=[]
    for frame,tx in zip(frames,reference_to_frames):
        transform=sitk.CompositeTransform([tx,mri_to_reference])
        resampled.append(sitk.Resample(frame,mri,transform,sitk.sitkLinear,0.,sitk.sitkFloat32))
        composed.append(transform)
    return sitk.JoinSeries(resampled,dynamic.GetOrigin()[3],dynamic.GetSpacing()[3]),composed


def motion_displacements(reference,transforms,radius_mm=50.):
    """Physical displacement QC on seven landmarks; not clinical motion accuracy.

    This is mean landmark displacement, not the Power framewise displacement
    statistic. Intraframe motion remains unobservable from reconstructed frames.
    """
    if not np.isfinite(radius_mm) or radius_mm<=0:raise ValueError('Positive finite landmark radius required')
    center=np.array(reference.TransformContinuousIndexToPhysicalPoint(((np.array(reference.GetSize())-1)/2).tolist()))
    points=[center]+[center+sign*radius_mm*np.eye(3)[axis] for axis in range(3) for sign in [-1,1]]
    previous=None;rows=[]
    for i,transform in enumerate(transforms):
        moved=np.array([transform.TransformPoint(p.tolist()) for p in points])
        absolute=np.linalg.norm(moved-np.array(points),axis=1)
        relative=None if previous is None else float(np.mean(np.linalg.norm(moved-previous,axis=1)))
        rows.append(dict(frame=i,mean_landmark_displacement_mm=float(absolute.mean()),
                         maximum_landmark_displacement_mm=float(absolute.max()),relative_to_previous_mm=relative))
        previous=moved
    return rows


def preprocess_dynamic(dynamic,mri,durations,reference_index=None,blank_policy='raise'):
    """Estimate on temporary corrected images; resample original PET once at end."""
    corrected,frame_transforms,qc=motion_correct(dynamic,reference_index,blank_policy)
    reference=_duration_mean(corrected,durations)
    registration=register_rigid(mri,reference,multimodal=True)
    final,composed=apply_composed_transforms(dynamic,mri,frame_transforms,registration.transform)
    landmarks=motion_displacements(split_dynamic(dynamic)[0],frame_transforms)
    for index,(row,displacement) in enumerate(zip(qc,landmarks)):
        row.update(displacement)
        if row.get('motion_estimated') is False:
            row['mean_landmark_displacement_mm']=None
            row['maximum_landmark_displacement_mm']=None
            row['relative_to_previous_mm']=None
        if index>0 and qc[index-1].get('motion_estimated') is False: row['relative_to_previous_mm']=None
    return dict(image=final,motion_corrected=corrected,frame_transforms=frame_transforms,
                mri_to_reference=registration.transform,composed_transforms=composed,motion_qc=qc,
                registration_metric=registration.metric,registration_stop=registration.stop,registration_candidates=registration.candidates,
                final_interpolations_per_frame=1,
                limitations=['Rigid interframe motion only','No intraframe or attenuation mismatch correction',
                             'Human alignment requires independent landmarks/segmentation review'])


def landmark_registration_error(fixed_points,moving_points,transform):
    """Measure TRE using independent corresponding physical LPS landmarks.

    This function does not create anatomical truth or invent correspondences.
    The supplied transform maps fixed points to the moving image coordinates.
    """
    fixed=np.asarray(fixed_points,float);moving=np.asarray(moving_points,float)
    if fixed.ndim!=2 or fixed.shape[1]!=3 or fixed.shape!=moving.shape or len(fixed)<3 or not np.isfinite(fixed).all() or not np.isfinite(moving).all():
        raise ValueError('At least three finite corresponding 3D physical landmarks required')
    if transform.GetDimension()!=3:raise ValueError('Three-dimensional transform required')
    mapped=np.array([transform.TransformPoint(p.tolist()) for p in fixed])
    error=np.linalg.norm(mapped-moving,axis=1);before=np.linalg.norm(fixed-moving,axis=1)
    return dict(landmarks=len(fixed),errors_mm=error.tolist(),mean_error_mm=float(error.mean()),
                rms_error_mm=float(np.sqrt(np.mean(error**2))),maximum_error_mm=float(error.max()),
                unregistered_rms_error_mm=float(np.sqrt(np.mean(before**2))),
                convention='Fixed-to-moving physical LPS millimeters',
                limitation='Accuracy depends on independent anatomical landmark quality')
