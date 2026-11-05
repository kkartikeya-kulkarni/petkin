"""Reproduce synthetic native preprocessing checks; no human accuracy claim.

Usage: python examples/validate_native_preprocessing.py EMPTY_OUTPUT_DIRECTORY
"""
import json
from pathlib import Path
import sys
import numpy as np
import SimpleITK as sitk
from scipy.ndimage import gaussian_filter
from petkin.preprocessing import segment_t1, n4_correct, gtm_correct, normalize_anatomy, preprocess_native


def main(output):
    out = Path(output)
    if out.exists() and any(out.iterdir()): raise ValueError('Output directory must be empty')
    out.mkdir(parents=True, exist_ok=True)
    z,y,x = np.mgrid[:40,:40,:40]
    radius = ((x-19)/15)**2 + ((y-20)/14)**2 + ((z-19)/13)**2
    truth = np.where(radius<.2,3,np.where(radius<.6,2,np.where(radius<1,1,0)))
    rng = np.random.default_rng(620)
    data = np.choose(truth,[0.,30.,65.,100.]) + rng.normal(0,.8,truth.shape)*(truth>0)
    def image(a):
        im=sitk.GetImageFromArray(a);im.SetSpacing([2.,2.,2.]);im.SetOrigin([-38.,-40.,-38.]);return im
    mri=image(data.astype(np.float32));mask=image((truth>0).astype(np.uint8));labels=image(truth.astype(np.uint8))
    classes,_,segqc=segment_t1(mri,mask)
    accuracy=float(np.mean(sitk.GetArrayFromImage(classes)[truth>0]==truth[truth>0]))
    biased=image((data*np.exp(.6*(x-20)/40)).astype(np.float32))
    corrected,_=n4_correct(biased,mask,iterations=(20,10))
    within=truth==2;before=sitk.GetArrayFromImage(biased)[within];after=sitk.GetArrayFromImage(corrected)[within]
    cv_before=float(before.std()/before.mean());cv_after=float(after.std()/after.mean())
    fwhm=np.array([5.,6.,7.]);sigma=(fwhm/np.sqrt(8*np.log(2))/2)[::-1]
    activities=np.array([[2.,3.,2.5],[8.,10.,9.],[20.,25.,24.],[40.,35.,32.]])
    frames=[image(gaussian_filter(activities[:,t][truth],sigma,mode='constant').astype(np.float32)) for t in range(3)]
    gtm=gtm_correct(sitk.JoinSeries(frames),labels,fwhm)
    shifted=sitk.Image(mri);shifted.SetOrigin(np.array(mri.GetOrigin())+[4.,-2.,3.])
    affine,affine_qc=normalize_anatomy(mri,shifted,deformable=False)
    points=np.array([[0.,0.,0.],[-10.,4.,8.],[12.,-6.,-8.]])
    tre=np.linalg.norm(np.array([affine.TransformPoint(p.tolist()) for p in points])-(points+[4.,-2.,3.]),axis=1)
    # Anatomical contrast is intentionally different from PET.
    pet_array=np.sqrt(np.maximum(data,0))*2;pet_array[17:22,21:25,14:18]+=10
    pet=image(pet_array.astype(np.float32));pet_frames=[pet];poses=[sitk.Transform(3,sitk.sitkIdentity)]
    for offset in ([3.,-2.,1.],[-1.,2.,-2.]):
        pose=sitk.TranslationTransform(3,offset);poses.append(pose)
        pet_frames.append(sitk.Resample(pet,pet,pose.GetInverse(),sitk.sitkLinear,0,sitk.sitkFloat32))
    dynamic=sitk.JoinSeries(pet_frames)
    for name,im in [('pet',dynamic),('mri',mri),('mask',mask),('labels',labels)]:sitk.WriteImage(im,str(out/(name+'.nii.gz')))
    meta=dict(FrameTimesStart=[0,30,60],FrameDuration=[30,30,30],Units='Bq/mL',
              TimeZero='12:00:00',InjectionStart=0,ImageDecayCorrected=True,ImageDecayCorrectionTime=0)
    (out/'pet.json').write_text(json.dumps(meta,indent=2))
    result=preprocess_native(out/'pet.nii.gz',out/'mri.nii.gz',out/'pet.json',out/'native',
                             brain_mask=out/'mask.nii.gz',labels=out/'labels.nii.gz',reference_index=0,n4_iterations=(10,5))
    motion_errors=[]
    for i,pose in enumerate(poses):
        estimated=sitk.ReadTransform(str(out/'native'/f'reference_to_frame_{i:03d}.tfm'))
        errors=np.linalg.norm(np.array([estimated.TransformPoint(p.tolist()) for p in points])-
                              np.array([pose.TransformPoint(p.tolist()) for p in points]),axis=1)
        motion_errors.extend(errors.tolist())
    report=dict(scope='Synthetic development checks only; not independent human validation',
                seed=620,tissue_voxel_accuracy=accuracy,
                n4_within_tissue_cv_before=cv_before,n4_within_tissue_cv_after=cv_after,
                gtm_maximum_absolute_error=float(np.max(np.abs(gtm['corrected']-activities))),
                gtm_uncorrected_maximum_absolute_error=float(np.max(np.abs(gtm['observed']-activities))),
                affine_landmark_maximum_error_mm=float(tre.max()),
                rigid_motion_landmark_maximum_error_mm=float(max(motion_errors)),
                native_status=result['status'],registration_note='Motion and affine tests use known synthetic transforms; not measured human truth')
    (out/'summary.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if accuracy<.99 or cv_after>=cv_before or report['gtm_maximum_absolute_error']>1e-3 or tre.max()>1 or max(motion_errors)>1:
        raise RuntimeError('Synthetic preprocessing acceptance check failed; inspect summary')


if __name__=='__main__':main(sys.argv[1])
