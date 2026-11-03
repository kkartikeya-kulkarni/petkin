"""Known-transform verification using MRI-derived synthetic PET; not human PET truth."""
from pathlib import Path
import json,sys
import numpy as np
import SimpleITK as sitk
from petkin.imaging import register_rigid,landmark_registration_error


def run(root,out):
    root,out=Path(root),Path(out);protocol=json.loads((out/'PROTOCOL.json').read_text())
    rng=np.random.default_rng(protocol['holdout']['registration_seed']);records=[]
    for subject in ['s028','s029','s030']:
        mri=sitk.Cast(sitk.ReadImage(str(root/subject/'mri_4mm.nii.gz')),sitk.sitkFloat32)
        smooth=sitk.SmoothingRecursiveGaussian(sitk.Maximum(mri,0),2.5)
        pet=sitk.Pow(sitk.Maximum(smooth,0),.6);base=sitk.GetArrayFromImage(pet)
        noisy=rng.poisson(base/max(base.max(),1)*10000)*base.max()/10000
        pet=sitk.GetImageFromArray(noisy.astype('float32'));pet.CopyInformation(mri)
        center=np.array(mri.TransformContinuousIndexToPhysicalPoint(((np.array(mri.GetSize())-1)/2).tolist()))
        points=np.array([center]+[center+s*30*np.eye(3)[a] for a in range(3) for s in [-1,1]])
        for pose in range(protocol['holdout']['registration_poses_per_mri']):
            truth=sitk.Euler3DTransform();truth.SetCenter(center.tolist())
            angles=rng.uniform(-5,5,3)*np.pi/180;translation=rng.uniform(-8,8,3)
            truth.SetRotation(*map(float,angles));truth.SetTranslation(translation.tolist())
            moving=sitk.Resample(pet,pet,truth.GetInverse(),sitk.sitkLinear,0.,sitk.sitkFloat32)
            fit=register_rigid(mri,moving,multimodal=True)
            corresponding=np.array([truth.TransformPoint(p.tolist()) for p in points])
            errors=landmark_registration_error(points,corresponding,fit.transform)
            records.append(dict(subject=subject,pose=pose,rotation_radians=angles.tolist(),translation_mm=translation.tolist(),
                                errors=errors,candidates=fit.candidates,metric=fit.metric))
            print(subject,pose,errors['mean_error_mm'],errors['maximum_error_mm'],flush=True)
    summary=dict(cases=len(records),landmarks_per_case=7,mean_landmark_error_mm=float(np.mean([r['errors']['mean_error_mm'] for r in records])),
                 maximum_landmark_error_mm=max(r['errors']['maximum_error_mm'] for r in records),
                 measured_human_PET_truth=False,
                 scope='Intensity-transformed, blurred MRI-derived synthetic PET with Poisson noise and known rigid poses; not biological PET-MRI accuracy')
    (out/'registration_holdout.json').write_text(json.dumps(dict(summary=summary,cases=records),indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
