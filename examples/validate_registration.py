"""Physical-space registration accuracy against known synthetic rigid motion."""
import json,sys
from pathlib import Path
import numpy as np
import SimpleITK as sitk
from petkin.imaging import register_rigid,pet_to_mri


def run(destination):
    z,y,x=np.mgrid[:40,:40,:40]
    a=np.exp(-((x-18)**2/60+(y-21)**2/100+(z-20)**2/80))
    a+=.7*np.exp(-((x-27)**2+(y-14)**2+(z-24)**2)/12)
    fixed=sitk.GetImageFromArray(a.astype('float32')); fixed.SetSpacing([2.,2.,2.])
    true=sitk.Euler3DTransform(); true.SetCenter([40,40,40]); true.SetRotation(.025,-.035,.04); true.SetTranslation([3,-2,2])
    moving=sitk.Resample(fixed,fixed,true.GetInverse(),sitk.sitkLinear,0.)
    points=[[25,30,35],[40,40,40],[50,42,46]]; results={}
    for multimodal in [False,True]:
        im=sitk.Sqrt(moving)*4 if multimodal else moving
        r=register_rigid(fixed,im,multimodal)
        error=[float(np.linalg.norm(np.array(r.transform.TransformPoint(p))-true.TransformPoint(p))) for p in points]
        results['PET_MRI_contrast_surrogate' if multimodal else 'PET_motion']=dict(landmark_errors_mm=error,mean_error_mm=float(np.mean(error)),passed=bool(max(error)<1),stop=r.stop)
    # Exercise complete mean-PET-to-MRI workflow with the known transform.
    dynamic=sitk.JoinSeries([moving,moving*2])
    out,result=pet_to_mri(dynamic,fixed,[1.,2.])
    results['dynamic_pipeline']=dict(output_size=list(out.GetSize()),fixed_size=list(fixed.GetSize()),passed=out.GetSize()[:3]==fixed.GetSize())
    results['scope']='Asymmetric Gaussian phantom; contrast surrogate is not biological MRI validation'
    Path(destination).write_text(json.dumps(results,indent=2)); print(json.dumps(results,indent=2))
    assert all(r['passed'] for r in results.values() if isinstance(r,dict))
if __name__=='__main__': run(sys.argv[1])
