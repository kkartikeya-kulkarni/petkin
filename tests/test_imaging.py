import numpy as np
import SimpleITK as sitk
import pytest
from petkin.imaging import register_rigid, motion_correct, roi_tacs


def phantom():
    z,y,x=np.mgrid[:40,:40,:40]
    a=np.exp(-((x-18)**2/60+(y-21)**2/100+(z-20)**2/80))
    a+=.7*np.exp(-((x-27)**2+(y-14)**2+(z-24)**2)/12)
    im=sitk.GetImageFromArray(a.astype('float32')); im.SetSpacing([2.,2.,2.])
    return im


@pytest.mark.parametrize('multimodal',[False,True])
def test_known_rigid(multimodal):
    fixed=phantom()
    true=sitk.Euler3DTransform(); true.SetCenter([40,40,40]); true.SetRotation(.025,-.035,.04); true.SetTranslation([3,-2,2])
    moving=sitk.Resample(fixed,fixed,true.GetInverse(),sitk.sitkLinear,0.)
    if multimodal: moving=sitk.Sqrt(moving)*4
    r=register_rigid(fixed,moving,multimodal)
    points=[[25,30,35],[40,40,40],[50,42,46]]
    tre=np.mean([np.linalg.norm(np.array(r.transform.TransformPoint(p))-true.TransformPoint(p)) for p in points])
    assert tre<1.0


def test_dynamic_and_labels():
    fixed=phantom(); tx=sitk.TranslationTransform(3,[3.,-2.,2.])
    moving=sitk.Resample(fixed,fixed,tx.GetInverse())
    dynamic=sitk.JoinSeries([fixed,moving]); dynamic.SetSpacing([2,2,2,30])
    out,transforms,qc=motion_correct(dynamic,0)
    assert out.GetSpacing()==dynamic.GetSpacing()
    a=sitk.GetArrayFromImage(out)
    assert np.mean((a[0]-a[1])**2) < np.mean((sitk.GetArrayFromImage(fixed)-sitk.GetArrayFromImage(moving))**2)*.1
    labels=sitk.Cast(fixed>.3,sitk.sitkUInt8)
    assert len(roi_tacs(out,labels)[1])==2
    labels.SetOrigin([10,0,0])
    with pytest.raises(ValueError): roi_tacs(out,labels)


def test_blank_frame_policy():
    fixed=phantom(); blank=fixed*0
    dynamic=sitk.JoinSeries([blank,fixed])
    with pytest.raises(ValueError,match="Blank frame"):
        motion_correct(dynamic,1)
    corrected,transforms,qc=motion_correct(dynamic,1,blank_policy="identity")
    assert not qc[0]['motion_estimated']
    assert np.all(sitk.GetArrayFromImage(corrected)[0]==0)
