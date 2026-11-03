import json
from pathlib import Path
import numpy as np
import pytest
import SimpleITK as sitk
from scipy.integrate import solve_ivp
from petkin.workflow import extract_regions,analyze_regions


def observations(tmp_path):
    edges=np.arange(0,61,1.);time=np.linspace(0,60,6001);plasma=8000*time*np.exp(-time/.7)+1000*(1-np.exp(-time/.1))*np.exp(-time/25)
    def rhs(t,y):return [.2*np.interp(t,time,plasma)-.15*y[0],y[0]]
    solution=solve_ivp(rhs,[0,60],[0,0],rtol=1e-9,atol=1e-9,max_step=.03,dense_output=True)
    tac=np.diff(solution.sol(edges)[1])/np.diff(edges)
    np.savetxt(tmp_path/'input.csv',np.c_[time,plasma],delimiter=',',header='time,plasma',comments='')
    return edges,tac


def test_image_to_regions_to_analysis_recovers_independent_known_VT(tmp_path):
    edges,tac=observations(tmp_path);labels=np.ones((4,4,4),np.int16);labels[:2]=2
    frames=[]
    for value in tac:
        image=sitk.GetImageFromArray(np.where(labels==1,value,value*.5).astype('float32'));image.SetSpacing([2,2,2]);frames.append(image)
    pet=sitk.JoinSeries(frames);lab=sitk.GetImageFromArray(labels);lab.CopyInformation(frames[0]);sitk.WriteImage(pet,str(tmp_path/'pet.nii.gz'));sitk.WriteImage(lab,str(tmp_path/'labels.nii.gz'))
    timing=dict(FrameTimesStart=(edges[:-1]*60).tolist(),FrameDuration=[60]*60,Units='Bq/mL',ImageDecayCorrected=True,ImageDecayCorrectionTime=0)
    (tmp_path/'timing.json').write_text(json.dumps(timing))
    manifest=extract_regions(tmp_path/'pet.nii.gz',tmp_path/'labels.nii.gz',tmp_path/'timing.json',tmp_path/'regions')
    assert len(manifest['regions'])==2 and all(x['voxels']==32 for x in manifest['regions'])
    result=analyze_regions(tmp_path/'regions/regions.json',tmp_path/'input.csv',tmp_path/'analysis',model='1t',n_starts=2)
    assert result['completed']==2 and result['local_qc_passed']==2
    assert result['records'][0]['estimate']==pytest.approx(.2/.15,rel=.01)
    assert result['records'][1]['estimate']==pytest.approx(.1/.15,rel=.01)
    assert (tmp_path/'analysis/index.html').exists() and (tmp_path/'analysis/results.csv').exists()
    assert result['records'][0]['held_frame_prediction']['normalized_RMSE']<.003
    with pytest.raises(ValueError,match='must be empty'):analyze_regions(tmp_path/'regions/regions.json',tmp_path/'input.csv',tmp_path/'analysis')


def test_one_bad_region_preserved_without_losing_good_result(tmp_path):
    edges,tac=observations(tmp_path);np.savetxt(tmp_path/'good.csv',np.c_[edges[:-1],edges[1:],tac],delimiter=',',header='start,end,activity',comments='')
    (tmp_path/'regions.json').write_text(json.dumps({'regions':[{'id':1,'name':'Cortex <left>','path':'good.csv'},{'id':2,'name':'Missing','path':'absent.csv'}]}))
    result=analyze_regions(tmp_path/'regions.json',tmp_path/'input.csv',tmp_path/'out',model='1t',n_starts=2)
    assert result['completed']==1 and result['failed']==1 and len(result['records'])==2
    assert 'Cortex &lt;left&gt;' in (tmp_path/'out/index.html').read_text()
    assert result['records'][1]['status']=='failed'


def test_rejects_geometry_mismatch_before_writing_TACS(tmp_path):
    pet=sitk.JoinSeries([sitk.Image([4,4,4],sitk.sitkFloat32)]*5);labels=sitk.Image([4,4,4],sitk.sitkUInt8)+1;labels.SetOrigin([10,0,0]);sitk.WriteImage(pet,str(tmp_path/'pet.nii.gz'));sitk.WriteImage(labels,str(tmp_path/'labels.nii.gz'))
    (tmp_path/'timing.json').write_text(json.dumps(dict(FrameTimesStart=[0,60,120,180,240],FrameDuration=[60]*5,Units='Bq/mL',ImageDecayCorrected=True,ImageDecayCorrectionTime=0)))
    with pytest.raises(ValueError,match='grid'):extract_regions(tmp_path/'pet.nii.gz',tmp_path/'labels.nii.gz',tmp_path/'timing.json',tmp_path/'regions')
    assert not (tmp_path/'regions').exists()
