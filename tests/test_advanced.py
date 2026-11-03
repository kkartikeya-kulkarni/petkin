import json
import numpy as np
import pytest
import SimpleITK as sitk
from petkin import Frames,InputFunction,predict,fit_compartment,logan,ma1,patlak,profile_likelihood,parametric_bootstrap
from petkin.graphical import linear_integral
from petkin.imaging import apply_composed_transforms,motion_displacements
from petkin.bids import frames_from_bids,validate_pet_metadata,input_from_bids


def sample():
    t=np.linspace(0,90,9001);cp=8*t*np.exp(-t/.7)+1.5*np.exp(-t/25)
    edges=np.r_[np.arange(0,5,.25),np.arange(5,20,1),np.arange(20,91,2)]
    return Frames(edges[:-1],edges[1:]),InputFunction(t,cp,.8*cp)


def test_exact_linear_input_integral():
    t=np.array([0.,1.,3.]);v=2*t+1
    q=np.array([-1.,0.,.1,1.5,3.])
    np.testing.assert_allclose(linear_integral(t,v,q),[0,0,.11,3.75,12])
    with pytest.raises(ValueError):linear_integral(t,v,[4.])
    assert linear_integral(t,v,np.array([4.]),tail='hold')[0]==19


def test_vascular_delay_recovery():
    f,inp=sample();p=[.2,.18,.07,.06]
    y=predict(f,inp,p,vb=.08,delay=.18)
    result=fit_compartment(f,inp,y,fit_vb=True,fit_delay=True,input_tail='hold',n_starts=4)
    assert result['success'] and result['usable_for_quantification']
    assert result['vb']==pytest.approx(.08,abs=.003)
    assert result['delay']==pytest.approx(.18,abs=.003)
    np.testing.assert_allclose(list(result['parameters'].values()),p,rtol=.03)


def test_delay_coverage_and_known_noise_contract():
    f,inp=sample()
    with pytest.raises(ValueError,match='cover'):predict(f,inp,[.1,.2],'1t',delay=-.5)
    with pytest.raises(ValueError,match='supplied SDs'):fit_compartment(f,inp,np.ones(len(f.start)),absolute_sigma=True)
    assert np.isfinite(predict(f,inp,[.1,.2],'1t',delay=-.5,input_tail='hold')).all()


def test_bound_estimate_is_not_released_as_valid_VT():
    f,inp=sample();y=predict(f,inp,[.15,.2,.05],model='2t_irreversible')
    result=fit_compartment(f,inp,y,model='2t',n_starts=3)
    assert result['at_bound']
    assert result['VT'] is None and result['raw_VT']>0
    assert not result['usable_for_quantification']
    assert 'parameter_at_bound' in result['qc_reasons']


def test_graphical_methods_known_solutions():
    f,inp=sample();p=[.2,.15]
    y=predict(f,inp,p,'1t',vb=.05)
    assert ma1(f,y,inp,30,vb=.05)['VT']==pytest.approx(.2/.15,rel=.02)
    assert logan(f,y,inp,30,vb=.05)['VT']==pytest.approx(.2/.15,rel=.02)
    irr=predict(f,inp,[.2,.2,.08],'2t_irreversible')
    assert patlak(f,irr,inp,40)['Ki']==pytest.approx(.2*.08/.28,rel=.025)


def test_graphical_zero_weights_and_no_silent_gap():
    f,inp=sample();y=predict(f,inp,[.2,.15],'1t')
    w=np.ones(len(y));w[-1]=0
    result=ma1(f,y,inp,30,weights=w)
    assert len(y)-1 not in result['frame_indices']
    gap=Frames(f.start[1:],f.end[1:])
    with pytest.raises(ValueError,match='contiguous'):ma1(gap,y[1:],inp,30)


def test_profile_closes_on_identifiable_rate():
    f,inp=sample();y=predict(f,inp,[.2,.15],'1t')
    sigma=np.full(len(y),.01);y+=np.random.default_rng(9).normal(0,sigma)
    fit=fit_compartment(f,inp,y,'1t',sigma,absolute_sigma=True,n_starts=2)
    result=profile_likelihood(f,inp,y,fit,'k2',np.linspace(.1,.2,17),sigma=sigma)
    assert result['closed']
    ci=result['interval_components'][0]
    assert ci['lower']<.15<ci['upper']


def test_bootstrap_deterministic_and_retains_all_draws():
    f,inp=sample();sigma=np.full(len(f.start),.01)
    y=predict(f,inp,[.2,.15],'1t')+np.random.default_rng(5).normal(0,sigma)
    fit=fit_compartment(f,inp,y,'1t',sigma,absolute_sigma=True,n_starts=2)
    a=parametric_bootstrap(f,inp,y,fit,sigma=sigma,n=4,seed=8,n_starts=1)
    b=parametric_bootstrap(f,inp,y,fit,sigma=sigma,n=4,seed=8,n_starts=1)
    assert a==b and len(a['draws'])==4


def test_composition_order_and_one_interpolation():
    z,y,x=np.mgrid[:24,:24,:24]
    im=sitk.GetImageFromArray(np.exp(-((x-12)**2+(y-11)**2+(z-13)**2)/20).astype('float32'))
    frame=sitk.Euler3DTransform();frame.SetCenter([12,12,12]);frame.SetRotation(.1,-.08,.12);frame.SetTranslation([1.,2.,0.])
    mri=sitk.TranslationTransform(3,[-2.,1.,.5])
    dynamic=sitk.JoinSeries([im]);dynamic.SetSpacing([1,1,1,30])
    out,txs=apply_composed_transforms(dynamic,im,[frame],mri)
    point=[5.,8.,9.]
    np.testing.assert_allclose(txs[0].TransformPoint(point),frame.TransformPoint(mri.TransformPoint(point)))
    expected=sitk.Resample(im,im,txs[0],sitk.sitkLinear,0.,sitk.sitkFloat32)
    np.testing.assert_array_equal(sitk.GetArrayFromImage(out)[0],sitk.GetArrayFromImage(expected))
    assert out.GetSpacing()[3]==30


def test_displacement_in_physical_units():
    ref=sitk.Image([10,10,10],sitk.sitkFloat32);ref.SetSpacing([2,3,4])
    rows=motion_displacements(ref,[sitk.Transform(3,sitk.sitkIdentity),sitk.TranslationTransform(3,[3.,4.,0.])])
    assert rows[1]['mean_landmark_displacement_mm']==pytest.approx(5)
    assert rows[1]['relative_to_previous_mm']==pytest.approx(5)


def test_bids_units_timing_and_explicit_assumptions(tmp_path):
    meta=dict(FrameTimesStart=[0,30],FrameDuration=[30,60],Units='Bq/mL',ImageDecayCorrected=True,ImageDecayCorrectionTime=0)
    assert validate_pet_metadata(meta)['frames'].end[-1]==1.5
    with pytest.raises(ValueError):validate_pet_metadata(dict(meta,ImageDecayCorrected=False))
    with pytest.raises(ValueError):validate_pet_metadata(dict(meta,ImageDecayCorrectionTime='n/a'))
    path=tmp_path/'blood.tsv';path.write_text('time\tplasma_radioactivity\n10\t2\n60\t1\n')
    bm={'time':{'Units':'s'},'plasma_radioactivity':{'Units':'kBq/mL'}}
    with pytest.raises(ValueError,match='time zero'):input_from_bids(path,bm)
    inp,qc=input_from_bids(path,bm,initial_policy='zero',terminal_time=1.5,terminal_policy='hold')
    assert inp.plasma[1]==2000 and inp.time[-1]==1.5
    assert qc['assumptions']==['initial_zero','terminal_hold']


def test_srtm_independent_ODE_and_reference_logan():
    from scipy.integrate import solve_ivp
    from petkin import ReferenceCurve,predict_srtm,fit_srtm,logan_reference
    t=np.linspace(0,90,9001);cr=5*(np.exp(-t*.025)-np.exp(-t*.3))
    reference=ReferenceCurve(t,cr);edges=np.linspace(0,90,181);f=Frames(edges[:-1],edges[1:])
    r1,k2,bp=1.2,.2,1.5;ka=k2/(1+bp)
    def rhs(time,y):
        c=np.interp(time,t,cr)
        return [c-ka*y[0],r1*c+(k2-r1*ka)*y[0]]
    sol=solve_ivp(rhs,(0,90),[0,0],rtol=1e-10,atol=1e-12,max_step=.05,dense_output=True)
    y=(sol.sol(f.end)[1]-sol.sol(f.start)[1])/f.duration
    np.testing.assert_allclose(predict_srtm(f,reference,[r1,k2,bp]),y,rtol=.002,atol=2e-4)
    result=fit_srtm(f,reference,y,n_starts=3)
    assert result['BPnd']==pytest.approx(bp,rel=.01)
    assert logan_reference(f,reference,y,30,k2/r1)['DVR']==pytest.approx(1+bp,rel=.01)


def test_reference_models_reject_invalid_region_data():
    from petkin import ReferenceCurve,predict_srtm,logan_reference
    with pytest.raises(ValueError):ReferenceCurve([0,1],[-1,2])
    ref=ReferenceCurve([0,1],[0,1]);f=Frames([0],[1])
    with pytest.raises(ValueError):predict_srtm(f,ref,[1,.1,-1])
    with pytest.raises(ValueError):logan_reference(f,ref,[1],0,-1)


def test_rank_deficient_fit_does_not_claim_kinetic_estimate():
    f=Frames(np.arange(6.),np.arange(1.,7.));inp=InputFunction([0,6],[0,0])
    fit=fit_compartment(f,inp,np.ones(6),'1t',n_starts=2)
    assert fit['success'] and fit['VT'] is None
    assert 'rank_deficient_jacobian' in fit['qc_reasons']
    json.dumps(fit,allow_nan=False)


def test_custom_rate_bounds_above_one_are_supported():
    t=np.linspace(0,3,301);inp=InputFunction(t,np.ones_like(t));edges=np.linspace(0,3,31);f=Frames(edges[:-1],edges[1:])
    y=predict(f,inp,[2.7,2.5],'1t')
    fit=fit_compartment(f,inp,y,'1t',rate_bounds=(2.,4.),n_starts=4)
    np.testing.assert_allclose(list(fit['parameters'].values()),[2.7,2.5],rtol=.01)


def test_html_withholds_unreliable_value_and_escapes_title(tmp_path):
    from petkin.reporting import write_qc_html,file_manifest
    path=tmp_path/'report.html'
    write_qc_html(dict(model='2t',VT=None,raw_VT=10000,usable_for_quantification=False,qc_reasons=['parameter_at_bound']),path,'<script>alert(1)</script>')
    text=path.read_text()
    assert '<script>' not in text and 'Withheld' in text and 'Review required' in text
    f=tmp_path/'input.txt';f.write_text('abc')
    manifest=file_manifest([f])
    assert manifest['inputs'][0]['sha256']=='ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'


def test_landmark_QC_requires_physical_correspondence():
    from petkin.imaging import landmark_registration_error
    tx=sitk.TranslationTransform(3,[3,4,0]);fixed=np.array([[0,0,0],[20,0,0],[0,20,0]],float);moving=fixed+[3,4,0]
    result=landmark_registration_error(fixed,moving,tx)
    assert result['rms_error_mm']==0 and result['unregistered_rms_error_mm']==5
    with pytest.raises(ValueError):landmark_registration_error(fixed[:2],moving[:2],tx)
