import copy,json
import numpy as np
import pytest
from petkin import Frames,InputFunction,fit_compartment,profile_macroparameter,analyze_research_tac,ResearchPolicy,qualification_status
from petkin.validation import REQUIRED_ASSUMPTIONS,assess_research_fit


def example(model='1t'):
    from petkin import predict
    time=np.linspace(0,90,9001);inp=InputFunction(time,8*time*np.exp(-time/.7)+1.5*np.exp(-time/25))
    edges=np.r_[np.arange(0,5,.5),np.arange(5,91,1.)];frames=Frames(edges[:-1],edges[1:])
    p={'1t':[.2,.15],'2t':[.2,.2,.06,.05],'2t_irreversible':[.2,.2,.06]}[model]
    y=predict(frames,inp,p,model);sigma=.02*np.maximum(y,.05*y.max());y+=np.random.default_rng(11).normal(0,sigma)
    return frames,inp,y,sigma


@pytest.mark.parametrize('model',['1t','2t','2t_irreversible'])
def test_macro_profile_preserves_exact_constraint_and_closes(model):
    f,inp,y,sigma=example(model);fit=fit_compartment(f,inp,y,model,sigma,n_starts=3,absolute_sigma=True)
    key='Ki' if model=='2t_irreversible' else 'VT';value=fit['raw_'+key]
    profile=profile_macroparameter(f,inp,y,fit,np.geomspace(value/5,value*5,21),sigma=sigma)
    assert profile['parameter']==key and profile['closed']
    for row in profile['grid']:
        p=row['parameters']
        macro=p['K1']*p['k3']/(p['k2']+p['k3']) if model=='2t_irreversible' else p['K1']/p['k2']*(1+p['k3']/p['k4'] if model=='2t' else 1)
        assert macro==pytest.approx(row['value'],rel=1e-10)
    interval=profile['interval_components'][0]
    assert interval['lower']<value<interval['upper']


def test_missing_assumptions_withhold_even_precise_fit():
    f,inp,y,sigma=example();fit=fit_compartment(f,inp,y,'1t',sigma,n_starts=2,absolute_sigma=True)
    point=fit['raw_VT'];profile=profile_macroparameter(f,inp,y,fit,np.geomspace(point/5,point*5,21),sigma=sigma)
    assessed=assess_research_fit(fit,profile,{})
    assert assessed['VT'] is None and not assessed['research_screen']['passed']
    assert assessed['research_screen']['unconfirmed_assumptions']==list(REQUIRED_ASSUMPTIONS)


def test_wrong_target_and_open_profiles_cannot_release_estimate():
    f,inp,y,sigma=example();fit=fit_compartment(f,inp,y,'1t',sigma,n_starts=2,absolute_sigma=True)
    point=fit['raw_VT'];profile=profile_macroparameter(f,inp,y,fit,np.geomspace(point/5,point*5,21),sigma=sigma)
    assumptions={k:True for k in REQUIRED_ASSUMPTIONS}
    wrong=copy.deepcopy(profile);wrong['parameter']='k2'
    assert assess_research_fit(fit,wrong,assumptions)['VT'] is None
    wrong=copy.deepcopy(profile);wrong['closed']=False
    assert assess_research_fit(fit,wrong,assumptions)['VT'] is None


def test_gate_report_fails_closed_for_missing_or_blocked_evidence():
    protocol={'protocol_id':'test','required_gates':[{'id':'software'},{'id':'human_truth'}]}
    report=qualification_status(protocol,[{'id':'software','status':'passed'}])
    assert not report['qualified_for_declared_research_scope'] and report['missing_gates']==['human_truth']
    report=qualification_status(protocol,[{'id':'software','status':'passed'},{'id':'human_truth','status':'blocked'}])
    assert report['overall_status']=='not_qualified'
    with pytest.raises(ValueError):qualification_status(protocol,[{'id':'x','status':'assumed_passed'}])
    with pytest.raises(ValueError):qualification_status(protocol,[{'id':'x','status':'passed'},{'id':'x','status':'passed'}])


def test_wide_interval_and_correlated_residuals_require_review():
    f,inp,y,sigma=example();fit=fit_compartment(f,inp,y,'1t',sigma,n_starts=2,absolute_sigma=True)
    point=fit['raw_VT'];profile=profile_macroparameter(f,inp,y,fit,np.geomspace(point/5,point*5,21),sigma=sigma)
    assumptions={k:True for k in REQUIRED_ASSUMPTIONS}
    result=assess_research_fit(fit,profile,assumptions)
    assert result['research_screen']['passed']
    wide=copy.deepcopy(profile);wide['interval_components'][0].update(lower=point/10,upper=point*10)
    assert not assess_research_fit(fit,wide,assumptions)['research_screen']['passed']
    serial=copy.deepcopy(fit);serial['residual_diagnostics']['lag1_correlation']=.9
    result=assess_research_fit(serial,profile,assumptions)
    assert result['VT'] is None and 'structured_residuals_require_review' in result['research_screen']['reasons']


def test_macro_interval_crossings_are_refitted_not_linearly_interpolated():
    f,inp,y,sigma=example();fit=fit_compartment(f,inp,y,'1t',sigma,n_starts=2,absolute_sigma=True)
    point=fit['raw_VT'];profile=profile_macroparameter(f,inp,y,fit,np.geomspace(point/20,point*20,11),sigma=sigma)
    assert profile['closed']
    component=profile['interval_components'][0]
    for endpoint in [component['lower'],component['upper']]:
        closest=min(profile['grid'],key=lambda r:abs(r['value']-endpoint))
        assert closest['value']==pytest.approx(endpoint,rel=1e-8)
        assert closest['likelihood_ratio']==pytest.approx(profile['threshold'],abs=.002)
