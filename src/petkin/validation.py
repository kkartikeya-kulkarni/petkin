"""Fail-closed research release screening and explicit qualification gates."""
from dataclasses import asdict,dataclass
import copy,hashlib,json
import numpy as np
from .fitting import fit_compartment
from .uncertainty import profile_macroparameter


@dataclass(frozen=True)
class ResearchPolicy:
    confidence: float=.95
    profile_points: int=41
    profile_span: float=20.
    max_relative_profile_half_width: float=.5
    max_absolute_residual_lag1: float=.6
    minimum_residual_frames: int=12
    def __post_init__(self):
        if not 0<self.confidence<1 or self.profile_points<5 or self.profile_span<=1 or self.max_relative_profile_half_width<=0 or not 0<self.max_absolute_residual_lag1<1:
            raise ValueError('Invalid research precision policy')
    @property
    def sha256(self):return hashlib.sha256(json.dumps(asdict(self),sort_keys=True).encode()).hexdigest()


REQUIRED_ASSUMPTIONS=('calibration_matched','decay_reference_matched','parent_correction_appropriate',
                      'frame_timing_verified','input_tissue_timing_verified','vascular_model_justified','ROI_geometry_verified')


def assess_research_fit(fit,profile,assumptions,policy=None):
    """Release an estimate only with local QC, a precise macro profile and inputs.

    This is a per-result screen, not toolbox-wide external validation. Caller
    assertions about inputs must be backed by study records, not inferred here.
    """
    policy=policy or ResearchPolicy();result=copy.deepcopy(fit);macro='Ki' if fit['model']=='2t_irreversible' else 'VT'
    reasons=[]
    if not fit.get('usable_for_quantification',False):reasons.extend('local:'+r for r in fit.get('qc_reasons',['local_QC_not_passed']))
    unconfirmed=[name for name in REQUIRED_ASSUMPTIONS if assumptions.get(name) is not True]
    if unconfirmed:reasons.append('unconfirmed_input_assumptions')
    if profile is None:reasons.append('macro_profile_not_available')
    else:
        if profile.get('parameter')!=macro:reasons.append('wrong_profile_target')
        if not np.isclose(profile.get('confidence',0),policy.confidence):reasons.append('profile_confidence_mismatch')
        if not profile.get('closed',False):reasons.append('macro_profile_open_or_incomplete')
        if profile.get('profile_found_lower_objective',False):reasons.append('profile_improved_on_full_fit')
        components=profile.get('interval_components',[])
        if len(components)!=1:reasons.append('macro_profile_not_single_component')
        elif components[0].get('lower') is not None and components[0].get('upper') is not None:
            lo,hi=components[0]['lower'],components[0]['upper'];point=fit['raw_'+macro]
            if not lo<=point<=hi:reasons.append('point_not_in_profile_interval')
            if (hi-lo)/(2*point)>policy.max_relative_profile_half_width:reasons.append('macro_profile_too_wide')
    lag=fit.get('residual_diagnostics',{}).get('lag1_correlation')
    if len(fit['prediction'])>=policy.minimum_residual_frames and lag is not None and abs(lag)>policy.max_absolute_residual_lag1:
        reasons.append('structured_residuals_require_review')
    result['local_qc_passed']=bool(fit.get('usable_for_quantification',False))
    result['research_screen']=dict(passed=not reasons,reasons=reasons,unconfirmed_assumptions=unconfirmed,
                                   policy=asdict(policy),policy_sha256=policy.sha256,
                                   scope='Per-result conditional screening; not toolbox validation or biological accuracy certification')
    result['usable_for_quantification']=not reasons;result['status']='research_screen_passed' if not reasons else 'review_required'
    result[macro]=fit['raw_'+macro] if not reasons else None
    return result


def analyze_research_tac(frames,input_function,tissue,model='2t',sigma=None,*,assumptions=None,policy=None,**fit_options):
    policy=policy or ResearchPolicy();fit=fit_compartment(frames,input_function,tissue,model,sigma,**fit_options)
    macro='Ki' if model=='2t_irreversible' else 'VT';profile=None
    if fit['usable_for_quantification']:
        point=fit['raw_'+macro]
        values=np.geomspace(max(1e-12,point/policy.profile_span),point*policy.profile_span,policy.profile_points)
        profile=profile_macroparameter(frames,input_function,tissue,fit,values,sigma=sigma,confidence=policy.confidence)
    result=assess_research_fit(fit,profile,assumptions or {},policy)
    result['macro_profile']=profile
    return result


def qualification_status(protocol,gates):
    """Missing, failed or unavailable mandatory evidence can never yield passed."""
    valid={'passed','failed','blocked','not_run'}
    if any(r.get('status') not in valid for r in gates):raise ValueError('Unknown qualification gate status')
    identifiers=[r['id'] for r in gates]
    if len(set(identifiers))!=len(identifiers):raise ValueError('Duplicate qualification gate')
    observed={r['id']:r for r in gates};required=[r['id'] for r in protocol['required_gates']]
    missing=[key for key in required if key not in observed]
    incomplete=[key for key in required if observed.get(key,{}).get('status')!='passed']
    return dict(protocol_id=protocol['protocol_id'],qualified_for_declared_research_scope=not incomplete,
                overall_status='passed' if not incomplete else 'not_qualified',missing_gates=missing,
                incomplete_gates=incomplete,gates=gates,
                claim='All mandatory gates passed' if not incomplete else 'Full research validation is not established')
