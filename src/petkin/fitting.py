"""Bounded multistart fitting with separate optimizer and quantification status."""
import numpy as np
from scipy.optimize import least_squares
from .kinetics import vector, PARAMS, predict


def _macro(p, model):
    if model == '1t': return float(p[0]/p[1])
    if model == '2t': return float(p[0]/p[1]*(1+p[2]/p[3]))
    return float(p[0]*p[2]/(p[1]+p[2]))


def _macro_gradient(p,model):
    g=np.zeros(len(p));k1,k2=p[:2]
    if model=='1t':g[:2]=[1/k2,-k1/k2**2]
    elif model=='2t':
        k3,k4=p[2:4];vt=_macro(p,model)
        g[:4]=[vt/k1,-vt/k2,k1/k2/k4,-k1*k3/k2/k4**2]
    else:
        k3=p[2];g[:3]=[k3/(k2+k3),-k1*k3/(k2+k3)**2,k1*k2/(k2+k3)**2]
    return g


def fit_compartment(frames, input_function, tissue, model='2t', sigma=None, vb=0., dt=.025,
                    *, fit_vb=False, fit_delay=False, delay=0., delay_bounds=(-.5,.5),
                    input_tail='error', n_starts=8, seed=0, rate_bounds=(1e-6,3.),
                    vb_bounds=(0.,.3), absolute_sigma=False, max_nfev=800):
    """Gaussian WLS using log-rate coordinates and deterministic multistarts.

    Positive delay means input arrives later in tissue: Cp(t-delay). Plasma and
    whole blood share that delay. Negative delay requires input coverage beyond
    scan end, or an explicit terminal-hold assumption. sigma defaults to relative
    weights; absolute_sigma=True treats supplied SDs as known calibrated SDs.

    VT/Ki is null when automatic checks fail; raw_VT/raw_Ki retains the numerical
    estimate for audit. QC thresholds are engineering flags, not clinical limits.
    """
    if model not in PARAMS:raise ValueError('Unknown model')
    y=vector(tissue,'tissue')
    if y.shape!=frames.start.shape:raise ValueError('Tissue must match frame count')
    if np.all(np.abs(y)<1e-15):raise ValueError('No tissue signal to estimate kinetics')
    s=np.ones_like(y) if sigma is None else vector(sigma,'sigma')
    if s.shape!=y.shape or np.any(s<=0):raise ValueError('sigma must be positive and match frames')
    if absolute_sigma and sigma is None:raise ValueError('absolute_sigma requires supplied SDs')
    if not isinstance(n_starts,(int,np.integer)) or n_starts<1:raise ValueError('n_starts must be a positive integer')
    if not isinstance(max_nfev,(int,np.integer)) or max_nfev<1:raise ValueError('max_nfev must be positive')
    nr=len(PARAMS[model]);names=list(PARAMS[model])
    lo,hi=map(float,rate_bounds)
    if not np.isfinite([lo,hi]).all() or not 0<lo<hi:raise ValueError('Positive ordered rate bounds required')
    lower=[lo]*nr;upper=[hi]*nr
    if not 0<=vb<1 or not np.isfinite(delay):raise ValueError('Invalid fixed vb/delay')
    if fit_vb:
        if not np.isfinite(vb_bounds).all() or not 0<=vb_bounds[0]<vb_bounds[1]<1:raise ValueError('Invalid vb bounds')
        names.append('vb');lower.append(vb_bounds[0]);upper.append(vb_bounds[1])
    if (vb>0 or fit_vb) and input_function.blood is None:raise ValueError('Whole blood is required for vascular modeling')
    if fit_delay:
        if not np.isfinite(delay_bounds).all() or delay_bounds[0]>=delay_bounds[1]:raise ValueError('Invalid delay bounds')
        names.append('delay');lower.append(delay_bounds[0]);upper.append(delay_bounds[1])
    lower,upper=np.array(lower,float),np.array(upper,float)
    if len(y)<=len(names):raise ValueError('Insufficient observations for free parameters')
    # Validate input coverage at the entire requested delay range up front.
    check_delay=delay_bounds[0] if fit_delay else delay
    predict(frames,input_function,[.1]*nr,model,vb,dt,delay=check_delay,input_tail=input_tail)
    def encode(p):
        z=np.array(p,float).copy();z[:nr]=np.log(z[:nr]);return z
    def decode(z):
        p=np.array(z,float).copy();p[:nr]=np.exp(p[:nr]);return p
    def unpack(p):
        v=float(p[names.index('vb')]) if fit_vb else vb
        d=float(p[names.index('delay')]) if fit_delay else delay
        return v,d
    def residual(z):
        p=decode(z);v,d=unpack(p)
        return (predict(frames,input_function,p[:nr],model,v,dt,delay=d,input_tail=input_tail)-y)/s
    base={'1t':[[.1,.1],[.4,.3]],'2t':[[.1,.15,.05,.03],[.3,.3,.15,.1],[.2,.1,.01,.01]],
          '2t_irreversible':[[.1,.1,.05],[.3,.3,.15]]}[model]
    rng=np.random.default_rng(seed);starts=[]
    for i in range(n_starts):
        sample_lo,sample_hi=max(lo,.003),min(hi,1.)
        if sample_lo>=sample_hi:sample_lo,sample_hi=lo,hi
        rates=np.array(base[i]) if i<len(base) else np.exp(rng.uniform(np.log(sample_lo),np.log(sample_hi),nr))
        p=list(rates)
        if fit_vb:p.append(.05 if i==0 else rng.uniform(*vb_bounds))
        if fit_delay:p.append(np.clip(delay,*delay_bounds) if i==0 else rng.uniform(*delay_bounds))
        margin=(upper-lower)*1e-8
        starts.append(np.clip(p,lower+margin,upper-margin))
    candidates=[least_squares(residual,encode(p),bounds=(encode(lower),encode(upper)),
                              x_scale='jac',max_nfev=max_nfev,ftol=1e-9,xtol=1e-9,gtol=1e-9) for p in starts]
    best=min(candidates,key=lambda f:f.cost);p=decode(best.x);v,d=unpack(p)
    prediction=predict(frames,input_function,p[:nr],model,v,dt,delay=d,input_tail=input_tail)
    sse=float(2*best.cost);n=len(y);npfree=len(p)
    dzdp=np.ones(npfree);dzdp[:nr]=1/p[:nr]
    j=best.jac*dzdp
    rank=int(np.linalg.matrix_rank(j))
    covariance_z=np.linalg.pinv(best.jac.T@best.jac)
    dpdz=1/dzdp;covariance=covariance_z*np.outer(dpdz,dpdz)
    if not absolute_sigma:covariance*=sse/(n-npfree)
    scales=np.r_[np.ones(nr),([.05] if fit_vb else []),([.1] if fit_delay else [])]
    condition=float(np.linalg.cond(best.jac*scales))
    boundary=(best.active_mask!=0)|((p-lower)<np.maximum(1e-10,(upper-lower)*1e-7))|((upper-p)<np.maximum(1e-10,(upper-lower)*1e-7))
    value=_macro(p,model);gradient=_macro_gradient(p,model)
    macro_se=float(np.sqrt(max(0.,gradient@covariance@gradient)))
    reasons=[]
    if not best.success:reasons.append('optimizer_not_converged')
    if np.any(boundary):reasons.append('parameter_at_bound')
    if rank<npfree:reasons.append('rank_deficient_jacobian')
    if not np.isfinite(condition) or condition>1e8:reasons.append('ill_conditioned_scaled_jacobian')
    if macro_se/max(abs(value),1e-15)>.5:reasons.append('large_local_macro_uncertainty')
    if not np.isfinite(value) or value<=0:reasons.append('nonpositive_or_nonfinite_macroparameter')
    if sse<=np.finfo(float).tiny:reasons.append('residual_variance_not_estimable' if not absolute_sigma else 'near_exact_fit')
    # A noiseless fit is numerically valid; zero estimated variance is not evidence
    # of accurate real-data uncertainty and is documented without suppressing it.
    blocking=[r for r in reasons if r not in ('near_exact_fit','residual_variance_not_estimable')]
    macro='Ki' if model=='2t_irreversible' else 'VT'
    variance_parameters=0 if absolute_sigma else 1;k=npfree+variance_parameters
    deviance=sse+n*np.log(2*np.pi)+2*np.log(s).sum() if absolute_sigma else n*(np.log(max(sse/n,1e-300))+1+np.log(2*np.pi))+2*np.log(s).sum()
    result=dict(model=model,parameters=dict(zip(PARAMS[model],map(float,p[:nr]))),
                success=bool(best.success),message=best.message,status='ok' if not blocking else 'review_required',
                usable_for_quantification=not blocking,qc_reasons=reasons,
                qc_scope='Local engineering checks only; profile, bootstrap, input and biological validation remain necessary',
                prediction=prediction.tolist(),residuals=(prediction-y).tolist(),weighted_sse=sse,
                rmse=float(np.sqrt(np.mean((prediction-y)**2))),
                standard_errors=dict(zip(names,np.sqrt(np.maximum(0,np.diag(covariance))).tolist())),
                covariance=covariance.tolist(),free_parameters=names,jacobian_rank=rank,
                condition_number=float(np.linalg.cond(j)) if np.isfinite(np.linalg.cond(j)) else None,
                scaled_condition_number=condition if np.isfinite(condition) else None,
                at_bound=bool(np.any(boundary)),bound_parameters=[name for name,b in zip(names,boundary) if b],
                vb=float(v),delay=float(d),macro_standard_error=macro_se,
                macro_local_95_interval=[float(value-1.96*macro_se),float(value+1.96*macro_se)],
                AIC=float(deviance+2*k),AICc=float(deviance+2*k+2*k*(k+1)/(n-k-1)) if n>k+1 else None,
                BIC=float(deviance+k*np.log(n)),
                settings=dict(dt=dt,fit_vb=fit_vb,fit_delay=fit_delay,delay_bounds=list(delay_bounds),
                              input_tail=input_tail,n_starts=n_starts,seed=seed,rate_bounds=list(rate_bounds),
                              vb_bounds=list(vb_bounds),absolute_sigma=absolute_sigma,max_nfev=max_nfev),
                multistart=dict(attempted=n_starts,converged=sum(bool(f.success) for f in candidates),
                                weighted_sse=[float(2*f.cost) for f in candidates]))
    weighted_residual=(prediction-y)/s
    residual_power=float(weighted_residual@weighted_residual)
    lag=None
    if len(y)>2 and np.std(weighted_residual[:-1])>1e-14 and np.std(weighted_residual[1:])>1e-14:
        lag=float(np.corrcoef(weighted_residual[:-1],weighted_residual[1:])[0,1])
    result['residual_diagnostics']=dict(lag1_correlation=lag,
        durbin_watson=float(np.sum(np.diff(weighted_residual)**2)/residual_power) if residual_power>1e-30 else None,
        weighted_sse_per_residual_dof=sse/(n-npfree),
        interpretation='Serial structure can reflect model/input mismatch or correlated noise; no calibrated pass/fail threshold is asserted')
    result['raw_'+macro]=value;result[macro]=value if not blocking else None
    return result
