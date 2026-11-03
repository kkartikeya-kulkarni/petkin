"""Conditional profile likelihood and parametric-bootstrap diagnostics."""
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
from .kinetics import vector,PARAMS,predict
from .fitting import fit_compartment,_macro


def profile_likelihood(frames,input_function,tissue,fit,parameter,values,*,sigma=None,confidence=.95):
    """Grid profile with nuisance refitting; missing crossings are reported open.

    Known Gaussian SDs: delta SSE. Unknown residual scale/relative weights:
    N*log(SSE/SSEmin). Chi-square thresholds are asymptotic and can be unreliable
    at boundaries or in small samples. Grid endpoints are never called closed CIs.
    """
    if not 0<confidence<1:raise ValueError('confidence must be between zero and one')
    names=fit['free_parameters'];model=fit['model'];nr=len(PARAMS[model])
    if parameter not in names:raise ValueError('Profile parameter must have been fitted')
    settings=fit['settings'];y=vector(tissue,'tissue');s=np.ones_like(y) if sigma is None else vector(sigma,'sigma')
    if y.shape!=frames.start.shape or s.shape!=y.shape or np.any(s<=0):raise ValueError('Invalid profile observations or sigma')
    if settings['absolute_sigma'] and sigma is None:raise ValueError('Known-SD profile needs the fitted SD vector')
    p=np.array([fit['parameters'][name] if name in fit['parameters'] else fit[name] for name in names])
    lo=[settings['rate_bounds'][0]]*nr;hi=[settings['rate_bounds'][1]]*nr
    if settings['fit_vb']:lo.append(settings['vb_bounds'][0]);hi.append(settings['vb_bounds'][1])
    if settings['fit_delay']:lo.append(settings['delay_bounds'][0]);hi.append(settings['delay_bounds'][1])
    lo,hi=np.array(lo),np.array(hi);index=names.index(parameter)
    grid=np.unique(np.r_[vector(values,'values'),p[index]])
    if np.any(grid<lo[index]) or np.any(grid>hi[index]):raise ValueError('Profile grid must be within fitted bounds')
    free=np.arange(len(names))!=index
    def encode(q):
        z=q.copy();z[:nr]=np.log(z[:nr]);return z
    def decode(z):
        q=z.copy();q[:nr]=np.exp(q[:nr]);return q
    initial=encode(p);zl,zh=encode(lo),encode(hi)
    rows=[]
    for value in grid:
        zfixed=float(np.log(value)) if index<nr else float(value)
        def residual(zfree):
            z=initial.copy();z[index]=zfixed;z[free]=zfree;q=decode(z)
            vb=float(q[names.index('vb')]) if 'vb' in names else fit['vb']
            delay=float(q[names.index('delay')]) if 'delay' in names else fit['delay']
            return (predict(frames,input_function,q[:nr],model,vb,settings['dt'],delay=delay,input_tail=settings['input_tail'])-y)/s
        trials=[]
        for start in [initial[free],(initial[free]+(zl[free]+zh[free])/2)/2]:
            trials.append(least_squares(residual,start,bounds=(zl[free],zh[free]),x_scale='jac',max_nfev=settings['max_nfev']))
        best=min(trials,key=lambda f:f.cost)
        rows.append(dict(value=float(value),weighted_sse=float(2*best.cost),converged=bool(best.success)))
    minimum=min(fit['weighted_sse'],min(r['weighted_sse'] for r in rows))
    for r in rows:
        r['likelihood_ratio']=max(0.,r['weighted_sse']-minimum) if settings['absolute_sigma'] else float(len(y)*np.log(max(r['weighted_sse'],1e-300)/max(minimum,1e-300)))
    cutoff=float(chi2.ppf(confidence,1));accepted=np.array([r['likelihood_ratio']<=cutoff and r['converged'] for r in rows])
    # Report all accepted grid components. Crossings are linearly interpolated,
    # never extrapolated outside the evaluated grid or optimizer parameter box.
    components=[]
    starts=np.flatnonzero(accepted & ~np.r_[False,accepted[:-1]])
    ends=np.flatnonzero(accepted & ~np.r_[accepted[1:],False])
    for a,b in zip(starts,ends):
        components.append(dict(lower=None if a==0 or not rows[a-1]['converged'] else rows[a-1]['value'],upper=None if b==len(rows)-1 or not rows[b+1]['converged'] else rows[b+1]['value'],
                               open_lower=bool(a==0 or not rows[a-1]['converged']),open_upper=bool(b==len(rows)-1 or not rows[b+1]['converged']),accepted_grid_range=[rows[a]['value'],rows[b]['value']]))
    return dict(parameter=parameter,confidence=confidence,threshold=cutoff,grid=rows,interval_components=components,
                closed=bool(components and all(r['converged'] for r in rows) and all(not r['open_lower'] and not r['open_upper'] for r in components)),
                all_nuisance_fits_converged=all(r['converged'] for r in rows),
                profile_found_lower_objective=bool(minimum<fit['weighted_sse']-max(1e-8,fit['weighted_sse']*1e-6)),
                relative_weight_scale_profiled=not settings['absolute_sigma'],
                interval_method='Conservative outer crossing brackets on sampled grid; no linear interpolation of likelihood ratio',
                warning='Conditional on input, model, vascular fraction and preprocessing; asymptotic chi-square approximation')


def parametric_bootstrap(frames,input_function,tissue,fit,*,sigma=None,n=100,seed=0,n_starts=3):
    """Gaussian conditional bootstrap; retains invalid draws in raw intervals."""
    if not isinstance(n,int) or n<2:raise ValueError('At least two bootstrap samples required')
    y=vector(tissue,'tissue');s=np.ones_like(y) if sigma is None else vector(sigma,'sigma')
    if y.shape!=frames.start.shape or s.shape!=y.shape or np.any(s<=0):raise ValueError('Invalid bootstrap observations or SDs')
    settings=fit['settings'].copy()
    if settings['absolute_sigma'] and sigma is None:raise ValueError('Known-SD bootstrap needs supplied SDs')
    scale=1. if settings['absolute_sigma'] else np.sqrt(fit['weighted_sse']/(len(y)-len(fit['free_parameters'])))
    rng=np.random.default_rng(seed);draws=[];key='Ki' if fit['model']=='2t_irreversible' else 'VT'
    settings['n_starts']=n_starts
    for i in range(n):
        observed=np.array(fit['prediction'])+rng.normal(0,s*scale)
        settings['seed']=seed+i+1
        result=fit_compartment(frames,input_function,observed,fit['model'],s,fit['vb'],
                               delay=fit['delay'],**settings)
        draws.append(dict(raw_estimate=result['raw_'+key],usable=result['usable_for_quantification'],
                          success=result['success'],at_bound=result['at_bound'],qc_reasons=result['qc_reasons'],parameters=result['parameters']))
    values=[r['raw_estimate'] for r in draws]
    return dict(method='conditional_parametric_Gaussian_bootstrap',n=n,seed=seed,draws=draws,
                raw_percentile_95_interval=np.quantile(values,[.025,.975]).tolist(),
                unusable_fraction=float(np.mean([not r['usable'] for r in draws])),
                warning='Percentiles retain boundary/invalid estimates; conditional noise model excludes input and registration uncertainty')


def profile_macroparameter(frames,input_function,tissue,fit,values,*,sigma=None,confidence=.95):
    """Profile VT/Ki directly, refitting rates under the exact macro constraint.

    K1 is eliminated algebraically. SLSQP enforces its original rate bounds;
    all other original fitted bounds remain active. Open, multimodal, failed or
    box-limited intervals are reported explicitly. Conditional/asymptotic only.
    """
    from scipy.optimize import minimize
    if not 0<confidence<1:raise ValueError('confidence must lie between zero and one')
    model=fit['model'];nr=len(PARAMS[model]);names=fit['free_parameters'];settings=fit['settings']
    y=vector(tissue,'tissue');s=np.ones_like(y) if sigma is None else vector(sigma,'sigma')
    if y.shape!=frames.start.shape or s.shape!=y.shape or np.any(s<=0):raise ValueError('Invalid observations or SDs')
    if settings['absolute_sigma'] and sigma is None:raise ValueError('Known-SD macro profile requires supplied SDs')
    base=np.array([fit['parameters'][name] if name in fit['parameters'] else fit[name] for name in names])
    lo=[settings['rate_bounds'][0]]*nr;hi=[settings['rate_bounds'][1]]*nr
    if settings['fit_vb']:lo.append(settings['vb_bounds'][0]);hi.append(settings['vb_bounds'][1])
    if settings['fit_delay']:lo.append(settings['delay_bounds'][0]);hi.append(settings['delay_bounds'][1])
    lo,hi=np.array(lo,float),np.array(hi,float)
    macro='Ki' if model=='2t_irreversible' else 'VT';point=fit['raw_'+macro]
    grid=np.unique(np.r_[vector(values,'macro values'),point])
    if np.any(grid<=0):raise ValueError('Positive macroparameter grid required')
    def encode(q):
        z=q.copy();z[:nr]=np.log(z[:nr]);return z[1:]
    def decode(z,value):
        q=base.copy();q[1:]=z;q[1:nr]=np.exp(q[1:nr])
        if model=='1t':q[0]=value*q[1]
        elif model=='2t':q[0]=value*q[1]/(1+q[2]/q[3])
        else:q[0]=value*(q[1]+q[2])/q[2]
        return q
    bounds=list(zip(encode(lo),encode(hi)));rows=[];normalizer=max(fit['weighted_sse'],1e-6)
    cache={}
    def evaluate(value):
        value=float(value)
        if value in cache:return cache[value]
        def residual(z):
            q=decode(z,value)
            v=float(q[names.index('vb')]) if settings['fit_vb'] else fit['vb']
            d=float(q[names.index('delay')]) if settings['fit_delay'] else fit['delay']
            return (predict(frames,input_function,q[:nr],model,v,settings['dt'],delay=d,input_tail=settings['input_tail'])-y)/s
        def objective(z):
            r=residual(z);return float(r@r/normalizer)
        def constraint(z):
            q=decode(z,value);return np.array([q[0]-lo[0],hi[0]-q[0]])
        adjusted=base.copy()
        if model in ('1t','2t'):adjusted[1]=np.clip(base[1]*point/value,lo[1]*(1+1e-8),hi[1]*(1-1e-8))
        else:adjusted[1]=np.clip(base[2]*(base[0]/value-1),lo[1]*(1+1e-8),hi[1]*(1-1e-8))
        starts=[encode(adjusted),encode(base)];trials=[]
        for start in starts:
            result=minimize(objective,start,method='SLSQP',bounds=bounds,
                            constraints=[dict(type='ineq',fun=constraint)],
                            options=dict(maxiter=settings['max_nfev'],ftol=1e-9))
            q=decode(result.x,value);feasible=bool(lo[0]-max(1e-12,abs(lo[0])*1e-6)<=q[0]<=hi[0]+max(1e-12,abs(hi[0])*1e-6))
            trials.append((result,feasible))
        valid=[t for t in trials if t[1] and np.isfinite(t[0].fun)]
        best,feasible=min(valid or trials,key=lambda t:t[0].fun)
        q=decode(best.x,value)
        row=dict(value=float(value),weighted_sse=float(best.fun*normalizer),
                 converged=bool(best.success and feasible),message=str(best.message),
                 parameters=dict(zip(names,map(float,q))))
        cache[value]=row
        return row
    rows=[evaluate(float(value)) for value in grid]
    minimum=min(fit['weighted_sse'],min(r['weighted_sse'] for r in rows if r['converged'])) if any(r['converged'] for r in rows) else fit['weighted_sse']
    for r in rows:
        r['likelihood_ratio']=float(max(0.,r['weighted_sse']-minimum)) if settings['absolute_sigma'] else float(max(0.,len(y)*np.log(max(r['weighted_sse'],1e-300)/max(minimum,1e-300))))
    cutoff=float(chi2.ppf(confidence,1));accepted=np.array([r['converged'] and r['likelihood_ratio']<=cutoff for r in rows])
    first=np.flatnonzero(accepted&~np.r_[False,accepted[:-1]]);last=np.flatnonzero(accepted&~np.r_[accepted[1:],False]);components=[]
    def cross(i,j):
        from scipy.optimize import brentq
        def objective(log_value):
            row=evaluate(float(np.exp(log_value)))
            if not row['converged']:raise ValueError('Nuisance optimization failed during crossing refinement')
            score=max(0.,row['weighted_sse']-minimum) if settings['absolute_sigma'] else max(0.,len(y)*np.log(max(row['weighted_sse'],1e-300)/max(minimum,1e-300)))
            row['likelihood_ratio']=float(score)
            return score-cutoff
        try:
            return float(np.exp(brentq(objective,np.log(rows[i]['value']),np.log(rows[j]['value']),xtol=1e-6)))
        except (ValueError,RuntimeError):return None
    for a,b in zip(first,last):
        open_lo=a==0 or not rows[a-1]['converged'];open_hi=b==len(rows)-1 or not rows[b+1]['converged']
        lower=None if open_lo else cross(a-1,a);upper=None if open_hi else cross(b,b+1)
        components.append(dict(lower=lower,upper=upper,open_lower=bool(open_lo or lower is None),open_upper=bool(open_hi or upper is None),
                               accepted_grid_range=[rows[a]['value'],rows[b]['value']]))
    rows=sorted(cache.values(),key=lambda r:r['value'])
    for row in rows:
        if 'likelihood_ratio' not in row:
            row['likelihood_ratio']=float(max(0.,row['weighted_sse']-minimum)) if settings['absolute_sigma'] else float(max(0.,len(y)*np.log(max(row['weighted_sse'],1e-300)/max(minimum,1e-300))))
    all_converged=all(r['converged'] for r in rows)
    best_discovered=min([minimum]+[r['weighted_sse'] for r in rows if r['converged']])
    return dict(parameter=macro,confidence=confidence,threshold=cutoff,grid=rows,interval_components=components,
                closed=bool(components and all_converged and all(not c['open_lower'] and not c['open_upper'] for c in components)),
                all_nuisance_fits_converged=all_converged,
                profile_found_lower_objective=bool(best_discovered<fit['weighted_sse']-max(1e-6,fit['weighted_sse']*1e-5)),
                interval_method='Nuisance-refitted log-coordinate Brent crossing refinement; no endpoint extrapolation',
                relative_weight_scale_profiled=not settings['absolute_sigma'],
                warning='Exact macro constraint; conditional, asymptotic likelihood interval, not a guarantee of biological accuracy')
