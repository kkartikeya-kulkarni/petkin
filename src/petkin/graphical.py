"""Plasma Logan, Ichise MA1 and Patlak with shared exact input integration."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from .kinetics import vector


def linear_integral(time,values,query,*,tail='error'):
    """Exact integral of the piecewise-linear input, zero before input time zero."""
    q=np.asarray(query,float)
    if tail not in ('error','hold'):raise ValueError('Invalid input tail policy')
    if np.any(q>time[-1]+1e-10) and tail=='error':raise ValueError('Input does not cover integration times')
    clipped=np.clip(q,0,time[-1]);i=np.clip(np.searchsorted(time,clipped,side='right')-1,0,len(time)-2)
    delta=clipped-time[i];slope=(values[i+1]-values[i])/(time[i+1]-time[i])
    prefix=cumulative_trapezoid(values,time,initial=0)
    ans=prefix[i]+values[i]*delta+.5*slope*delta**2
    if tail=='hold':ans+=np.maximum(q-time[-1],0)*values[-1]
    return np.where(q<=0,0.,ans)


def _components(frames,tissue,inp,vb,delay,input_tail):
    y=vector(tissue,'tissue').copy()
    if y.shape!=frames.start.shape:raise ValueError('Tissue must match frame count')
    if frames.start[0]!=0 or not np.allclose(frames.start[1:],frames.end[:-1],rtol=0,atol=1e-9):
        raise ValueError('Graphical analysis requires contiguous tissue frames from zero; do not silently fill missing intervals')
    if not 0<=vb<1 or not np.isfinite(delay):raise ValueError('Invalid vb/delay')
    if inp.time[-1]<max(frames.end[-1],frames.end[-1]-delay) and input_tail=='error':raise ValueError('Input must cover scan and delay')
    if vb:
        if inp.blood is None:raise ValueError('Whole blood is required for vascular correction')
        blood=(linear_integral(inp.time,inp.blood,frames.end-delay,tail=input_tail)-linear_integral(inp.time,inp.blood,frames.start-delay,tail=input_tail))/frames.duration
        y=(y-vb*blood)/(1-vb)
    it=np.cumsum(y*frames.duration)-.5*y*frames.duration
    ip=linear_integral(inp.time,inp.plasma,frames.mid-delay,tail=input_tail)
    return y,it,ip


def _regress(design,y,weights):
    root=np.sqrt(weights);a=design*root[:,None];b=y*root
    if np.linalg.matrix_rank(a)<design.shape[1]:raise ValueError('Degenerate graphical design')
    beta=np.linalg.lstsq(a,b,rcond=None)[0];residual=y-design@beta
    dof=len(y)-len(beta)
    covariance=np.linalg.pinv(a.T@a)*np.sum(weights*residual**2)/dof
    total=np.sum(weights*(y-np.average(y,weights=weights))**2)
    return beta,covariance,residual,float(1-np.sum(weights*residual**2)/total) if total>0 else None,float(np.linalg.cond(a))


def _weights(frames,weights,tstar):
    if not np.isfinite(tstar) or tstar<0:raise ValueError('tstar must be finite and nonnegative')
    w=np.ones(len(frames.start)) if weights is None else vector(weights,'weights')
    if w.shape!=frames.start.shape or np.any(w<0):raise ValueError('Weights must be nonnegative and match frames')
    return w,(frames.mid>=tstar)&(w>0)


def logan(frames,tissue,input_function,tstar,vb=0.,*,weights=None,delay=0.,input_tail='error'):
    """Plasma Logan; errors on both axes can bias OLS. tstar is never auto-tuned."""
    ct,it,ip=_components(frames,tissue,input_function,vb,delay,input_tail)
    w,valid=_weights(frames,weights,tstar);valid&=ct>0
    if valid.sum()<3:raise ValueError('At least three positive late frames required')
    x,z=ip[valid]/ct[valid],it[valid]/ct[valid]
    b,cov,residual,r2,condition=_regress(np.c_[x,np.ones_like(x)],z,w[valid])
    reasons=[]
    if b[0]<=0:reasons.append('nonpositive_VT')
    if condition>1e12:reasons.append('ill_conditioned_graphical_fit')
    return dict(method='logan',VT=float(b[0]) if not reasons else None,raw_VT=float(b[0]),
                intercept=float(b[1]),tstar=float(tstar),vb=vb,delay=delay,n=int(valid.sum()),
                status='ok' if not reasons else 'review_required',qc_reasons=reasons,
                usable_for_quantification=not reasons,r_squared=r2,condition_number=condition,
                slope_standard_error=float(np.sqrt(max(0,cov[0,0]))),
                uncertainty_caveat='Conditional linear error only; correlated errors in both Logan axes are ignored',
                x=x.tolist(),y=z.tolist(),residuals=residual.tolist(),frame_indices=np.flatnonzero(valid).tolist())


def ma1(frames,tissue,input_function,tstar,vb=0.,*,weights=None,delay=0.,input_tail='error'):
    """Ichise MA1: Ct = beta_plasma*integral(Cp)+beta_tissue*integral(Ct).

    VT = -beta_plasma/beta_tissue. Reduces a source of Logan denominator noise;
    neither method is unbiased for every tracer, scan duration or tstar.
    """
    ct,it,ip=_components(frames,tissue,input_function,vb,delay,input_tail)
    w,valid=_weights(frames,weights,tstar)
    if valid.sum()<4:raise ValueError('At least four late frames required for MA1')
    design=np.c_[ip[valid],it[valid]]
    b,cov,residual,r2,condition=_regress(design,ct[valid],w[valid])
    raw=float(-b[0]/b[1]) if abs(b[1])>1e-15 else None
    reasons=[]
    if raw is None or raw<=0:reasons.append('nonpositive_or_undefined_VT')
    if b[1]+1.96*np.sqrt(max(0,cov[1,1]))>=0:reasons.append('unstable_MA1_denominator')
    if condition>1e12:reasons.append('ill_conditioned_graphical_fit')
    se=None
    if raw is not None:
        gradient=np.array([-1/b[1],b[0]/b[1]**2]);se=float(np.sqrt(max(0,gradient@cov@gradient)))
    return dict(method='ma1',VT=raw if not reasons else None,raw_VT=raw,coefficients=b.tolist(),
                tstar=float(tstar),vb=vb,delay=delay,n=int(valid.sum()),status='ok' if not reasons else 'review_required',
                usable_for_quantification=not reasons,qc_reasons=reasons,r_squared=r2,condition_number=condition,
                macro_standard_error=se,uncertainty_caveat='Conditional regression approximation; input and integral uncertainty ignored',
                frame_indices=np.flatnonzero(valid).tolist(),prediction=(design@b).tolist(),residuals=residual.tolist())


def patlak(frames,tissue,input_function,tstar,vb=0.,*,weights=None,delay=0.,input_tail='error'):
    """Plasma Patlak for effectively irreversible trapping after chosen tstar."""
    ct,it,ip=_components(frames,tissue,input_function,vb,delay,input_tail)
    w,valid=_weights(frames,weights,tstar)
    cp=np.interp(frames.mid-delay,input_function.time,input_function.plasma,left=0.)
    valid&=cp>0
    if valid.sum()<3:raise ValueError('At least three positive plasma late frames required')
    x,z=ip[valid]/cp[valid],ct[valid]/cp[valid]
    b,cov,residual,r2,condition=_regress(np.c_[x,np.ones_like(x)],z,w[valid])
    reasons=[]
    if b[0]<=0:reasons.append('nonpositive_Ki')
    if condition>1e12:reasons.append('ill_conditioned_graphical_fit')
    return dict(method='patlak',Ki=float(b[0]) if not reasons else None,raw_Ki=float(b[0]),
                intercept=float(b[1]),tstar=float(tstar),vb=vb,delay=delay,n=int(valid.sum()),
                status='ok' if not reasons else 'review_required',usable_for_quantification=not reasons,
                qc_reasons=reasons,r_squared=r2,slope_standard_error=float(np.sqrt(max(0,cov[0,0]))),
                x=x.tolist(),y=z.tolist(),residuals=residual.tolist(),frame_indices=np.flatnonzero(valid).tolist())


def graphical_stability(frames,tissue,input_function,tstars,method='logan',**options):
    """Report every requested tstar, without picking the prettiest result."""
    methods={'logan':logan,'ma1':ma1,'patlak':patlak}
    if method not in methods:raise ValueError('Unknown graphical method')
    results=[]
    for t in vector(tstars,'tstars'):
        try:results.append(methods[method](frames,tissue,input_function,float(t),**options))
        except ValueError as e:results.append(dict(tstar=float(t),status='failed',reason=str(e)))
    key='Ki' if method=='patlak' else 'VT'
    values=[r[key] for r in results if r.get(key) is not None]
    span=float(np.ptp(values)/np.median(values)) if values and np.median(values)>0 else None
    return dict(method=method,results=results,relative_span=span,selection='No automatic tstar selection')
