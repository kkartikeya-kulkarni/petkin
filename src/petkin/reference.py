"""Reference-tissue methods; a biologically valid reference region is required.

Never use cerebellum as a universal reference. TSPO/PBR28 has no generally
valid receptor-free reference region; the human PBR28 benchmark is arterial.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares
from .kinetics import InputFunction,vector,predict
from .graphical import linear_integral,_regress,_weights


@dataclass(frozen=True)
class ReferenceCurve:
    time: np.ndarray
    activity: np.ndarray
    def __post_init__(self):
        checked=InputFunction(self.time,self.activity)
        object.__setattr__(self,'time',checked.time)
        object.__setattr__(self,'activity',checked.plasma)


def predict_srtm(frames,reference,parameters,dt=.025):
    """Frame-averaged SRTM: R1*Cr + (k2-R1*k2a)*(Cr convolved exp(-k2a*t))."""
    p=vector(parameters,'SRTM parameters')
    if len(p)!=3 or p[0]<=0 or p[1]<=0 or p[2]<0:raise ValueError('SRTM requires R1>0, k2>0, BPnd>=0')
    r1,k2,bp=p;ka=k2/(1+bp)
    avg=(linear_integral(reference.time,reference.activity,frames.end)-linear_integral(reference.time,reference.activity,frames.start))/frames.duration
    conv=predict(frames,InputFunction(reference.time,reference.activity),[1.,ka],'1t',dt=dt)
    return r1*avg+(k2-r1*ka)*conv


def fit_srtm(frames,reference,tissue,sigma=None,*,dt=.025,n_starts=8,seed=0):
    """Regional multistart SRTM; reference measurement uncertainty is ignored."""
    y=vector(tissue,'tissue');s=np.ones_like(y) if sigma is None else vector(sigma,'sigma')
    if y.shape!=frames.start.shape or s.shape!=y.shape or np.any(s<=0) or len(y)<=3:raise ValueError('Invalid SRTM data or weights')
    if not isinstance(n_starts,int) or n_starts<1:raise ValueError('n_starts must be positive')
    rng=np.random.default_rng(seed);lo=np.array([.01,1e-6,0.]);hi=np.array([3.,3.,15.])
    starts=[np.array([1.,.2,1.])]+[np.array([rng.uniform(.3,2),np.exp(rng.uniform(np.log(.01),np.log(.6))),rng.uniform(.1,5)]) for _ in range(n_starts-1)]
    def residual(p):return (predict_srtm(frames,reference,p,dt)-y)/s
    trials=[least_squares(residual,p,bounds=(lo,hi),x_scale='jac',max_nfev=800) for p in starts]
    fit=min(trials,key=lambda r:r.cost);p=fit.x
    covariance=np.linalg.pinv(fit.jac.T@fit.jac)*2*fit.cost/(len(y)-3)
    rank=int(np.linalg.matrix_rank(fit.jac));condition=float(np.linalg.cond(fit.jac*np.maximum(p,[.01,.01,.1])))
    bound=bool(np.any(fit.active_mask));reasons=[]
    if not fit.success:reasons.append('optimizer_not_converged')
    if bound:reasons.append('parameter_at_bound')
    if rank<3 or condition>1e8:reasons.append('ill_conditioned_SRTM')
    if np.sqrt(max(0,covariance[2,2]))/max(p[2],1e-10)>.5:reasons.append('large_local_binding_uncertainty')
    return dict(method='srtm',parameters=dict(zip(['R1','k2','BPnd'],map(float,p))),
                BPnd=float(p[2]) if not reasons else None,raw_BPnd=float(p[2]),DVR=float(1+p[2]) if not reasons else None,
                success=bool(fit.success),at_bound=bound,status='ok' if not reasons else 'review_required',
                usable_for_quantification=not reasons,qc_reasons=reasons,prediction=predict_srtm(frames,reference,p,dt).tolist(),
                weighted_sse=float(2*fit.cost),standard_errors=np.sqrt(np.maximum(0,np.diag(covariance))).tolist(),
                covariance=covariance.tolist(),jacobian_rank=rank,scaled_condition_number=condition,
                assumptions=['Reference has negligible specific binding','Comparable nondisplaceable distribution volume',
                             'One-tissue equivalent kinetics','Reference uncertainty not propagated'],
                settings=dict(dt=dt,n_starts=n_starts,seed=seed))


def logan_reference(frames,reference,tissue,tstar,k2prime,*,weights=None):
    """Reference Logan including Cr/k2prime term; k2prime is min^-1."""
    if not np.isfinite(k2prime) or k2prime<=0:raise ValueError('Positive finite k2prime required')
    y=vector(tissue,'tissue')
    if y.shape!=frames.start.shape:raise ValueError('Tissue must match frames')
    if frames.start[0]!=0 or not np.allclose(frames.start[1:],frames.end[:-1],rtol=0,atol=1e-9):raise ValueError('Contiguous tissue frames from zero required')
    w,valid=_weights(frames,weights,tstar);valid&=y>0
    if valid.sum()<3:raise ValueError('At least three positive late frames required')
    it=np.cumsum(y*frames.duration)-.5*y*frames.duration
    ir=linear_integral(reference.time,reference.activity,frames.mid)
    cr=np.interp(frames.mid,reference.time,reference.activity)
    x=(ir[valid]+cr[valid]/k2prime)/y[valid];z=it[valid]/y[valid]
    b,cov,residual,r2,condition=_regress(np.c_[x,np.ones_like(x)],z,w[valid])
    reasons=[]
    if b[0]<=0:reasons.append('nonpositive_DVR')
    if b[0]<1:reasons.append('negative_binding_potential')
    if condition>1e12:reasons.append('ill_conditioned_reference_Logan')
    return dict(method='logan_reference',DVR=float(b[0]) if not reasons else None,raw_DVR=float(b[0]),
                BPnd=float(b[0]-1) if not reasons else None,raw_BPnd=float(b[0]-1),tstar=tstar,k2prime=k2prime,
                status='ok' if not reasons else 'review_required',usable_for_quantification=not reasons,qc_reasons=reasons,
                r_squared=r2,frame_indices=np.flatnonzero(valid).tolist(),x=x.tolist(),y=z.tolist(),
                assumptions=['Valid reference region required','Fixed k2prime uncertainty ignored','No plasma or vascular correction applied'])
