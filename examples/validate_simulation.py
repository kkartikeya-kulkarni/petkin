"""Reproducible synthetic held-out cohort, independent ODE ground truth."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from petkin import Frames, InputFunction, fit_compartment, logan


def run(out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261102)
    t=np.linspace(0,120,12001)
    edges=np.r_[np.arange(0,5,.25),np.arange(5,20,1),np.arange(20,121,2)]
    frames=Frames(edges[:-1],edges[1:]); rows=[]
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for i in range(12):
        p=np.array([rng.uniform(.12,.3),rng.uniform(.15,.35),rng.uniform(.05,.12),rng.uniform(.05,.1)])
        cp=8*t*np.exp(-t/.7)+rng.uniform(1,2)*np.exp(-t/25)
        inp=InputFunction(t,cp)
        def rhs(time,y):
            return [p[0]*np.interp(time,t,cp)-(p[1]+p[2])*y[0]+p[3]*y[1],p[2]*y[0]-p[3]*y[1],y[0]+y[1]]
        sol=solve_ivp(rhs,(0,120),[0,0,0],rtol=1e-9,atol=1e-11,max_step=.03,dense_output=True)
        truth=(sol.sol(frames.end)[2]-sol.sol(frames.start)[2])/frames.duration
        sigma=.01*np.maximum(truth,.05*truth.max())
        observed=truth+rng.normal(0,sigma)
        fit=fit_compartment(frames,inp,observed,sigma=sigma)
        graphical=logan(frames,observed,inp,tstar=50)
        vt=p[0]/p[1]*(1+p[2]/p[3])
        estimated=np.array(list(fit['parameters'].values()))
        rows.append(dict(subject=i,truth_parameters=p.tolist(),fit=fit,logan=graphical,truth_VT=vt,
                         parameter_relative_error=((estimated-p)/p).tolist(),
                         fit_VT_relative_error=fit['VT']/vt-1,logan_VT_relative_error=graphical['VT']/vt-1))
        np.savetxt(out/f'synthetic_{i:02d}_tac.csv',np.c_[frames.start,frames.end,observed,sigma],delimiter=',',header='start,end,activity,sigma',comments='')
        np.savetxt(out/f'synthetic_{i:02d}_input.csv',np.c_[t,cp],delimiter=',',header='time,plasma',comments='')
        if i==0:
            np.savetxt(out/'input.csv',np.c_[t,cp],delimiter=',',header='time,plasma',comments='')
            axes[0].errorbar(frames.mid,observed,yerr=sigma,fmt='.',label='Simulated observations')
            axes[0].plot(frames.mid,fit['prediction'],label='2-tissue fit'); axes[0].set(xlabel='Time (min)',ylabel='Activity (arbitrary units)'); axes[0].legend()
            x=np.array(graphical['x']); axes[1].scatter(x,graphical['y'],s=12)
            axes[1].plot(x,graphical['VT']*x+graphical['intercept']); axes[1].set(xlabel='Integrated plasma / tissue (min)',ylabel='Integrated tissue / tissue (min)',title='Logan, t* = 50 min')
    errors=np.array([r['parameter_relative_error'] for r in rows])
    summary=dict(cohort='12 synthetic subjects; fixed seed 20261102; no clinical data',
                 noise='Independent Gaussian SD = 1% of activity with 5%-peak floor; not a scanner noise model',
                 simulator='Adaptive ODE solver with integrated observation state, independent of fitting convolution',
                 parameter_median_absolute_relative_error=np.median(abs(errors),axis=0).tolist(),
                 parameter_max_absolute_relative_error=np.max(abs(errors),axis=0).tolist(),
                 fit_VT_median_absolute_relative_error=float(np.median([abs(r['fit_VT_relative_error']) for r in rows])),
                 logan_VT_median_absolute_relative_error=float(np.median([abs(r['logan_VT_relative_error']) for r in rows])),
                 all_fits_converged=all(r['fit']['success'] for r in rows),
                 acceptance='All converge; median parameter errors < 10%; median fit and Logan VT errors < 5%')
    summary['passed']=bool(summary['all_fits_converged'] and np.all(np.median(abs(errors),axis=0)<.1) and summary['fit_VT_median_absolute_relative_error']<.05 and summary['logan_VT_median_absolute_relative_error']<.05)
    (out/'synthetic_results.json').write_text(json.dumps(dict(summary=summary,subjects=rows),indent=2))
    fig.tight_layout(); fig.savefig(out/'synthetic_validation.png',dpi=150); plt.close(fig)
    print(json.dumps(summary,indent=2))
    if not summary['passed']: raise SystemExit('Validation thresholds failed')

if __name__=='__main__':
    import sys
    run(sys.argv[1] if len(sys.argv)>1 else 'validation')
