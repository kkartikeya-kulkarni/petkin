"""Independent ODE Monte Carlo across model families and three noise levels."""
import json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from petkin import Frames,InputFunction,fit_compartment,logan,ma1,patlak


def run(out,n_per_model=8):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(20261103)
    t=np.linspace(0,120,12001);edges=np.r_[np.arange(0,5,.25),np.arange(5,20,1),np.arange(20,121,2)]
    f=Frames(edges[:-1],edges[1:]);rows=[]
    for model in ['1t','2t','2t_irreversible']:
        for person in range(n_per_model):
            p=[rng.uniform(.1,.3),rng.uniform(.12,.35)]
            if model!='1t':p.append(rng.uniform(.03,.15))
            if model=='2t':p.append(rng.uniform(.025,.12))
            k1,k2=p[:2];k3=p[2] if len(p)>2 else 0;k4=p[3] if len(p)>3 else 0
            cp=8*t*np.exp(-t/rng.uniform(.4,1.2))+rng.uniform(.8,2)*np.exp(-t/rng.uniform(15,35))
            inp=InputFunction(t,cp,.8*cp);vb=.05;delay=.1
            def rhs(time,y):
                plasma=np.interp(time-delay,t,cp,left=0.)
                return [k1*plasma-(k2+k3)*y[0]+k4*y[1],k3*y[0]-k4*y[1],(1-vb)*(y[0]+y[1])+vb*.8*plasma]
            sol=solve_ivp(rhs,(0,120),[0,0,0],rtol=1e-9,atol=1e-11,max_step=.05,dense_output=True)
            truth=(sol.sol(f.end)[2]-sol.sol(f.start)[2])/f.duration
            target=k1*k3/(k2+k3) if model=='2t_irreversible' else k1/k2*(1+k3/k4 if model=='2t' else 1)
            for noise in [.01,.05,.15]:
                sigma=noise*np.maximum(truth,.05*truth.max());y=truth+rng.normal(0,sigma)
                fit=fit_compartment(f,inp,y,model,sigma,vb=vb,delay=delay,absolute_sigma=True,n_starts=6)
                key='Ki' if model=='2t_irreversible' else 'VT'
                row=dict(model=model,case=person,noise=noise,truth_parameters=p,truth_macro=target,
                         raw_relative_error=fit['raw_'+key]/target-1,fit=fit,
                         local_interval_covers_truth=fit['macro_local_95_interval'][0]<=target<=fit['macro_local_95_interval'][1])
                if model=='2t_irreversible':row['patlak']=patlak(f,y,inp,60,vb=vb,delay=delay)
                else:
                    row['logan']=logan(f,y,inp,50,vb=vb,delay=delay)
                    row['ma1']=ma1(f,y,inp,50,vb=vb,delay=delay)
                rows.append(row)
            print(model,person+1,'/',n_per_model,flush=True)
    summary=[]
    for model in ['1t','2t','2t_irreversible']:
        for noise in [.01,.05,.15]:
            group=[r for r in rows if r['model']==model and r['noise']==noise]
            errors=np.array([r['raw_relative_error'] for r in group]);accepted=[r for r in group if r['fit']['usable_for_quantification']]
            stats=dict(model=model,noise=noise,n=len(group),median_absolute_macro_error=float(np.median(abs(errors))),
                       median_macro_bias=float(np.median(errors)),max_absolute_macro_error=float(max(abs(errors))),
                       usable=len(accepted),converged=sum(r['fit']['success'] for r in group),
                       local_95_coverage=float(np.mean([r['local_interval_covers_truth'] for r in group])))
            for method in (['patlak'] if model=='2t_irreversible' else ['logan','ma1']):
                key='raw_Ki' if method=='patlak' else 'raw_VT'
                stats[method+'_median_absolute_error']=float(np.median([abs(r[method][key]/r['truth_macro']-1) for r in group if r[method].get(key) is not None]))
            summary.append(stats)
    result=dict(protocol='24 independent ODE TACs, 3 known Gaussian noise levels = 72 fits; seed 20261103; known fixed vb/delay',
                limitations='Small engineering Monte Carlo, not a scanner/reconstruction simulation or calibrated clinical coverage study',summary=summary,cases=rows)
    (out/'simulation_v02.json').write_text(json.dumps(result,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':run(sys.argv[1])
