"""Conditional uncertainty audits of a stable scan and the known difficult scan."""
from pathlib import Path
import json,sys
import numpy as np
from petkin import Frames,InputFunction,profile_likelihood,parametric_bootstrap


def run(data,out):
    data,out=Path(data),Path(out);records=json.loads((out/'real_v02.json').read_text());audits=[]
    for name,variant in [('cgyu_1','fixed'),('rbqc_1','fixed'),('rbqc_1','fitted_vascular_delay')]:
        row=next(r for r in records if r['scan']==name);fit=row['models'][variant]['full_fit']
        tac=np.genfromtxt(data/name/'tacs.csv',delimiter=',',names=True);tac=tac[tac['Duration']>0]
        f=Frames(tac['StartTime']/60,(tac['StartTime']+tac['Duration'])/60);b=np.genfromtxt(data/name/'input.csv',delimiter=',',names=True)
        time,cp,cb=b['Time'],b['AIF'],b['Blood']
        if time[-1]<f.end[-1]:time=np.r_[time,f.end[-1]];cp=np.r_[cp,cp[-1]];cb=np.r_[cb,cb[-1]]
        inp=InputFunction(time,cp,cb);valid=tac['Weights']>0;f=Frames(f.start[valid],f.end[valid]);y=tac['FC'][valid];sigma=1/np.sqrt(tac['Weights'][valid])
        profile=profile_likelihood(f,inp,y,fit,'k4',np.geomspace(1e-6,3.,51),sigma=sigma)
        bootstrap=parametric_bootstrap(f,inp,y,fit,sigma=sigma,n=100,seed=20261103,n_starts=3)
        audits.append(dict(scan=name,variant=variant,profile_k4=profile,bootstrap=bootstrap))
        (out/'uncertainty_audit.json').write_text(json.dumps(audits,indent=2))
        print(name,variant,'profile_closed',profile['closed'],'bootstrap invalid',bootstrap['unusable_fraction'],bootstrap['raw_percentile_95_interval'],flush=True)
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
