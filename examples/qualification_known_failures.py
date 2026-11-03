"""Regression audit of the two exposed v0.2 misses; not confirmatory evidence."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from petkin import Frames,InputFunction,analyze_research_tac
from petkin.validation import REQUIRED_ASSUMPTIONS
source,out=Path(sys.argv[1]),Path(sys.argv[2]);old=json.loads(source.read_text());records=[]
t=np.linspace(0,120,12001);cp=8*t*np.exp(-t/.7)+1.5*np.exp(-t/25);inp=InputFunction(t,cp,.8*cp)
for case in old['cases']:
    if not case['large_error_unflagged']:continue
    k4=case['k4_truth'];length=case['scan_minutes'];edges=np.arange(0,length+1.);f=Frames(edges[:-1],edges[1:])
    def rhs(time,y):
        plasma=np.interp(time,t,cp)
        return [.2*plasma-(.25+.1)*y[0]+k4*y[1],.1*y[0]-k4*y[1],.95*(y[0]+y[1])+.05*.8*plasma]
    sol=solve_ivp(rhs,(0,length),[0,0,0],rtol=1e-10,atol=1e-12,max_step=.05,dense_output=True)
    truth=(sol.sol(f.end)[2]-sol.sol(f.start)[2])/f.duration
    exposure=case['sensitivity']*f.duration*np.exp(-np.log(2)*f.mid/20.33);sigma=np.sqrt((truth+.02*truth.max())/exposure)
    observed=np.array(case['fit']['prediction'])-np.array(case['fit']['residuals'])
    result=analyze_research_tac(f,inp,observed,sigma=sigma,vb=.05,absolute_sigma=True,
                               assumptions={name:True for name in REQUIRED_ASSUMPTIONS})
    records.append(dict(truth_k4=k4,old_relative_error=case['raw_relative_error'],
                        new_relative_error=result['raw_VT']/case['truth_VT']-1,result=result,role='Exposed regression case'))
    print(k4,result['research_screen']['passed'],result['research_screen']['reasons'],flush=True)
(out/'known_failure_regression.json').write_text(json.dumps(records,indent=2))
