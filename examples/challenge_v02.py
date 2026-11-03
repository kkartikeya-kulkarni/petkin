"""Hard identifiability challenge: slow dissociation, short scans, count noise.

Synthetic counts are not a full reconstruction/scanner model. No QC threshold
is changed using these cases. All raw errors and accepted failures are retained.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import solve_ivp
from petkin import Frames,InputFunction,fit_compartment


def run(out):
    rng=np.random.default_rng(1742);t=np.linspace(0,120,12001)
    cp=8*t*np.exp(-t/.7)+1.5*np.exp(-t/25);inp=InputFunction(t,cp,.8*cp);rows=[]
    for k4 in [.0001,.001,.01,.05]:
        def rhs(time,y):
            plasma=np.interp(time,t,cp)
            return [.2*plasma-(.25+.1)*y[0]+k4*y[1],.1*y[0]-k4*y[1],.95*(y[0]+y[1])+.05*.8*plasma]
        sol=solve_ivp(rhs,(0,120),[0,0,0],rtol=1e-10,atol=1e-12,max_step=.05,dense_output=True)
        for length in [30,60,120]:
            edges=np.arange(0,length+.1,1.);f=Frames(edges[:-1],edges[1:])
            truth=(sol.sol(f.end)[2]-sol.sol(f.start)[2])/f.duration
            true_vt=.2/.25*(1+.1/k4)
            for sensitivity in [100,1000]:
                exposure=sensitivity*f.duration*np.exp(-np.log(2)*f.mid/20.33)
                background=.02*truth.max()
                counts=rng.poisson((truth+background)*exposure)
                y=counts/exposure-background;sigma=np.sqrt((truth+background)/exposure)
                fit=fit_compartment(f,inp,y,sigma=sigma,vb=.05,n_starts=6,absolute_sigma=True)
                error=fit['raw_VT']/true_vt-1
                rows.append(dict(k4_truth=k4,scan_minutes=length,sensitivity=sensitivity,truth_VT=true_vt,
                                 raw_relative_error=error,fit=fit,
                                 large_error_unflagged=bool(abs(error)>.5 and fit['usable_for_quantification'])))
    summary=dict(cases=len(rows),flagged=sum(not r['fit']['usable_for_quantification'] for r in rows),
                 raw_errors_over_50_percent=sum(abs(r['raw_relative_error'])>.5 for r in rows),
                 large_errors_not_flagged=sum(r['large_error_unflagged'] for r in rows),
                 median_absolute_raw_error=float(np.median([abs(r['raw_relative_error']) for r in rows])),
                 protocol='4 k4 values x 3 scan lengths x 2 count scales; seed 1742; C11-like decay; Poisson counts plus background subtraction')
    Path(out,'identifiability_challenge.json').write_text(json.dumps(dict(summary=summary,cases=rows),indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':run(sys.argv[1])
