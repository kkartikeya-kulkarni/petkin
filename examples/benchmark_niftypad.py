"""Independent ODE SRTM cases, compared with installed NiftyPAD 2.0.0.

NiftyPAD uses SRTM basis (1000 rates), second-grid discrete convolution; petkin uses minute-grid
trapezoidal convolution and explicit frame means. Small discretization
agreement differences are expected. This is synthetic, not a reference-region
validation on human PBR28 (where a reference model would be inappropriate).
"""
from pathlib import Path
import hashlib,inspect,json,sys,importlib.metadata
import numpy as np
from scipy.integrate import solve_ivp
import niftypad.models as nifty
from petkin import Frames,ReferenceCurve,fit_srtm,logan_reference


def run(out):
    out=Path(out);t=np.arange(3601)/60;cr=5*(np.exp(-t*.025)-np.exp(-t*.3))
    reference=ReferenceCurve(t,cr);edges=np.arange(0,60.1,.5);frames=Frames(edges[:-1],edges[1:])
    dt=np.array([frames.start*60,frames.end*60]).astype(int);rows=[]
    for r1 in [.7,1.2]:
        for k2 in [.1,.25]:
            for bp in [.5,1.5,3.]:
                ka=k2/(1+bp)
                def rhs(time,y):
                    c=5*(np.exp(-time*.025)-np.exp(-time*.3))
                    return [c-ka*y[0],r1*c+(k2-r1*ka)*y[0]]
                sol=solve_ivp(rhs,(0,60),[0,0],rtol=1e-10,atol=1e-12,max_step=.05,dense_output=True)
                y=(sol.sol(frames.end)[1]-sol.sol(frames.start)[1])/frames.duration
                own=fit_srtm(frames,reference,y,n_starts=4,dt=1/600)
                other=nifty.srtmb(y.copy(),dt,cr[:-1],[5e-5,.02],1000,None)
                own_graph=logan_reference(frames,reference,y,30,k2/r1)
                other_graph=nifty.logan_ref_k2p(y,dt,cr[:-1],k2/r1/60,None,1800,3600,False)
                row=dict(truth=dict(R1=r1,k2=k2,BPnd=bp),petkin=own,
                         niftypad=dict(R1=float(other['r1']),k2=float(other['k2']*60),BPnd=float(other['bp'])),
                         BPnd_relative_difference=own['raw_BPnd']/other['bp']-1,
                         petkin_BPnd_truth_error=own['raw_BPnd']/bp-1,
                         petkin_reference_Logan_DVR=own_graph['raw_DVR'],niftypad_reference_Logan_DVR=float(other_graph['bp']+1))
                rows.append(row);print(r1,k2,bp,row['BPnd_relative_difference'],flush=True)
    source=Path(inspect.getfile(nifty))
    result=dict(niftypad_method='srtmb basis; 1000 log-spaced rates [5e-5,.02] per second',niftypad_version=importlib.metadata.version('niftypad'),niftypad_models_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                cases=rows,summary=dict(cases=len(rows),median_absolute_BPnd_difference=float(np.median([abs(r['BPnd_relative_difference']) for r in rows])),
                                       maximum_absolute_BPnd_difference=float(max(abs(r['BPnd_relative_difference']) for r in rows))))
    (out/'niftypad_agreement.json').write_text(json.dumps(result,indent=2));print(json.dumps(result['summary'],indent=2))
if __name__=='__main__':run(sys.argv[1])
