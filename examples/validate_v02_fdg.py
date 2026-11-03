"""Exploratory FDG model/Patlak agreement on the same three real subjects."""
from pathlib import Path
import json,sys
import numpy as np
from petkin import Frames,InputFunction,fit_compartment,graphical_stability
root,out=Path(sys.argv[1]),Path(sys.argv[2]);rows=[]
for name in ['s028','s029','s030']:
    folder=list(root.glob(f'**/{name}/tac.csv'))
    if len(folder)!=1:raise ValueError('Need one frozen regional TAC per subject')
    folder=folder[0].parent;t=np.genfromtxt(folder/'tac.csv',delimiter=',',names=True);b=np.genfromtxt(folder/'input.csv',delimiter=',',names=True)
    frames=Frames(t['start'],t['end']);inp=InputFunction(b['time'],b['plasma']);valid=t['activity']!=0
    fit=fit_compartment(Frames(frames.start[valid],frames.end[valid]),inp,t['activity'][valid],model='2t_irreversible')
    graphical=graphical_stability(frames,t['activity'],inp,[15,20,25],method='patlak')
    rows.append(dict(subject=name,fit=fit,patlak=graphical,
                     limitations=['No whole blood: vb fixed at zero','Midpoint-sampled input and boundary assumptions',
                                  'No known Ki; method agreement is not accuracy','Input decay convention not independently verified']))
    print(name,fit['Ki'],[r.get('Ki') for r in graphical['results']],fit['residual_diagnostics'],flush=True)
(out/'fdg_v02.json').write_text(json.dumps(rows,indent=2))
