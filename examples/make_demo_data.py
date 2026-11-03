"""Write a reproducible simulated PET dataset for the README quick start."""
import argparse,json
from pathlib import Path
import numpy as np
from petkin import Frames,InputFunction,predict

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir',type=Path)
    args=parser.parse_args()
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    for name in ('tac.csv','input.csv','truth.json'):
        if (out/name).exists():raise FileExistsError(f'Refusing to overwrite {out/name}')
    edges=np.r_[np.arange(0,5,.5),np.arange(5,20,1.),np.arange(20,61,2.)]
    frames=Frames(edges[:-1],edges[1:]);time=np.linspace(0,60,6001)
    plasma=4000*(8*time*np.exp(-time/.7)+.8*(1-np.exp(-time/.1))*np.exp(-time/30))
    whole_blood=.8*plasma;inp=InputFunction(time,plasma,whole_blood)
    parameters=[.15,.2,.08,.05];vb=.05
    expected=predict(frames,inp,parameters,model='2t',vb=vb,dt=.01)
    sigma=.02*np.maximum(expected,.05*expected.max())
    activity=expected+np.random.default_rng(12007).normal(0,sigma)
    np.savetxt(out/'tac.csv',np.c_[frames.start,frames.end,activity,sigma],delimiter=',',header='start,end,activity,sigma',comments='')
    np.savetxt(out/'input.csv',np.c_[time,plasma,whole_blood],delimiter=',',header='time,plasma,blood',comments='')
    truth=dict(role='Simulated demonstration only',model='2t',parameters=dict(zip(['K1','k2','k3','k4'],parameters)),VT=parameters[0]/parameters[1]*(1+parameters[2]/parameters[3]),vb=vb,delay_minutes=0,seed=12007,time_units='minutes',concentration_units='Bq/mL')
    (out/'truth.json').write_text(json.dumps(truth,indent=2))
    print(f'Wrote simulated TAC, input and truth to {out}')

if __name__=='__main__':main()
