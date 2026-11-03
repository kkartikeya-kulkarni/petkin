"""One prospectively frozen research qualification run; failures are retained."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import beta
from petkin import Frames,InputFunction,analyze_research_tac
from petkin.validation import REQUIRED_ASSUMPTIONS


def run(out):
    out=Path(out);protocol=json.loads((out/'PROTOCOL.json').read_text())
    expected=(out/'PROTOCOL.sha256').read_text().split()[0]
    if hashlib.sha256((out/'PROTOCOL.json').read_bytes()).hexdigest()!=expected:raise ValueError('Frozen protocol changed')
    rng=np.random.default_rng(protocol['holdout']['seed']);t=np.linspace(0,120,12001);cases=[]
    def make_case(model,p,length,input_spec,kind,index):
        cp=input_spec['bolus']*t*np.exp(-t/input_spec['fast_tau'])+input_spec['tail']*np.exp(-t/input_spec['slow_tau'])
        inp=InputFunction(t,cp,.8*cp);vb=.05;delay=.1 if kind=='routine' else 0.
        k1,k2=p[:2];k3=p[2] if len(p)>2 else 0;k4=p[3] if len(p)>3 else 0
        def rhs(time,y):
            plasma=np.interp(time-delay,t,cp,left=0.)
            return [k1*plasma-(k2+k3)*y[0]+k4*y[1],k3*y[0]-k4*y[1],.95*(y[0]+y[1])+.05*.8*plasma]
        if kind=='routine':edges=np.r_[np.arange(0,5,.25),np.arange(5,20,1.),np.arange(20,length+1,2.)]
        else:edges=np.arange(0,length+1,1.)
        f=Frames(edges[:-1],edges[1:]);sol=solve_ivp(rhs,(0,length),[0,0,0],rtol=1e-9,atol=1e-11,max_step=.05,dense_output=True)
        truth=(sol.sol(f.end)[2]-sol.sol(f.start)[2])/f.duration
        spec={}
        if kind=='routine':
            noise=float(rng.uniform(.02,.08));sigma=noise*np.maximum(truth,.05*truth.max());observed=truth+rng.normal(0,sigma);spec['noise']=noise
        else:
            sensitivity=int(rng.choice([100,1000]));exposure=sensitivity*f.duration*np.exp(-np.log(2)*f.mid/20.33);background=.02*truth.max()
            observed=rng.poisson((truth+background)*exposure)/exposure-background;sigma=np.sqrt((truth+background)/exposure);spec['count_sensitivity']=sensitivity
        target=k1*k3/(k2+k3) if model=='2t_irreversible' else k1/k2*(1+k3/k4 if model=='2t' else 1)
        cases.append(dict(id=f'{kind}-{model}-{index:03d}',kind=kind,model=model,parameters=list(map(float,p)),
                          length=float(length),input_spec=input_spec,vb=vb,delay=delay,truth_macro=float(target),
                          start=f.start.tolist(),end=f.end.tolist(),tissue=observed.tolist(),sigma=sigma.tolist(),**spec))
    for model in ['1t','2t','2t_irreversible']:
        for i in range(protocol['holdout']['routine_cases_per_model']):
            p=[rng.uniform(.1,.3),rng.uniform(.12,.35)]
            if model!='1t':p.append(rng.uniform(.03,.15))
            if model=='2t':p.append(rng.uniform(.03,.12))
            spec=dict(bolus=8.,fast_tau=float(rng.uniform(.4,1.2)),tail=float(rng.uniform(.8,2.)),slow_tau=float(rng.uniform(15,35)))
            make_case(model,p,120,spec,'routine',i)
    for i in range(protocol['holdout']['stress_cases']):
        k4=float(np.exp(rng.uniform(np.log(.0001),np.log(.05))));length=int(rng.choice([30,60,120]))
        make_case('2t',[.2,.25,.1,k4],length,dict(bolus=8.,fast_tau=.7,tail=1.5,slow_tau=25.),'stress',i)
    data_path=out/'holdout_data.json'
    if data_path.exists():raise ValueError('Refusing to overwrite the frozen holdout dataset; choose a separate protocol directory')
    data_path.write_text(json.dumps(cases,indent=2));data_sha=hashlib.sha256(data_path.read_bytes()).hexdigest()
    source_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (out/'holdout_manifest.json').write_text(json.dumps(dict(protocol_sha256=expected,data_sha256=data_sha,generator_sha256=source_sha,
                                                           locally_frozen_not_public_preregistration=True),indent=2))
    records=[];assumptions={name:True for name in REQUIRED_ASSUMPTIONS}
    for number,case in enumerate(cases):
        spec=case['input_spec'];cp=spec['bolus']*t*np.exp(-t/spec['fast_tau'])+spec['tail']*np.exp(-t/spec['slow_tau'])
        inp=InputFunction(t,cp,.8*cp);f=Frames(case['start'],case['end'])
        result=analyze_research_tac(f,inp,case['tissue'],case['model'],case['sigma'],assumptions=assumptions,
                                    vb=case['vb'],delay=case['delay'],absolute_sigma=True,n_starts=8)
        key='Ki' if case['model']=='2t_irreversible' else 'VT';error=result['raw_'+key]/case['truth_macro']-1
        records.append(dict(id=case['id'],kind=case['kind'],model=case['model'],truth_macro=case['truth_macro'],
                            raw_relative_error=float(error),released=result['research_screen']['passed'],result=result))
        (out/'holdout_results.json').write_text(json.dumps(records,indent=2))
        print(number+1,'/',len(cases),case['id'],'released',result['research_screen']['passed'],flush=True)
    summaries={}
    for kind in ['routine','stress']:
        group=[r for r in records if r['kind']==kind];accepted=[r for r in group if r['released']]
        errors=np.array([abs(r['raw_relative_error']) for r in accepted]);large=int(np.sum(errors>.5));n=len(accepted)
        upper=1. if n==0 or large==n else float(beta.ppf(.95,large+1,n-large))
        summaries[kind]=dict(total=len(group),released=n,retention=n/len(group),large_errors_released=large,
                             median_absolute_error=float(np.median(errors)) if n else None,
                             q95_absolute_error=float(np.quantile(errors,.95)) if n else None,
                             one_sided_95_large_error_risk_upper=upper)
    (out/'holdout_summary.json').write_text(json.dumps(summaries,indent=2));print(json.dumps(summaries,indent=2))
if __name__=='__main__':run(sys.argv[1])
