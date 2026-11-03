"""Prespecified fixed versus fitted nuisance models on the exposed PBR28 cohort.

Development benchmark, not a new external held-out cohort. No settings are
changed after viewing scores. CV scores are on regional tissue observations.
"""
import csv,json,sys
from pathlib import Path
import numpy as np
from petkin import Frames,InputFunction,fit_compartment,predict,logan,ma1,graphical_stability


def run(data,out):
    data,out=Path(data),Path(out);out.mkdir(parents=True,exist_ok=True)
    records=[]
    for item in csv.DictReader((data/'scans.csv').open()):
        name=item['scan'];tac=np.genfromtxt(data/name/'tacs.csv',delimiter=',',names=True)
        tac=tac[tac['Duration']>0];f=Frames(tac['StartTime']/60,(tac['StartTime']+tac['Duration'])/60)
        b=np.genfromtxt(data/name/'input.csv',delimiter=',',names=True);time,cp,cb=b['Time'],b['AIF'],b['Blood']
        if time[-1]<f.end[-1]:time=np.r_[time,f.end[-1]];cp=np.r_[cp,cp[-1]];cb=np.r_[cb,cb[-1]]
        inp=InputFunction(time,cp,cb);y=tac['FC'];eligible=tac['Weights']>0;idx=np.arange(len(y))
        models={}
        for variant,options in [('fixed',{}),('fitted_vascular_delay',dict(fit_vb=True,fit_delay=True,input_tail='hold'))]:
            fit=fit_compartment(Frames(f.start[eligible],f.end[eligible]),inp,y[eligible],sigma=1/np.sqrt(tac['Weights'][eligible]),vb=.05,**options)
            cv=np.full(len(y),np.nan);folds=[]
            for fold in range(5):
                held=(idx%5==fold)&eligible;train=(~held)&eligible
                fitfold=fit_compartment(Frames(f.start[train],f.end[train]),inp,y[train],sigma=1/np.sqrt(tac['Weights'][train]),vb=.05,n_starts=4,**options)
                p=predict(f,inp,list(fitfold['parameters'].values()),vb=fitfold['vb'],delay=fitfold['delay'],input_tail=options.get('input_tail','error'))
                cv[held]=p[held]
                folds.append(dict(fold=fold,usable=fitfold['usable_for_quantification'],qc_reasons=fitfold['qc_reasons'],
                                  raw_VT=fitfold['raw_VT'],held_out_indices=idx[held].tolist(),
                                  nrmse=float(np.sqrt(np.mean((p[held]-y[held])**2))/np.mean(y[held]))))
            models[variant]=dict(full_fit=fit,folds=folds,cross_validated_nrmse=float(np.sqrt(np.mean((cv[eligible]-y[eligible])**2))/np.mean(y[eligible])))
        # Make the missing initial tissue interval explicit, matching the v0.1 audit.
        gap=f.start[0]
        if gap>0:gf=Frames(np.r_[0,f.start],np.r_[gap,f.end]);gy=np.r_[0,y]
        else:gf,gy=f,y
        graphical=dict(logan=graphical_stability(gf,gy,inp,[30,40,50],method='logan',vb=.05),
                       ma1=graphical_stability(gf,gy,inp,[30,40,50],method='ma1',vb=.05))
        records.append(dict(scan=name,subject=item['subject'],models=models,graphical=graphical,
                            initial_zero_interval_minutes=float(gap),role='Previously exposed development benchmark'))
        (out/'real_v02.json').write_text(json.dumps(records,indent=2))
        print(name,{k:(v['cross_validated_nrmse'],v['full_fit']['VT'],v['full_fit']['qc_reasons']) for k,v in models.items()},flush=True)
    summary={}
    for variant in ['fixed','fitted_vascular_delay']:
        summary[variant]=dict(median_cross_validated_nrmse=float(np.median([r['models'][variant]['cross_validated_nrmse'] for r in records])),
                              unusable_scans=[r['scan'] for r in records if not r['models'][variant]['full_fit']['usable_for_quantification']],
                              converged=sum(r['models'][variant]['full_fit']['success'] for r in records))
    (out/'real_v02_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
