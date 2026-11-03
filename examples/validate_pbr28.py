"""Real reversible-tracer check: 20 PBR28 scans / 10 people, frontal cortex.

No ground truth. Uses source processed metabolite-corrected AIF and whole blood,
source positive weights, fixed vb=.05 and no input-delay estimation. Logan's
missing initial PET interval is explicitly filled with zero for integration.
"""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from petkin import Frames,InputFunction,fit_compartment,predict,logan


def run(data_dir,out_dir):
    data_dir,out_dir=Path(data_dir),Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    scans=list(csv.DictReader((data_dir/'scans.csv').open()))
    results=[];fig,axes=plt.subplots(1,2,figsize=(11,4))
    for item in scans:
        name=str(item['scan']);folder=data_dir/name
        tac=np.genfromtxt(folder/'tacs.csv',delimiter=',',names=True)
        b=np.genfromtxt(folder/'input.csv',delimiter=',',names=True)
        valid=tac['Duration']>0;tac=tac[valid]
        frames=Frames(tac['StartTime']/60,(tac['StartTime']+tac['Duration'])/60)
        time,cp,cb=b['Time'],b['AIF'],b['Blood']
        extension=max(0.,frames.end[-1]-time[-1])
        if extension:
            time=np.r_[time,frames.end[-1]];cp=np.r_[cp,cp[-1]];cb=np.r_[cb,cb[-1]]
        inp=InputFunction(time,cp,cb);y=tac['FC'];eligible=tac['Weights']>0
        sigma=1/np.sqrt(tac['Weights'][eligible])
        full=fit_compartment(Frames(frames.start[eligible],frames.end[eligible]),inp,y[eligible],sigma=sigma,vb=.05)
        full_p=predict(frames,inp,list(full['parameters'].values()),vb=.05)
        folds=[];cv=np.full(len(y),np.nan);indices=np.arange(len(y))
        for k in range(5):
            held=(indices%5==k)&eligible;train=(~held)&eligible
            fit=fit_compartment(Frames(frames.start[train],frames.end[train]),inp,y[train],
                                sigma=1/np.sqrt(tac['Weights'][train]),vb=.05)
            p=predict(frames,inp,list(fit['parameters'].values()),vb=.05);cv[held]=p[held]
            folds.append(dict(fold=k,VT=fit['VT'],success=fit['success'],at_bound=fit['at_bound'],
                              held_out_indices=indices[held].tolist(),
                              normalized_rmse=float(np.sqrt(np.mean((p[held]-y[held])**2))/np.mean(y[held]))))
        gap=float(frames.start[0])
        if gap>1e-10:
            graphical_frames=Frames(np.r_[0,frames.start],np.r_[gap,frames.end]);graphical_y=np.r_[0,y]
        else:graphical_frames,graphical_y=frames,y
        graphical=[logan(graphical_frames,graphical_y,inp,tstar=t,vb=.05) for t in [30.,40.,50.]]
        g=graphical[1]
        record=dict(scan=name,subject=str(item['subject']),visit=int(item['visit']),region='FC',
                    full_fit=full,folds=folds,zero_weight_indices=indices[~eligible].tolist(),
                    cross_validated_nrmse=float(np.sqrt(np.mean((cv[eligible]-y[eligible])**2))/np.mean(y[eligible])),
                    logan=graphical,logan_40_vs_compartment_VT_relative_difference=g['raw_VT']/full.get('raw_VT',full['VT'])-1,
                    initial_unmeasured_interval_minutes=gap,terminal_input_hold_minutes=extension)
        results.append(record)
        print(name,record['cross_validated_nrmse'],full['VT'],g['VT'],full['at_bound'],flush=True)
        if not full['at_bound']: axes[0].scatter(full.get('raw_VT',full['VT']),g['raw_VT'])
        axes[1].plot([30,40,50],[v['VT'] for v in graphical],alpha=.6)
        (out_dir/'pbr28_results.json').write_text(json.dumps(results,indent=2))
    interior=[r for r in results if not r['full_fit']['at_bound']]
    limits=[min(r['full_fit'].get('raw_VT',r['full_fit']['VT']) for r in interior),max(r['full_fit'].get('raw_VT',r['full_fit']['VT']) for r in interior)]
    flagged=[r['scan'] for r in results if r['full_fit']['at_bound']]
    axes[0].text(.02,.98,'Boundary fits excluded from plot: '+', '.join(flagged),transform=axes[0].transAxes,va='top',fontsize=8,color='red')
    axes[0].plot(limits,limits,'k--',label='Agreement line');axes[0].legend()
    axes[0].set(xlabel='2-tissue VT (mL/mL)',ylabel='Logan VT, t*=40 min (mL/mL)',title='Real PBR28: method agreement, not ground truth')
    axes[1].set(xlabel='Logan t* (min)',ylabel='VT (mL/mL)',title='Sensitivity to graphical fit start')
    fig.tight_layout();fig.savefig(out_dir/'pbr28_validation.png',dpi=150)
    paired=[]
    for subject in sorted(set(r['subject'] for r in results)):
        rows=[r for r in results if r['subject']==subject];a,c=[r['full_fit'].get('raw_VT',r['full_fit']['VT']) for r in rows]
        paired.append(dict(subject=subject,VTs=[a,c],absolute_test_retest_difference=abs(a-c)/((a+c)/2)))
    summary=dict(scans=len(results),people=len(paired),region='Frontal cortex',fixed_vb=.05,
                 median_cross_validated_nrmse=float(np.median([r['cross_validated_nrmse'] for r in results])),
                 max_cross_validated_nrmse=float(max(r['cross_validated_nrmse'] for r in results)),
                 converged_full_fits=sum(r['full_fit']['success'] for r in results),
                 full_fits_at_bounds=[r['scan'] for r in results if r['full_fit']['at_bound']],
                 fold_fits_at_bounds=[dict(scan=r['scan'],fold=f['fold']) for r in results for f in r['folds'] if f['at_bound']],
                 median_absolute_logan_vs_compartment_difference=float(np.median([abs(r['logan_40_vs_compartment_VT_relative_difference']) for r in results])),
                 median_absolute_test_retest_difference=float(np.median([r['absolute_test_retest_difference'] for r in paired])),pairs=paired)
    (out_dir/'pbr28_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
