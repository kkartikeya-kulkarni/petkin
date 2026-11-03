"""Conditional kinetic cross-validation and assumption sensitivity on human FDG.

Uses frozen image preprocessing; frame splits do not establish independence of
image preprocessing. PET uptake ROIs are defined using training frames per fold.
No optimization settings are selected using the resulting scores.
"""
import json,sys
from pathlib import Path
import numpy as np
import SimpleITK as sitk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from petkin import Frames,InputFunction,fit_compartment,predict


def fit_predict(frames,inp,y,train,model='2t_irreversible',dt=.025):
    fit=fit_compartment(Frames(frames.start[train],frames.end[train]),inp,y[train],model=model,dt=dt)
    prediction=predict(frames,inp,list(fit['parameters'].values()),model=model,dt=dt)
    return fit,prediction


def nrmse(y,p,mask):
    return float(np.sqrt(np.mean((y[mask]-p[mask])**2))/np.mean(y[mask]))


def run(root,out,subjects=('s028','s029','s030')):
    root,out=Path(root),Path(out);out.mkdir(parents=True,exist_ok=True)
    records=[];fig,axes=plt.subplots(2,2,figsize=(11,8)); axes=axes.ravel()
    for j,subject in enumerate(subjects):
        candidates=list(root.glob(f'**/{subject}/pet_in_mri_4mm.nii.gz'))
        if len(candidates)!=1: raise ValueError(f'Need exactly one processed subject: {subject}')
        folder=candidates[0].parent
        data=np.genfromtxt(folder/'tac.csv',delimiter=',',names=True)
        frames=Frames(data['start'],data['end'])
        blood=np.genfromtxt(folder/'input.csv',delimiter=',',names=True)
        inp=InputFunction(blood['time'],blood['plasma'])
        stack=sitk.GetArrayFromImage(sitk.ReadImage(str(candidates[0])))
        eligible=np.std(stack,axis=(1,2,3))>1e-8
        indices=np.arange(len(eligible)); folds=[]; cv=np.full(len(eligible),np.nan); cv_y=cv.copy()
        for fold in range(5):
            held=(indices%5==fold)&eligible;train=(~held)&eligible
            # Define a different uptake ROI for each fold using training images only.
            mean=np.average(stack[train],axis=0,weights=frames.duration[train]);mask=mean>.35*mean.max()
            y=stack[:,mask].mean(axis=1)
            fit,p=fit_predict(frames,inp,y,train)
            baseline,bp=fit_predict(frames,inp,y,train,model='1t')
            folds.append(dict(fold=fold,held_out_indices=indices[held].tolist(),train_indices=indices[train].tolist(),
                              roi_voxels=int(mask.sum()),nrmse=nrmse(y,p,held),one_tissue_nrmse=nrmse(y,bp,held),
                              Ki=fit['Ki'],parameters=fit['parameters'],success=fit['success'],at_bound=fit['at_bound'],
                              condition_number=fit['condition_number']))
            cv[held]=p[held];cv_y[held]=y[held]
        train=(indices<27)&eligible;held=(indices>=27)&eligible
        mean=np.average(stack[train],axis=0,weights=frames.duration[train]);mask=mean>.35*mean.max()
        y=stack[:,mask].mean(axis=1)
        late_fit,late_p=fit_predict(frames,inp,y,train)
        sensitivity=[]
        mean=np.average(stack[eligible],axis=0,weights=frames.duration[eligible])
        for threshold in [.25,.35,.45]:
            y=stack[:,mean>threshold*mean.max()].mean(axis=1)
            fit,p=fit_predict(frames,inp,y,eligible)
            sensitivity.append(dict(roi_threshold=threshold,Ki=fit['Ki'],parameters=fit['parameters'],at_bound=fit['at_bound']))
        y=stack[:,mean>.35*mean.max()].mean(axis=1)
        reference,_=fit_predict(frames,inp,y,eligible)
        refined,_=fit_predict(frames,inp,y,eligible,dt=.0125)
        changed_plasma=inp.plasma.copy();changed_plasma[0]=changed_plasma[1]
        alternative,_=fit_predict(frames,InputFunction(inp.time,changed_plasma),y,eligible)
        kis=np.array([r['Ki'] for r in folds])
        record=dict(subject=subject,prior_exposure='new subject' if j>1 else 'previously evaluated subject',
                    excluded_blank_frames=indices[~eligible].tolist(),folds=folds,
                    pooled_cross_validated_nrmse=nrmse(cv_y,cv,eligible),
                    median_fold_nrmse=float(np.median([r['nrmse'] for r in folds])),
                    median_one_tissue_fold_nrmse=float(np.median([r['one_tissue_nrmse'] for r in folds])),
                    late_block_nrmse=nrmse(y=stack[:,mask].mean(axis=1),p=late_p,mask=held),
                    late_block_Ki=late_fit['Ki'],late_block_indices=indices[held].tolist(),
                    late_train_end_min=float(frames.end[train][-1]),late_test_start_min=float(frames.start[held][0]),late_test_end_min=float(frames.end[held][-1]),
                    fold_Ki_range=kis.tolist(),fold_Ki_span_relative_to_median=float(np.ptp(kis)/np.median(kis)),
                    roi_sensitivity=sensitivity,all_frame_Ki=reference['Ki'],
                    step_halving_Ki_relative_change=refined['Ki']/reference['Ki']-1,
                    initial_plasma_hold_Ki_relative_change=alternative['Ki']/reference['Ki']-1,
                    all_frame_parameters=reference['parameters'])
        records.append(record)
        ax=axes[j];ax.plot(frames.mid[eligible],cv_y[eligible],'.-',label='Observed (training-defined ROI)')
        ax.plot(frames.mid[eligible],cv[eligible],'x',label='Out-of-fold prediction')
        ax.set(title=f"{subject}: 5-fold NRMSE {record['pooled_cross_validated_nrmse']:.1%}",xlabel='Time (min)',ylabel='Activity (Bq/mL)');ax.legend(fontsize=8)
        print(subject,record['pooled_cross_validated_nrmse'],record['late_block_nrmse'],flush=True)
        (out/'stress_results.json').write_text(json.dumps(records,indent=2))
    for ax in axes[len(records):]: ax.axis('off')
    fig.tight_layout();fig.savefig(out/'cross_validation.png',dpi=150)
if __name__=='__main__': run(sys.argv[1],sys.argv[2])
