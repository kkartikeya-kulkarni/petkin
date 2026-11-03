"""Public FDG out-of-sample prediction check. No ground-truth parameter claim.

Download first using fetch_public.py. Run: python validate_public.py DATA_DIR OUT_DIR
Images are reduced to ~4-mm spacing for this demonstration. ROI is a thresholded
PET uptake region, not a validated anatomical segmentation.
"""
import gzip,json,shutil,sys
from pathlib import Path
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from petkin import Frames,InputFunction,fit_compartment,predict
from petkin.imaging import motion_correct,pet_to_mri,split_dynamic


def to_sitk(data,affine):
    im=sitk.GetImageFromArray(np.transpose(np.asarray(data,dtype=np.float32),(2,1,0)).copy())
    linear=np.diag([-1.,-1.,1.])@affine[:3,:3]
    spacing=np.linalg.norm(linear,axis=0); direction=linear/spacing
    if not np.allclose(direction.T@direction,np.eye(3),atol=1e-4): raise ValueError('Sheared NIfTI unsupported')
    im.SetSpacing(spacing.tolist()); im.SetDirection(direction.ravel().tolist())
    im.SetOrigin((np.diag([-1.,-1.,1.])@affine[:3,3]).tolist())
    return im


def reduce(im):
    spacing=np.array(im.GetSpacing()); target=np.maximum(spacing,4.)
    size=np.ceil((np.array(im.GetSize())-1)*spacing/target).astype(int)+1
    return sitk.Resample(im,size.tolist(),sitk.Transform(3,sitk.sitkIdentity),sitk.sitkLinear,
                         im.GetOrigin(),target.tolist(),im.GetDirection(),0.,sitk.sitkFloat32)


def run(data_dir,out_dir,subjects=('s028','s029')):
    data_dir,out_dir=Path(data_dir),Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    results=[]; fig,axes=plt.subplots(1,len(subjects),figsize=(5.5*len(subjects),4),squeeze=False); axes=axes[0]
    for j,subject in enumerate(subjects):
        prefix=f'sub-{subject}_ses-open'; print('Processing',subject,flush=True)
        metadata=json.loads((data_dir/f'{prefix}_task-rest_pet.json').read_text())
        assert metadata['Units']=='Bq/mL' and metadata['ImageDecayCorrected'] and metadata['ImageDecayCorrectionTime']==0
        start=np.array(metadata['FrameTimesStart'])/60; duration=np.array(metadata['FrameDuration'])/60
        frames=Frames(start,start+duration)
        gz=data_dir/f'{prefix}_task-rest_pet.nii.gz'; raw=data_dir/f'{prefix}_task-rest_pet.nii'
        if not raw.exists():
            with gzip.open(gz,'rb') as src,raw.open('wb') as dst: shutil.copyfileobj(src,dst)
        nii=nib.load(raw)
        small=[reduce(to_sitk(nii.dataobj[:,:,:,i],nii.affine)) for i in range(nii.shape[3])]
        pet=sitk.JoinSeries(small)
        m=nib.load(data_dir/f'{prefix}_T1w.nii.gz'); mri=reduce(to_sitk(m.get_fdata(dtype=np.float32),m.affine))
        corrected,transforms,qc=motion_correct(pet,reference_index=len(small)-1,blank_policy="identity")
        registered,reg=pet_to_mri(corrected,mri,duration)
        stack=np.stack([sitk.GetArrayFromImage(f) for f in split_dynamic(corrected)])
        mean=np.average(stack,axis=0,weights=duration)
        mask=mean>.35*mean.max()
        tac=stack[:,mask].mean(axis=1)
        blood=np.genfromtxt(data_dir/f'{prefix}_task-rest_recording-autosampler_blood.tsv',names=True,delimiter='\t')
        # Dataset supplies midpoint values, not a high-frequency measured bolus.
        # Explicit boundary assumptions: zero plasma at injection; terminal hold.
        time=np.r_[0,blood['time']/60,frames.end[-1]]
        plasma=np.r_[0,blood['plasma_radioactivity'],blood['plasma_radioactivity'][-1]]
        inp=InputFunction(time,plasma)
        def correlation(a,b):
            a=np.asarray(a,dtype=float).ravel();b=np.asarray(b,dtype=float).ravel()
            if np.std(a)<1e-8 or np.std(b)<1e-8: return None
            return float(np.corrcoef(a,b)[0,1])
        ref=sitk.GetArrayFromImage(small[-1])
        for i,q in enumerate(qc):
            q['correlation_before']=correlation(sitk.GetArrayFromImage(small[i]),ref)
            q['correlation_after']=correlation(stack[i],ref)
        eligible=np.std(stack,axis=(1,2,3))>1e-8
        held=(np.arange(len(tac))%5==4)&eligible; train=(~held)&eligible
        fit=fit_compartment(Frames(start[train],frames.end[train]),inp,tac[train],model='2t_irreversible')
        estimate=predict(frames,inp,list(fit['parameters'].values()),model='2t_irreversible')
        nrmse=float(np.sqrt(np.mean((estimate[held]-tac[held])**2))/np.mean(tac[held]))
        record=dict(subject=subject,role='Held out from simulation and implementation development; interleaved frames held out of fitting',
                    model='2t_irreversible',vb=0,fit=fit,excluded_blank_frames=np.flatnonzero(~eligible).tolist(),held_out_frame_indices=np.flatnonzero(held).tolist(),
                    held_out_normalized_rmse=nrmse,registration_metric=reg.metric,registration_stop=reg.stop,
                    motion_qc=qc,roi_voxels=int(mask.sum()),
                    limitations=['No true kinetic parameters','Frame-midpoint plasma interpolation and assumed boundaries',
                                 'No whole blood: fixed vb=0','Plasma decay convention not independently verified','PET uptake threshold ROI, not anatomical brain segmentation',
                                 '4-mm resampling; rigid registration requires visual review','FDG is not a reversible Logan validation'])
        results.append(record)
        target=out_dir/subject; target.mkdir(exist_ok=True)
        for i,tx in enumerate(transforms): sitk.WriteTransform(tx,str(target/f'reference_to_frame_{i:02d}.tfm'))
        sitk.WriteTransform(reg.transform,str(target/'mri_to_pet_reference.tfm'))
        sitk.WriteImage(registered,str(target/'pet_in_mri_4mm.nii.gz'))
        sitk.WriteImage(mri,str(target/'mri_4mm.nii.gz'))
        np.savetxt(target/'tac.csv',np.c_[start,frames.end,tac,estimate,held.astype(int)],delimiter=',',header='start,end,activity,predicted,held_out',comments='')
        np.savetxt(target/'input.csv',np.c_[time,plasma],delimiter=',',header='time,plasma',comments='')
        axes[j].plot(frames.mid,tac,'.',label='Observed'); axes[j].plot(frames.mid,estimate,label='Train-frame fit')
        axes[j].scatter(frames.mid[held],tac[held],facecolors='none',edgecolors='red',s=60,label='Held out')
        axes[j].set(title=f'{subject}: held-out NRMSE {nrmse:.1%}',xlabel='Time (min)',ylabel='Activity (Bq/mL)'); axes[j].legend()
        (out_dir/'public_results.json').write_text(json.dumps(results,indent=2))
        raw.unlink(); print(subject,'NRMSE',nrmse,flush=True)
    fig.tight_layout(); fig.savefig(out_dir/'public_validation.png',dpi=150)

if __name__=='__main__': run(sys.argv[1],sys.argv[2])
