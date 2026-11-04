"""Historical PS13 PET-threshold experiment; not regional brain validation.
At least one tested mask was subsequently found outside the estimated brain.
Use anatomical labels with petkin extract/analyze for regional brain analysis.

Original frozen external PS13 technical evaluation. See validation/external_ps13/PROTOCOL.json.
No individual kinetic truth; this does not establish full research qualification.
"""
import sys,json,hashlib,gzip,shutil
from pathlib import Path
import numpy as np
import nibabel as nib
from scipy.ndimage import label,binary_erosion
from petkin import Frames,InputFunction,fit_compartment,predict
from petkin.graphical import logan,ma1
from petkin.validation import analyze_research_tac,REQUIRED_ASSUMPTIONS

def run(data,out):
 data,out=Path(data),Path(out);out.mkdir(parents=True,exist_ok=True);results=[]
 paths=sorted(data.glob('sub-*/ses-*/pet/*_pet.nii.gz'))
 if len(paths)!=20:raise ValueError(f'Expected all 20 scans, got {len(paths)}')
 for path in paths:
  scan=path.name.removesuffix('_pet.nii.gz');print('Processing',scan,flush=True)
  record=dict(scan=scan,subject=scan.split('_')[0],session=scan.split('_')[1],released_for_quantification=False)
  try:
   meta=json.loads(path.with_name(scan+'_pet.json').read_text());assert meta['Units'] in ('BQML','Bq/mL');assert meta['ImageDecayCorrected'] and meta['ImageDecayCorrectionTime']==0
   starts=np.asarray(meta['FrameTimesStart'])/60;end=starts+np.asarray(meta['FrameDuration'])/60;frames=Frames(starts,end)
   held=np.arange(len(starts))%5==4;train=~held
   # Decompress to scratch to avoid repeated gzip decompression while reading frames.
   raw=out/(scan+'.nii')
   with gzip.open(path,'rb') as src,raw.open('wb') as dst:shutil.copyfileobj(src,dst)
   nii=nib.load(raw);mean=np.zeros(nii.shape[:3],dtype=np.float64)
   for i in np.flatnonzero(train):mean+=np.asarray(nii.dataobj[...,i],dtype=float)*frames.duration[i]
   mean/=frames.duration[train].sum();regions,count=label(mean>.35*mean.max());sizes=np.bincount(regions.ravel());sizes[0]=0
   mask=binary_erosion(regions==sizes.argmax());assert mask.sum()>100
   tac=np.array([np.asarray(nii.dataobj[...,i],float)[mask].mean() for i in range(len(starts))]);nib.save(nib.Nifti1Image(mask.astype('uint8'),nii.affine),str(out/(scan+'_uptake_mask.nii.gz')))
   del nii;raw.unlink()
   blood=np.genfromtxt(path.with_name(scan+'_recording-manual_blood.tsv'),names=True,delimiter='\t')
   t=blood['time']/60;cp=blood['plasma_radioactivity']*blood['metabolite_parent_fraction'];wb=blood['whole_blood_radioactivity']
   assert np.all((blood['metabolite_parent_fraction']>=0)&(blood['metabolite_parent_fraction']<=1));assert t[0]>0 and np.all(np.diff(t)>0)
   terminal=max(end[-1]+.5,t[-1]+1e-6);inp=InputFunction(np.r_[0,t,terminal],np.r_[0,cp,cp[-1]],np.r_[0,wb,wb[-1]])
   sigma=np.sqrt(np.maximum(tac,1)/frames.duration);opts=dict(model='2t',vb=.05,dt=.025,n_starts=8,seed=0,input_tail='hold')
   cv={}
   for mode in ('fixed_delay','fitted_delay'):
    fit=fit_compartment(Frames(starts[train],end[train]),inp,tac[train],sigma=sigma[train],fit_delay=mode=='fitted_delay',**opts)
    pred=predict(frames,inp,list(fit['parameters'].values()),model='2t',vb=.05,delay=fit.get('delay',0),input_tail='hold')
    cv[mode]=dict(fit=fit,held_frame_nrmse=float(np.sqrt(np.mean((pred[held]-tac[held])**2))/np.mean(tac[held])),prediction=pred.tolist())
   assumptions={k:False for k in REQUIRED_ASSUMPTIONS};assumptions.update(frame_timing_verified=True,parent_correction_appropriate=True)
   full=analyze_research_tac(frames,inp,tac,sigma=sigma,assumptions=assumptions,**opts)
   graphics={name:[fn(frames,tac,inp,tstar,vb=.05) for tstar in (20,30,40)] for name,fn in [('logan',logan),('ma1',ma1)]}
   np.savetxt(out/(scan+'_tac.csv'),np.c_[starts,end,tac,sigma,held.astype(int)],delimiter=',',header='start,end,activity,sigma,held',comments='')
   grid=np.linspace(0,terminal,int(np.ceil(terminal/.001))+1)
   np.savetxt(out/(scan+'_input.csv'),np.c_[grid,np.interp(grid,inp.time,inp.plasma),np.interp(grid,inp.time,inp.blood)],delimiter=',',header='time,plasma,blood',comments='')
   record.update(status='completed',mask_voxels=int(mask.sum()),held_indices=np.flatnonzero(held).tolist(),cross_validation=cv,full_research_fit=full,graphical=graphics,terminal_hold_seconds=max(0,(end[-1]-t[-1])*60))
   print(scan,'NRMSE',cv['fixed_delay']['held_frame_nrmse'],'screen',full['research_screen']['passed'],flush=True)
  except Exception as exc:record.update(status='failed',error=str(exc));print('FAILED',scan,str(exc),flush=True)
  results.append(record);(out/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False))
 return results
if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('data_dir');parser.add_argument('output_dir')
 parser.add_argument('--reproduce-historical-experiment',action='store_true')
 args=parser.parse_args()
 if not args.reproduce_historical_experiment:
  parser.error('This historical PET-threshold experiment is not brain-PET validation. Use anatomical labels with petkin extract/analyze. Pass --reproduce-historical-experiment only to reproduce the historical records.')
 run(args.data_dir,args.output_dir)

