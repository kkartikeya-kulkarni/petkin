"""Measured PET/MRI review pack; registration accuracy is NOT established without independent truth."""
import sys,json,gzip,shutil
from pathlib import Path
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from petkin.imaging import register_rigid
from validate_public import to_sitk,reduce

def run(data,out):
 data,out=Path(data),Path(out);out.mkdir(parents=True,exist_ok=True);results=[]
 for petpath in sorted(data.glob('sub-*/ses-*/pet/*_pet.nii.gz')):
  scan=petpath.name.removesuffix('_pet.nii.gz');subject=scan.split('_')[0];dest=out/scan;dest.mkdir(exist_ok=True)
  print('Registering',scan,flush=True);row=dict(scan=scan,accuracy_status='not_evaluated',reference_landmarks=None)
  try:
   mripath=next((data/subject/'ses-baseline'/'anat').glob('*_T1w.nii*'))
   raw=dest/'scratch.nii'
   with gzip.open(petpath,'rb') as src,raw.open('wb') as dst:shutil.copyfileobj(src,dst)
   nii=nib.load(raw);meta=json.loads(petpath.with_name(scan+'_pet.json').read_text());dur=np.asarray(meta['FrameDuration']);start=np.asarray(meta['FrameTimesStart']);eligible=start>=1200
   mean=np.zeros(nii.shape[:3],float)
   for i in np.flatnonzero(eligible):mean+=np.asarray(nii.dataobj[...,i],float)*dur[i]
   mean/=dur[eligible].sum();pet=reduce(to_sitk(mean,nii.affine));del nii;raw.unlink()
   m=nib.load(mripath);mri=reduce(to_sitk(m.get_fdata(dtype=np.float32),m.affine));reg=register_rigid(mri,pet,multimodal=True)
   aligned=sitk.Resample(pet,mri,reg.transform,sitk.sitkLinear,0,sitk.sitkFloat32)
   sitk.WriteImage(pet,str(dest/'pet_mean_native.nii.gz'))
   sitk.WriteTransform(reg.transform,str(dest/'mri_to_pet.tfm'));sitk.WriteImage(aligned,str(dest/'pet_mean_in_mri.nii.gz'));sitk.WriteImage(mri,str(dest/'mri_4mm.nii.gz'))
   a=sitk.GetArrayFromImage(mri);b=sitk.GetArrayFromImage(aligned);fig,axes=plt.subplots(1,3,figsize=(12,4))
   for ax,k in zip(axes,[.35,.5,.65]):
    z=int(a.shape[0]*k);ax.imshow(a[z],cmap='gray',vmin=0,vmax=np.percentile(a,99.5));im=b[z];threshold=np.max(b)*.15
    ax.imshow(np.ma.masked_where(im<threshold,im),cmap='inferno',alpha=.5,vmin=0,vmax=np.max(b));ax.set_title(f'Native plane 3, slice {z}');ax.axis('off')
   fig.suptitle(scan+' — automatic alignment; independent review pending');fig.tight_layout();fig.savefig(dest/'overlay.png',dpi=120);plt.close(fig)
   row.update(status='completed',metric=reg.metric,stop=reg.stop,transform='mri_to_pet.tfm',convention='fixed MRI physical LPS mm -> moving PET physical LPS mm')
  except Exception as exc:row.update(status='failed',error=str(exc));print('FAILED',scan,str(exc),flush=True)
  results.append(row);(out/'registrations.json').write_text(json.dumps(results,indent=2))
 return results
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
