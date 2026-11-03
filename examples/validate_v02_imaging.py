"""Exercise composed-transform pipeline on three real FDG PET/MRI scans at 4 mm.

This is functional processing and visual QC, not human anatomical ground truth.
Benchmark inputs are first downsampled; one final interpolation follows from
those inputs. The library can operate directly on original-resolution images.
"""
import gzip,json,shutil,sys
from pathlib import Path
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validate_public import to_sitk,reduce
from petkin.imaging import preprocess_dynamic,split_dynamic


def run(data,out):
    data,out=Path(data),Path(out);rows=[]
    for name in ['s028','s029','s030']:
        prefix=f'sub-{name}_ses-open';target=out/'imaging'/name;target.mkdir(parents=True,exist_ok=True)
        raw=data/f'{prefix}_v02_temp.nii'
        try:
            with gzip.open(data/f'{prefix}_task-rest_pet.nii.gz','rb') as src,raw.open('wb') as dst:shutil.copyfileobj(src,dst)
            nii=nib.load(raw);frames=[reduce(to_sitk(nii.dataobj[:,:,:,i],nii.affine)) for i in range(nii.shape[3])]
            dynamic=sitk.JoinSeries(frames)
            mr=nib.load(data/f'{prefix}_T1w.nii.gz');mri=reduce(to_sitk(mr.get_fdata(dtype=np.float32),mr.affine))
            metadata=json.loads((data/f'{prefix}_task-rest_pet.json').read_text());duration=np.array(metadata['FrameDuration'])
            result=preprocess_dynamic(dynamic,mri,duration,reference_index=len(frames)-1,blank_policy='identity')
            stack=sitk.GetArrayFromImage(result['image']);mean=np.average(stack,axis=0,weights=duration)
            image=sitk.GetImageFromArray(mean.astype('float32'));image.CopyInformation(mri)
            sitk.WriteImage(image,str(target/'mean_pet_in_mri.nii.gz'));sitk.WriteImage(mri,str(target/'mri_4mm.nii.gz'))
            for i,tx in enumerate(result['composed_transforms']):tx.FlattenTransform();sitk.WriteTransform(tx,str(target/f'mri_to_frame_{i:02d}.h5'))
            qc={k:result[k] for k in ['motion_qc','registration_metric','registration_stop','registration_candidates','final_interpolations_per_frame','limitations']}
            qc.update(subject=name,input_resolution='Downsampled to >=4 mm before pipeline',frames=len(frames),human_ground_truth=False)
            (target/'qc.json').write_text(json.dumps(qc,indent=2));rows.append(qc)
            anatomy=sitk.GetArrayFromImage(mri);fig,axes=plt.subplots(1,3,figsize=(10,4))
            for axis,ax in enumerate(axes):
                a=np.take(anatomy,anatomy.shape[axis]//2,axis=axis);p=np.take(mean,mean.shape[axis]//2,axis=axis)
                ax.imshow(a,cmap='gray',vmin=0,vmax=np.percentile(anatomy[anatomy>0],98),origin='lower')
                ax.imshow(np.ma.masked_less(p,.15*mean.max()),cmap='inferno',alpha=.55,vmin=0,vmax=np.percentile(mean[mean>0],98),origin='lower');ax.axis('off')
            fig.suptitle(f'{name}: composed-transform PET/MRI overlay (voxel views)');fig.tight_layout();fig.savefig(target/'review.png',dpi=150);plt.close(fig)
            print(name,'processed',len(frames),'frames; candidate metrics',[r.get('metric') for r in result['registration_candidates']],flush=True)
        finally:
            if raw.exists():raw.unlink()
    (out/'imaging_summary.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
