"""Create MRI/PET orthogonal overlays for human review, not an accuracy score."""
from pathlib import Path
import sys
import numpy as np
import SimpleITK as sitk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run(folder):
    folder=Path(folder)
    m=sitk.GetArrayFromImage(sitk.ReadImage(str(folder/'mri_4mm.nii.gz')))
    p=sitk.GetArrayFromImage(sitk.ReadImage(str(folder/'pet_in_mri_4mm.nii.gz'))).mean(axis=0)
    fig,axes=plt.subplots(1,3,figsize=(10,4))
    for axis,ax in enumerate(axes):
        m_slice=np.take(m,m.shape[axis]//2,axis=axis)
        p_slice=np.take(p,p.shape[axis]//2,axis=axis)
        ax.imshow(m_slice,cmap='gray',vmin=0,vmax=np.percentile(m[m>0],98),origin='lower')
        ax.imshow(np.ma.masked_less(p_slice,.15*np.max(p)),cmap='inferno',alpha=.55,vmin=0,vmax=np.percentile(p[p>0],98),origin='lower')
        ax.axis('off'); ax.set_title(['Axial','Coronal','Sagittal'][axis])
    fig.suptitle('Rigid PET–MRI overlay (voxel views; orientation labels omitted)')
    fig.tight_layout(); fig.savefig(folder/'registration_overlay.png',dpi=150); plt.close(fig)
if __name__=='__main__': run(sys.argv[1])
