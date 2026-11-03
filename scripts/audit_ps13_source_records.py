"""Public-record checks; no independent anatomical or calibration certification."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np,nibabel as nib
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('data_dir');parser.add_argument('audit_dir');parser.add_argument('masks_dir');args=parser.parse_args()
root=Path(args.data_dir);out=Path(args.audit_dir);out.mkdir(parents=True,exist_ok=True)
results=[]
for path in sorted(root.glob('sub-*/ses-*/pet/*_pet.json')):
 meta=json.loads(path.read_text());stem=path.name[:-5];scan=stem.removesuffix('_pet');nii_path=path.with_suffix('.nii.gz');blood_path=path.with_name(scan+'_recording-manual_blood.tsv');bm_path=path.with_name(scan+'_recording-manual_blood.json');bm=json.loads(bm_path.read_text());blood=np.genfromtxt(blood_path,names=True,delimiter='\t');start=np.asarray(meta['FrameTimesStart'],float);duration=np.asarray(meta['FrameDuration'],float);time=blood['time'];n=len(time)
 flags=[];checks={}
 checks['PET_units_are_Bq_per_mL']=meta.get('Units') in ('BQML','Bq/mL')
 checks['blood_plasma_units_are_Bq_per_mL']=bm.get('plasma_radioactivity',{}).get('Units')=='Bq/mL'
 checks['blood_time_units_are_seconds']=bm.get('time',{}).get('Units')=='s'
 checks['PET_decay_corrected_to_time_zero']=meta.get('ImageDecayCorrected') is True and meta.get('ImageDecayCorrectionTime')==0
 checks['frame_intervals_contiguous_from_zero']=bool(start[0]==0 and np.allclose(start[1:],(start+duration)[:-1],atol=1e-8) and np.all(duration>0))
 checks['blood_sampling_times_increase']=bool(np.isfinite(time).all() and np.all(np.diff(time)>0))
 checks['blood_concentrations_nonnegative_finite']=bool(all(np.isfinite(blood[k]).all() and np.all(blood[k]>=0) for k in ['plasma_radioactivity','whole_blood_radioactivity']))
 checks['parent_fractions_in_zero_one']=bool(np.isfinite(blood['metabolite_parent_fraction']).all() and np.all((blood['metabolite_parent_fraction']>=0)&(blood['metabolite_parent_fraction']<=1)))
 row=dict(scan=scan,image_present=nii_path.exists(),blood_samples=n,first_blood_seconds=float(time[0]),last_blood_seconds=float(time[-1]),scan_end_seconds=float((start+duration)[-1]),terminal_sample_gap_seconds=float((start+duration)[-1]-time[-1]),checks=checks,flags=flags,source_json_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_blood_sha256=hashlib.sha256(blood_path.read_bytes()).hexdigest())
 if time[-1]<(start+duration)[-1]:flags.append('terminal_input_hold_assumed')
 if bm.get('DispersionCorrected') is False:flags.append('source_blood_not_dispersion_corrected')
 if bm.get('MetaboliteRecoveryCorrectionApplied') is False:flags.append('source_metabolite_recovery_not_corrected')
 # Matching concentration units do not demonstrate cross-calibration or shared decay reference.
 d=np.asarray(meta.get('DecayFactor',[]),float);ref=np.asarray(meta.get('FrameReferenceTime',[]),float);half=meta.get('RadionuclideHalfLife')
 if half and len(d)==len(ref) and len(d):
  offset=np.log(d)*float(half)/np.log(2)-ref
  row['diagnostic_decay_factor_to_frame_reference_offset_seconds']=dict(median=float(np.median(offset)),range=[float(offset.min()),float(offset.max())],caveat='Implementation-dependent legacy DICOM fields; do not shift or rescale without source confirmation.')
 row['unresolved']=['PET_blood_cross_calibration','blood_decay_reference','input_to_tissue_timing_alignment','vascular_fraction_justification','expert_anatomical_ROI_verification']
 if nii_path.exists():
  nii=nib.load(nii_path);checks['image_frames_match_timing']=nii.ndim==4 and nii.shape[3]==len(start);checks['physical_affine_invertible']=bool(np.isfinite(nii.affine).all() and abs(np.linalg.det(nii.affine[:3,:3]))>1e-8);row.update(shape=list(nii.shape),voxel_spacing_mm=list(map(float,nii.header.get_zooms()[:3])),native_orientation=list(nib.aff2axcodes(nii.affine)))
  maskpath=Path(args.masks_dir)/(scan+'_uptake_mask.nii.gz')
  mask=nib.load(maskpath);volume=float(np.count_nonzero(np.asarray(mask.dataobj))*abs(np.linalg.det(mask.affine[:3,:3]))/1000)
  row['exploratory_ROI_volume_mL']=volume
 else:flags.append('metadata_without_PET_image');checks['image_frames_match_timing']=None
 results.append(row)
summary=dict(role='Source-record audit, not independent human validation',metadata_scans=len(results),available_image_scans=sum(x['image_present'] for x in results),missing_image_scans=[x['scan'] for x in results if not x['image_present']],check_failures=[dict(scan=x['scan'],check=k) for x in results for k,v in x['checks'].items() if v is False],independent_calibration_records_available=False,independent_landmark_pairs_available=0,qualified=False)
(out/'records.json').write_text(json.dumps(results,indent=2));(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
