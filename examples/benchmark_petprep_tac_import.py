"""Execute a pinned, unmodified PETPrep extraction interface and verify import.

Usage: python benchmark_petprep_tac_import.py PINNED_TACS_PY OUTPUT_DIR
Requires nipype, pandas, nibabel. This runs the TAC interface on synthetic images,
not full PETPrep preprocessing or human anatomical validation.
"""
import argparse,hashlib,importlib.util,json
from importlib.metadata import version
from pathlib import Path
import numpy as np,nibabel as nib
from petkin.petprep import import_petprep_tacs


def run(source,output):
    source=Path(source).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=True)
    module_spec=importlib.util.spec_from_file_location('upstream_petprep_tacs',source);module=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
    labels=np.ones((6,6,6),dtype=np.int16);labels[3:]=2;time=np.arange(12.);values=[10+4*np.exp(-time/3),3+2*np.exp(-time/4)]
    data=np.stack([np.where(labels==1,values[0][i],values[1][i]) for i in range(12)],axis=-1).astype('float32')
    nib.save(nib.Nifti1Image(data,np.diag([2.,2.,2.,1.])),out/'pet.nii.gz');nib.save(nib.Nifti1Image(labels,np.diag([2.,2.,2.,1.])),out/'labels.nii.gz')
    (out/'labels.tsv').write_text('index\tname\n1\tCortex left\n2\tThalamus\n')
    (out/'pet.json').write_text(json.dumps(dict(TimeZero='12:00:00',InjectionStart=0,FrameTimesStart=(time*60).tolist(),FrameDuration=[60]*12,Units='Bq/mL',ImageDecayCorrected=True,ImageDecayCorrectionTime=0)))
    run_dir=out/'upstream';run_dir.mkdir(exist_ok=True)
    extraction=module.ExtractTACs(in_file=str(out/'pet.nii.gz'),segmentation=str(out/'labels.nii.gz'),dseg_tsv=str(out/'labels.tsv'),metadata=str(out/'pet.json')).run(cwd=str(run_dir))
    result=import_petprep_tacs(extraction.outputs.out_file,out/'pet.json',out/'imported')
    comparisons=[]
    for index,region in enumerate(result['regions']):
        imported=np.genfromtxt(out/'imported'/region['path'],names=True,delimiter=',');np.testing.assert_allclose(imported['start'],time,rtol=0,atol=1e-8);np.testing.assert_allclose(imported['activity'],values[index],rtol=1e-5)
        comparisons.append(dict(region=region['name'],max_absolute_activity_difference=float(np.max(np.abs(imported['activity']-values[index])))))
    summary=dict(interface='Unmodified PETPrep ExtractTACs',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),dependencies={x:version(x) for x in ['nipype','pandas','nibabel']},regions=2,frames_per_region=12,comparisons=comparisons,passed=True,scope='Synthetic interface compatibility only; full PETPrep was not run; no human registration or biological accuracy claim')
    (out/'comparison.json').write_text(json.dumps(summary,indent=2));return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output');a=p.parse_args();print(json.dumps(run(a.source,a.output),indent=2))
