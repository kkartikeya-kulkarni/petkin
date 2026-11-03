"""Import PETPrep regional TAC derivatives without re-running image processing."""
import csv,json
from pathlib import Path
import numpy as np
from .bids import validate_pet_metadata
from .reporting import file_manifest


def import_petprep_tacs(tsv_path,pet_json,output_dir):
    """Convert documented frame_start/frame_end TSVs to regional model inputs.

    PETPrep's frame columns are in BIDS seconds. Exact agreement with the source
    PET JSON is required, preventing accidental seconds/minutes or series mixes.
    Region columns keep their names; n/a/nonfinite observations remain explicit
    invalid regions for analyze to report, never silently replaced by zero.
    """
    metadata=json.loads(Path(pet_json).read_text())
    if metadata.get('Units')=='BQML':metadata=dict(metadata,Units='Bq/mL')
    checked=validate_pet_metadata(metadata);frames=checked['frames']
    with Path(tsv_path).open(newline='') as file:
        reader=csv.reader(file,delimiter='\t');header=next(reader,None);rows=list(reader)
    if not header or len(header)!=len(set(header)) or not {'frame_start','frame_end'}<=set(header):raise ValueError('PETPrep TSV needs distinct headers including frame_start,frame_end')
    if len(rows)!=len(frames.start) or any(len(row)!=len(header) for row in rows):raise ValueError('PETPrep TAC rows must match PET timing and column count')
    start=np.array([float(row[header.index('frame_start')]) for row in rows])/60
    end=np.array([float(row[header.index('frame_end')]) for row in rows])/60
    if not np.allclose(start,frames.start,rtol=0,atol=1e-8) or not np.allclose(end,frames.end,rtol=0,atol=1e-8):raise ValueError('PETPrep frame timing does not match this PET JSON (expected seconds relative to TimeZero)')
    names=[name for name in header if name not in ('frame_start','frame_end')]
    if not names:raise ValueError('No regional curves found')
    curves=[]
    for name in names:
        column=header.index(name);values=[]
        for row in rows:
            text=row[column].strip()
            try:values.append(float(text) if text.lower() not in ('','n/a','na','nan') else np.nan)
            except ValueError as exc:raise ValueError('Non-numeric regional activity in '+name) from exc
        curves.append(np.array(values))
    out=Path(output_dir)
    if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be empty')
    out.mkdir(parents=True,exist_ok=True);regions=[]
    for index,(name,values) in enumerate(zip(names,curves)):
        filename=f'region_{index:04d}.csv';np.savetxt(out/filename,np.c_[start,end,values],delimiter=',',header='start,end,activity',comments='')
        regions.append(dict(id=name,name=name,path=filename,nonfinite_frames=np.flatnonzero(~np.isfinite(values)).tolist()))
    result=dict(schema_version=1,source='PETPrep regional TAC TSV',units=checked['units'],time_origin=checked['time_origin'],regions=regions,
        provenance=file_manifest([tsv_path,pet_json]),scope='Imports supplied derivatives; does not certify segmentation or registration accuracy')
    result['pet_metadata_sha256']=result['provenance']['inputs'][-1]['sha256']
    result['decay_reference_seconds']=checked['decay_reference_seconds']
    (out/'regions.json').write_text(json.dumps(result,indent=2,allow_nan=False));return result
