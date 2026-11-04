"""Strict PET-BIDS timing and concentration checks; no implicit correction."""
import json,csv
from pathlib import Path
import numpy as np
from .kinetics import Frames,InputFunction


def frames_from_bids(metadata):
    """BIDS seconds -> toolkit minutes, preserving the BIDS TimeZero origin."""
    if not isinstance(metadata,dict):metadata=json.loads(Path(metadata).read_text())
    if 'FrameTimesStart' not in metadata or 'FrameDuration' not in metadata:
        raise ValueError('PET-BIDS FrameTimesStart and FrameDuration are required')
    starts=np.asarray(metadata['FrameTimesStart'],float);durations=np.asarray(metadata['FrameDuration'],float)
    if starts.shape!=durations.shape:raise ValueError('BIDS timing arrays must match')
    return Frames(starts/60,(starts+durations)/60)


def validate_pet_metadata(metadata):
    """Fail instead of assuming concentration units or an unknown decay origin."""
    frames=frames_from_bids(metadata)
    source_units=metadata.get('Units')
    aliases={'Bq/mL':'Bq/mL','Bq/ml':'Bq/mL','BQML':'Bq/mL','kBq/mL':'kBq/mL','kBq/ml':'kBq/mL'}
    if source_units not in aliases:
        raise ValueError('Supported PET concentration units are Bq/mL or kBq/mL')
    if metadata.get('ImageDecayCorrected') is not True:
        raise ValueError('Decay-corrected PET is required; perform and document correction upstream')
    origin=metadata.get('ImageDecayCorrectionTime')
    try:
        if origin is None or not np.isscalar(origin):raise ValueError
        origin=float(origin)
        if not np.isfinite(origin):raise ValueError
    except (TypeError,ValueError):
        raise ValueError('Explicit finite ImageDecayCorrectionTime required') from None
    return dict(frames=frames,units=aliases[source_units],source_units=source_units,
                unit_alias_applied=source_units!=aliases[source_units],decay_reference_seconds=float(origin),
                time_origin=metadata.get('TimeZero'),tracer=metadata.get('TracerName'))


def input_from_bids(tsv_path,metadata,*,plasma_column='plasma_radioactivity',
                    blood_column='whole_blood_radioactivity',target_units='Bq/mL',
                    initial_policy='error',terminal_time=None,terminal_policy='error',
                    plasma_mode='processed',parent_fraction_column='metabolite_parent_fraction'):
    """Read input with an explicit processed-plasma or total-plasma convention.

    Missing initial samples require explicit zero/hold policy. Terminal extension
    also requires an explicit hold policy. All assumptions are returned for audit.
    Caller must establish parent correction and matching decay-reference times.
    """
    if plasma_mode not in ('processed','total'):raise ValueError("plasma_mode must be 'processed' or 'total'")
    if initial_policy not in ('error','zero','hold') or terminal_policy not in ('error','hold'):raise ValueError('Invalid boundary policy')
    if terminal_time is not None and (not np.isfinite(terminal_time) or terminal_time<0):raise ValueError('terminal_time must be finite and nonnegative')
    if not isinstance(metadata,dict):metadata=json.loads(Path(metadata).read_text())
    with Path(tsv_path).open(newline='') as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    if not rows:raise ValueError('Empty input TSV')
    if 'time' not in metadata or metadata['time'].get('Units')!='s':raise ValueError('BIDS input time must explicitly use seconds')
    units=metadata.get(plasma_column,{}).get('Units')
    factors={'Bq/mL':1.,'Bq/ml':1.,'kBq/mL':1000.,'kBq/ml':1000.}
    if units not in factors or target_units not in factors:raise ValueError('Explicit supported plasma concentration units required')
    factor=factors[units]/factors[target_units]
    time=np.array([float(r['time']) for r in rows])/60
    plasma=np.array([float(r[plasma_column]) for r in rows])*factor
    parent_correction='Caller declares plasma is already processed'
    if plasma_mode=='total':
        if metadata.get('MetaboliteAvail') is False or parent_fraction_column not in rows[0]:raise ValueError('Total plasma requires measured parent fractions')
        fraction=np.array([float(r[parent_fraction_column]) for r in rows])
        if not np.isfinite(fraction).all() or np.any((fraction<0)|(fraction>1)):raise ValueError('Parent fractions must be finite and between zero and one')
        plasma=plasma*fraction
        parent_correction='Total plasma multiplied by '+parent_fraction_column
    if metadata.get('WholeBloodAvail') is True and blood_column not in rows[0]:raise ValueError('WholeBloodAvail is true but the blood column is missing')
    blood=None
    if blood_column in rows[0]:
        bu=metadata.get(blood_column,{}).get('Units')
        if bu not in factors:raise ValueError('Explicit supported whole-blood units required')
        blood=np.array([float(r[blood_column]) for r in rows])*factors[bu]/factors[target_units]
    assumptions=[]
    if time[0]<0:raise ValueError('Pre-zero samples require explicit upstream clock alignment')
    if time[0]>0:
        if initial_policy not in ('zero','hold'):raise ValueError('Input lacks time zero: choose an explicit initial policy')
        time=np.r_[0,time];plasma=np.r_[0 if initial_policy=='zero' else plasma[0],plasma]
        if blood is not None:blood=np.r_[0 if initial_policy=='zero' else blood[0],blood]
        assumptions.append('initial_'+initial_policy)
    if terminal_time is not None and time[-1]<terminal_time:
        if terminal_policy!='hold':raise ValueError('Input does not cover terminal_time; explicit terminal hold required')
        time=np.r_[time,terminal_time];plasma=np.r_[plasma,plasma[-1]]
        if blood is not None:blood=np.r_[blood,blood[-1]]
        assumptions.append('terminal_hold')
    return InputFunction(time,plasma,blood),dict(assumptions=assumptions,units=target_units,
                                               parent_correction=parent_correction,plasma_mode=plasma_mode,
                                               dispersion_corrected=metadata.get('DispersionCorrected'),
                                               recovery_corrected=metadata.get('MetaboliteRecoveryCorrectionApplied'),
                                               decay_reference='Caller must verify')


def prepare_bids_input(tsv_path,blood_json,pet_json,output_dir,*,plasma_kind,
                       initial_policy='error',terminal_policy='error'):
    """Create model-ready CSV in the PET concentration units, preserving TimeZero.

    No clock offset, PET/blood calibration, dispersion correction, or recovery
    correction is inferred from legacy fields or a plasma-kind declaration.
    """
    from .reporting import file_manifest
    metadata=json.loads(Path(pet_json).read_text());source_units=metadata.get('Units')
    # Explicit legacy DICOM alias found in public source data; concentration is unchanged.
    if source_units=='BQML':metadata=dict(metadata,Units='Bq/mL')
    checked=validate_pet_metadata(metadata)
    if not isinstance(checked['time_origin'],str) or not checked['time_origin']:raise ValueError('PET-BIDS TimeZero must be supplied')
    injection=metadata.get('InjectionStart')
    if injection is None or not np.isscalar(injection) or not np.isfinite(float(injection)):raise ValueError('Explicit finite InjectionStart is required for model input preparation')
    if float(injection)<0:raise ValueError('Injection predates TimeZero: align blood and tissue to an earlier common origin before fitting; zero initial tissue cannot be assumed at scan start')
    if not np.isclose(checked['decay_reference_seconds'],0,rtol=0,atol=1e-9):raise ValueError('PET is not decay-corrected to TimeZero; align decay references upstream before preparing input')
    if plasma_kind not in ('total','parent'):raise ValueError('Choose total or parent plasma explicitly')
    inp,qc=input_from_bids(tsv_path,blood_json,target_units=checked['units'],
        initial_policy=initial_policy,terminal_time=checked['frames'].end[-1],terminal_policy=terminal_policy,
        plasma_mode='total' if plasma_kind=='total' else 'processed')
    qc.update(time_origin=checked['time_origin'],pet_decay_reference_seconds=checked['decay_reference_seconds'],
        PET_blood_cross_calibration_verified=False,blood_decay_reference_verified=False,
        timing_source='FrameTimesStart and FrameDuration relative to TimeZero; no legacy FrameReferenceTime correction',
        pet_source_units=source_units,provenance=file_manifest([tsv_path,blood_json,pet_json]))
    qc['pet_metadata_sha256']=qc['provenance']['inputs'][-1]['sha256']
    out=Path(output_dir)
    if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be empty')
    out.mkdir(parents=True,exist_ok=True)
    values=[inp.time,inp.plasma];names=['time','plasma']
    if inp.blood is not None:values.append(inp.blood);names.append('blood')
    np.savetxt(out/'input.csv',np.column_stack(values),delimiter=',',header=','.join(names),comments='')
    (out/'input.json').write_text(json.dumps(qc,indent=2,allow_nan=False))
    return qc
