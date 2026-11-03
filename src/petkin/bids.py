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
    if metadata.get('Units') not in ('Bq/mL','kBq/mL','kBq/ml','Bq/ml'):
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
    return dict(frames=frames,units=metadata['Units'],decay_reference_seconds=float(origin),
                time_origin=metadata.get('TimeZero'),tracer=metadata.get('TracerName'))


def input_from_bids(tsv_path,metadata,*,plasma_column='plasma_radioactivity',
                    blood_column='whole_blood_radioactivity',target_units='Bq/mL',
                    initial_policy='error',terminal_time=None,terminal_policy='error'):
    """Read already processed input; never invent metabolite/delay correction.

    Missing initial samples require explicit zero/hold policy. Terminal extension
    also requires an explicit hold policy. All assumptions are returned for audit.
    Caller must establish parent correction and matching decay-reference times.
    """
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
                                               parent_correction='Caller must verify',decay_reference='Caller must verify')
