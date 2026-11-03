"""Regional image extraction and complete, auditable scan-level analysis."""
import csv,html,json
from pathlib import Path
import numpy as np
from .kinetics import Frames,InputFunction,predict
from .fitting import fit_compartment
from .graphical import logan,ma1,patlak
from .reporting import file_manifest,write_qc_html


def _new_output(path):
    path=Path(path)
    if path.exists() and any(path.iterdir()):raise ValueError(f'Output directory must be empty: {path}')
    path.mkdir(parents=True,exist_ok=True)
    return path


def extract_regions(pet_path,labels_path,timing_path,output_dir):
    """Extract named-by-ID regions from supplied anatomical integer labels.

    Labels must already share the PET spatial grid. There is no intensity-based
    segmentation, registration inference, or silent label interpolation here.
    """
    import SimpleITK as sitk
    from .bids import validate_pet_metadata
    from .imaging import roi_tacs,split_dynamic
    pet=sitk.ReadImage(str(pet_path));labels=sitk.ReadImage(str(labels_path))
    metadata=json.loads(Path(timing_path).read_text());checked=validate_pet_metadata(metadata);frames=checked['frames']
    parts=split_dynamic(pet)
    if len(parts)!=len(frames.start):raise ValueError('PET frame count does not match timing metadata')
    if any(not np.isfinite(sitk.GetArrayViewFromImage(f)).all() for f in parts):raise ValueError('PET contains nonfinite voxels')
    curves=roi_tacs(pet,labels)
    if not curves:raise ValueError('Labels contain no positive region IDs')
    out=_new_output(output_dir);lab=sitk.GetArrayFromImage(labels);volume=np.prod(labels.GetSpacing())/1000
    regions=[]
    for key,tac in curves.items():
        name=f'region_{key}';filename=f'{name}.csv';np.savetxt(out/filename,np.c_[frames.start,frames.end,tac],delimiter=',',header='start,end,activity',comments='')
        regions.append(dict(id=key,name=name,path=filename,voxels=int((lab==key).sum()),volume_mL=float((lab==key).sum()*volume)))
    manifest=dict(schema_version=1,role='Regional observations from supplied labels',units=checked['units'],decay_reference_seconds=checked['decay_reference_seconds'],regions=regions,provenance=file_manifest([pet_path,labels_path,timing_path]))
    manifest['pet_metadata_sha256']=manifest['provenance']['inputs'][-1]['sha256']
    manifest['time_origin']=checked['time_origin']
    (out/'regions.json').write_text(json.dumps(manifest,indent=2,allow_nan=False))
    return manifest


def _read_input(path):
    data=np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True))
    if not {'time','plasma'}<=set(data.dtype.names or []):raise ValueError('Input CSV requires time,plasma and optional blood')
    return InputFunction(data['time'],data['plasma'],data['blood'] if 'blood' in data.dtype.names else None)


def _read_tac(path):
    tac=np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True));names=set(tac.dtype.names or [])
    if not {'start','end','activity'}<=names:raise ValueError('TAC CSV requires start,end,activity')
    frames=Frames(tac['start'],tac['end']);y=np.asarray(tac['activity'],float)
    if not np.isfinite(y).all():raise ValueError('TAC activity must be finite')
    if {'weight','sigma'}<=names:raise ValueError('Supply weight or sigma, not both')
    include=np.ones(len(y),dtype=bool);sigma=None
    if 'weight' in names:
        w=np.asarray(tac['weight'],float)
        if not np.isfinite(w).all() or np.any(w<0):raise ValueError('Weights must be finite and nonnegative')
        include=w>0;sigma=np.ones(len(y));sigma[include]=1/np.sqrt(w[include])
    elif 'sigma' in names:
        sigma=np.asarray(tac['sigma'],float)
        if not np.isfinite(sigma).all() or np.any(sigma<=0):raise ValueError('Sigma must be finite and positive')
    return frames,y,sigma,include


def analyze_regions(manifest_path,input_path,output_dir,*,model='2t',vb=0.,tstar=30.,fit_delay=False,input_tail='error',n_starts=8,seed=0,hold_every=5):
    """Write all region fits, predictions, plots, summary CSV and an HTML index.

    Model selection and regional labels are supplied by the analyst. Point
    estimates, QC and held-frame prediction are distinct, explicitly named fields.
    One failed region never disappears and does not suppress successful regions.
    """
    if isinstance(hold_every,bool) or not isinstance(hold_every,int) or hold_every<2:raise ValueError('hold_every must be an integer >=2')
    if model not in ('1t','2t','2t_irreversible'):raise ValueError('Unknown model')
    if not np.isfinite(tstar) or tstar<0:raise ValueError('tstar must be finite and nonnegative')
    inp=_read_input(input_path);manifest_path=Path(manifest_path);manifest=json.loads(manifest_path.read_text());regions=manifest.get('regions',[])
    if not regions or len({str(x['id']) for x in regions})!=len(regions):raise ValueError('Manifest needs distinct region IDs')
    input_sidecar=Path(input_path).with_suffix('.json');identity='No paired PET metadata identity supplied'
    if input_sidecar.exists():
        input_info=json.loads(input_sidecar.read_text())
        if manifest.get('pet_metadata_sha256') and input_info.get('pet_metadata_sha256'):
            if manifest['pet_metadata_sha256']!=input_info['pet_metadata_sha256']:raise ValueError('Input and regions were prepared from different PET metadata files')
            identity='Matched PET metadata SHA256'
        if manifest.get('units') and input_info.get('units') and manifest['units']!=input_info['units']:raise ValueError('Input and regional activity concentration units differ')
        if manifest.get('time_origin') and input_info.get('time_origin') and manifest['time_origin']!=input_info['time_origin']:raise ValueError('Input and regional TimeZero origins differ')
    out=_new_output(output_dir);records=[]
    # Import once, before processing, so a missing plotting dependency gives a clear error.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for index,region in enumerate(regions):
        folder=out/f'region_{index:04d}';folder.mkdir();record=dict(id=region['id'],name=str(region.get('name',region['id'])),folder=folder.name,model=model)
        try:
            path=manifest_path.parent/region['path'];frames,y,sigma,eligible=_read_tac(path)
            if not eligible.any():raise ValueError('No positive-weight frames remain')
            fit_frames=Frames(frames.start[eligible],frames.end[eligible]);s=None if sigma is None else sigma[eligible]
            options=dict(model=model,vb=vb,fit_delay=fit_delay,input_tail=input_tail,n_starts=n_starts,seed=seed)
            fit=fit_compartment(fit_frames,inp,y[eligible],sigma=s,**options)
            # The reported full curve uses the full-fit parameters, not fold-selected parameters.
            prediction=predict(frames,inp,list(fit['parameters'].values()),model,vb=fit['vb'],delay=fit['delay'],input_tail=input_tail)
            held=(np.arange(len(y))%hold_every==hold_every-1)&eligible;train=eligible&~held;cv=dict(status='not_run',reason='Too few training or held-out observations')
            if held.any() and train.sum()>len(fit['free_parameters']):
                try:
                    trainfit=fit_compartment(Frames(frames.start[train],frames.end[train]),inp,y[train],sigma=None if sigma is None else sigma[train],**options)
                    cvpred=predict(frames,inp,list(trainfit['parameters'].values()),model,vb=trainfit['vb'],delay=trainfit['delay'],input_tail=input_tail)
                    scale=float(np.mean(np.abs(y[held])))
                    cv=dict(status='completed',held_indices=np.flatnonzero(held).tolist(),normalized_RMSE=float(np.sqrt(np.mean((cvpred[held]-y[held])**2))/scale) if scale>0 else None,training_fit_local_qc_passed=trainfit['usable_for_quantification'],prediction=cvpred.tolist())
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:cv=dict(status='failed',reason=str(exc))
            methods={}
            for name,fn in ([('patlak',patlak)] if model=='2t_irreversible' else [('logan',logan),('ma1',ma1)]):
                try:methods[name]=fn(frames,y,inp,tstar,vb=fit['vb'],delay=fit['delay'],input_tail=input_tail)
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:methods[name]=dict(status='not_available',reason=str(exc))
            fit['provenance']=file_manifest([path,input_path,manifest_path],dict(**options,tstar=tstar,hold_every=hold_every))
            (folder/'fit.json').write_text(json.dumps(fit,indent=2,allow_nan=False));write_qc_html(fit,folder/'fit.html',title=record['name'])
            np.savetxt(folder/'predictions.csv',np.c_[frames.start,frames.end,y,prediction,y-prediction,eligible.astype(int),held.astype(int)],delimiter=',',header='start,end,observed,full_fit_prediction,observed_minus_prediction,included,held_for_prediction_check',comments='')
            fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True,gridspec_kw={'height_ratios':[3,1]})
            axes[0].plot(frames.mid,y,'o',ms=3,label='Measured region');axes[0].plot(frames.mid,prediction,label='Full-data fit')
            if cv['status']=='completed':axes[0].plot(frames.mid,cv['prediction'],'--',alpha=.7,label='Fit excluding held frames');axes[0].scatter(frames.mid[held],y[held],facecolors='none',edgecolors='red',s=50,label='Held frames')
            axes[0].set_ylabel('Activity (input-matched units)');axes[0].legend(fontsize=8);axes[0].set_title(record['name'])
            axes[1].axhline(0,color='gray',lw=1);axes[1].plot(frames.mid,y-prediction,'.');axes[1].set(xlabel='Minutes',ylabel='Residual');fig.tight_layout();fig.savefig(folder/'fit.png',dpi=140);plt.close(fig)
            macro='Ki' if model=='2t_irreversible' else 'VT'
            record.update(status='completed',local_qc_passed=bool(fit['usable_for_quantification']),qc_reasons=fit['qc_reasons'],macroparameter=macro,estimate=fit[macro],raw_estimate=fit['raw_'+macro],rmse=fit['rmse'],parameters=fit['parameters'],delay_minutes=fit['delay'],excluded_frame_indices=np.flatnonzero(~eligible).tolist(),held_frame_prediction=cv,graphical=methods)
        except (ValueError,RuntimeError,OSError,KeyError,np.linalg.LinAlgError) as exc:record.update(status='failed',error=str(exc),local_qc_passed=False)
        records.append(record)
    summary=dict(schema_version=1,input_identity_check=identity,regions=len(records),completed=sum(r['status']=='completed' for r in records),failed=sum(r['status']=='failed' for r in records),local_qc_passed=sum(r['local_qc_passed'] for r in records),records=records,scope='Regional analysis and local QC; does not establish biological accuracy or full research qualification',provenance=file_manifest([manifest_path,input_path]))
    (out/'results.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    with (out/'results.csv').open('w',newline='') as file:
        fields=['id','name','status','local_qc_passed','macroparameter','estimate','raw_estimate','rmse','held_frame_normalized_RMSE','K1','k2','k3','k4','delay_minutes','qc_reasons','error'];writer=csv.DictWriter(file,fieldnames=fields);writer.writeheader()
        for r in records:
            row={key:r.get(key,'') for key in fields};row.update(r.get('parameters',{}));row['held_frame_normalized_RMSE']=r.get('held_frame_prediction',{}).get('normalized_RMSE','');row['qc_reasons']='; '.join(r.get('qc_reasons',[]));writer.writerow(row)
    e=lambda x:html.escape(str(x));rows=[]
    for r in records:
        if r['status']=='failed':rows.append(f"<tr><td>{e(r['name'])}</td><td colspan='4'>Failed: {e(r['error'])}</td></tr>");continue
        value='Withheld' if r['estimate'] is None else f"{r['estimate']:.5g}";metric=r['held_frame_prediction'].get('normalized_RMSE');cvtext='Unavailable' if metric is None else f'{metric:.1%}'
        rows.append(f"<tr><td><a href='{r['folder']}/fit.html'>{e(r['name'])}</a></td><td>{value}</td><td>{r['raw_estimate']:.5g}</td><td>{cvtext}</td><td>{e('; '.join(r['qc_reasons']) or 'Local checks passed')}</td></tr>")
    images=''.join(f"<h2>{e(r['name'])}</h2><img width='800' style='max-width:100%' src='{r['folder']}/fit.png'>" for r in records if r['status']=='completed')
    (out/'index.html').write_text("<!doctype html><meta charset='utf-8'><title>Regional PET analysis</title><style>body{font:16px system-ui;max-width:1050px;margin:40px auto;padding:0 20px}table{border-collapse:collapse;width:100%}td,th{border-bottom:1px solid #ddd;padding:12px;text-align:left}</style><h1>Regional PET analysis</h1>"+f"<p>{summary['completed']}/{summary['regions']} regions processed; {summary['local_qc_passed']} pass local checks.</p><p>Held-frame error measures prediction, not parameter accuracy. Raw estimates and flagged fits remain visible for review. Local checks do not establish biological validity.</p><p><a href='results.csv'>Download results CSV</a> · <a href='results.json'>Complete results</a></p><table><tr><th>Region</th><th>Estimate</th><th>Raw estimate</th><th>Held-frame error</th><th>Quality checks</th></tr>"+''.join(rows)+'</table>'+images)
    return summary
