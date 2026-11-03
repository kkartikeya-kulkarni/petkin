"""Portable, escaped QC reports and reproducible file manifests."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,html,importlib.metadata,json
from . import __version__


def file_manifest(paths,settings=None):
    files=[]
    for value in paths:
        path=Path(value);digest=hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
        files.append(dict(path=str(path.resolve()),bytes=path.stat().st_size,sha256=digest.hexdigest()))
    versions={}
    for name in ['numpy','scipy','SimpleITK']:
        versions[name]=importlib.metadata.version(name)
    return dict(petkin_version=__version__,created_utc=datetime.now(timezone.utc).isoformat(),
                inputs=files,settings=settings or {},dependencies=versions)


def write_qc_html(result,path,title='PET kinetic analysis'):
    """Human-readable summary. Raw estimates are always labeled as unvalidated."""
    escape=lambda x:html.escape(str(x))
    reasons=list(dict.fromkeys(result.get('qc_reasons',[])+result.get('research_screen',{}).get('reasons',[])))
    status='Review required' if reasons or not result.get('usable_for_quantification',False) else 'Passes automated local checks'
    key='BPnd' if result.get('method') in ('srtm','logan_reference') else ('Ki' if result.get('model')=='2t_irreversible' or result.get('method')=='patlak' else 'VT')
    estimate=result.get(key)
    display='Withheld — review QC' if estimate is None else f'{estimate:.5g}'
    if 'interval_components' in result:
        key='Profile interval';display=str(result['interval_components'])
        status='Closed on evaluated grid' if result['closed'] else 'Open or incomplete profile — review required'
    elif result.get('method')=='conditional_parametric_Gaussian_bootstrap':
        key='Raw bootstrap interval';display=str(result['raw_percentile_95_interval'])
        status=f"Conditional bootstrap; {result['unusable_fraction']:.1%} unusable draws"
    rows=''.join(f'<tr><th>{escape(name)}</th><td>{escape(value)}</td></tr>' for name,value in
                 [('Status',status),(key,display),('Method',result.get('model',result.get('method','unknown'))),
                  ('Vascular fraction',result.get('vb')),('Input delay (minutes)',result.get('delay',0)),
                  ('Optimizer converged',result.get('success','Not applicable')),('QC reasons',', '.join(reasons) or 'None triggered')])
    params=''.join(f'<tr><th>{escape(k)}</th><td>{escape(f"{v:.6g}")}</td></tr>' for k,v in result.get('parameters',{}).items())
    raw=result.get('raw_'+key)
    document=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escape(title)}</title><style>body{{font:16px system-ui;color:#183047;background:#f2f5f8;margin:0}}main{{max-width:900px;margin:36px auto;padding:28px;background:white;border-radius:12px}}h1{{font-size:28px}}table{{border-collapse:collapse;width:100%;margin:20px 0}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #ddd}}.note{{padding:16px;background:#fff3da}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}}a{{color:#135fa4}}</style><main>
<h1>{escape(title)}</h1><p>Research analysis · petkin {escape(__version__)}</p><table>{rows}</table>
<p class="note">Automated checks do not establish biological accuracy. Parameter bounds, uncertain input timing, missing blood measurements and image alignment can invalidate otherwise plausible curves.</p>
<h2>Fitted rate constants</h2><table>{params}</table><p>Rates use minutes. VT is mL plasma/mL tissue; Ki is mL plasma/mL tissue/min.</p>
<details><summary>Raw numerical estimate and complete audit record</summary><p>Raw {escape(key)}: {escape(raw)}. A raw value may be unsuitable for inference.</p><pre>{escape(json.dumps(result,indent=2,allow_nan=False))}</pre></details></main></html>'''
    Path(path).write_text(document)
