"""Evaluate the frozen acceptance gates; blocked evidence never becomes passed."""
from pathlib import Path
import hashlib,html,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from petkin import qualification_status


def run(folder):
    folder=Path(folder);load=lambda n:json.loads((folder/n).read_text())
    protocol=load('PROTOCOL.json');expected=(folder/'PROTOCOL.sha256').read_text().split()[0]
    if hashlib.sha256((folder/'PROTOCOL.json').read_bytes()).hexdigest()!=expected:raise ValueError('Protocol checksum mismatch')
    stats=load('holdout_summary.json');registration=load('registration_holdout.json')['summary'];external=load('EXTERNAL_EVIDENCE_STATUS.json')
    logs=[folder/'test_results.txt',folder/'wheel_test_results.txt']
    root=folder.parents[1];freeze=load('IMPLEMENTATION_FREEZE.json');source_matches=[]
    for record in freeze['files']:
        content=(root/record['path']).read_bytes();expected_hash=record['sha256']
        if record['path']=='src/petkin/__init__.py':content=content.replace(b"__version__='0.3.0rc1'",b"__version__='0.2.0'")
        if record['path']=='examples/qualification_registration.py':expected_hash=load('REGISTRATION_GENERATOR_AMENDMENT.json')['sha256']
        source_matches.append(hashlib.sha256(content).hexdigest()==expected_hash)
    freeze_matches=all(source_matches)
    sw=freeze_matches and all(p.exists() and ' passed' in p.read_text() and 'failed' not in p.read_text() and 'error' not in p.read_text().lower() for p in logs)
    routine,stress=stats['routine'],stats['stress'];limits={g['id']:g.get('limits',{}) for g in protocol['required_gates']}
    numeric=routine['retention']>=limits['NUM-01']['min_retention'] and routine['median_absolute_error'] is not None and routine['median_absolute_error']<=limits['NUM-01']['max_median_error'] and routine['q95_absolute_error']<=limits['NUM-01']['max_q95_error']
    risk=routine['one_sided_95_large_error_risk_upper']<=limits['NUM-02']['max_binomial_upper']
    robust=stress['large_errors_released']<=limits['ROB-01']['max_large_errors'] and stress['retention']>=limits['ROB-01']['min_retention']
    geometric=registration['mean_landmark_error_mm']<=limits['REG-01']['max_mean_mm'] and registration['maximum_landmark_error_mm']<=limits['REG-01']['max_error_mm']
    gates=[dict(id='SW-01',status='passed' if sw else ('failed' if all(p.exists() for p in logs) else 'not_run'),frozen_source_matches=freeze_matches,evidence=[str(p.name) for p in logs]),
           dict(id='NUM-01',status='passed' if numeric else 'failed',observed=routine),
           dict(id='NUM-02',status='passed' if risk else 'failed',observed=routine['one_sided_95_large_error_risk_upper']),
           dict(id='ROB-01',status='passed' if robust else 'failed',observed=stress),
           dict(id='REG-01',status='passed' if geometric else 'failed',observed=registration),
           dict(id='EXT-01',status='blocked',evidence=external['external_kinetic']),
           dict(id='REG-02',status='blocked',evidence=external['measured_registration'])]
    report=qualification_status(protocol,gates);report['protocol_sha256']=expected
    report['scope']=protocol['scope'];report['status_caveat']='Local software-generated qualification record, not an independent certification or public preregistration'
    (folder/'qualification.json').write_text(json.dumps(report,indent=2))
    known=load('known_failure_regression.json');results=load('holdout_results.json')
    text=f'''# Research qualification record

## Overall decision: NOT QUALIFIED for the full declared scope

The requested use is **regional brain-PET research analysis**. Numerical and
synthetic-registration verification were evaluated against a locally frozen
protocol. Untouched external kinetic validation and measured-human PET/MRI
registration truth remain unavailable. Neither a passing simulation nor a
new release label substitutes for those mandatory gates.

Protocol: `{protocol['protocol_id']}`. SHA256: `{expected}`.
The criteria were frozen before the new holdout outcomes were generated.
This is a local prospective freeze, not a public preregistration or independent
certification. Implementation/generator hashes and a documented setup-only
registration-generator correction are included.

| Gate | Result | Requirement |
|---|---|---|
'''
    for gate in gates:
        requirement=next(g for g in protocol['required_gates'] if g['id']==gate['id'])
        text+=f"| {gate['id']} | {gate['status'].upper()} | {requirement['title']} |\n"
    text+=f'''
## Fresh kinetic tests

The new seed was {protocol['holdout']['seed']}. All {len(results)} cases and
outcomes are retained in `holdout_data.json` and `holdout_results.json`.
Ninety routine cases span 1T, reversible 2T and irreversible 2T with known
processed inputs, timing, vascular fraction and Gaussian SDs. Forty hard
reversible cases vary scan duration, slow dissociation and count sensitivity,
with C11-like decay and Poisson/background-subtracted observations. Ground
truth comes from independent adaptive ODE integration of frame activity.
These are numerical/model tests, not independent human biology.

| Set | Returned estimates | Retention | Median absolute macro error among returned | 95th percentile error | Returned errors >50% |
|---|---:|---:|---:|---:|---:|
'''
    for name,value in [('Routine',routine),('Hard challenge',stress)]:
        text+=f"| {name} | {value['released']}/{value['total']} | {value['retention']:.1%} | {value['median_absolute_error']:.2%} | {value['q95_absolute_error']:.2%} | {value['large_errors_released']} |\n"
    text+=f'''
The stricter screen deliberately withholds uncertain results. Rejection is not
successful recovery of the missing kinetic information. Retention floors in
the frozen criteria prevent passing merely by rejecting everything.

The one-sided 95% binomial upper summary for >50% error among returned routine
estimates is **{routine['one_sided_95_large_error_risk_upper']:.2%}**. In the hard
set it is **{stress['one_sided_95_large_error_risk_upper']:.2%}** because only
{stress['released']} estimates were returned. Zero observed failures in that
small subset does not establish a low real-world failure risk. Binomial
interpretation assumes independent/exchangeable returned trials and is limited
to the tested simulation case mix, not tracers or sites generally. Results are
not pooled across the routine and stress sets to hide this uncertainty.

## Changes to inference and confidence limits

The research entry point now profiles **VT or Ki itself**, refitting nuisance
rates under the exact algebraic macroparameter constraint. Confidence crossings
are refined by new nuisance optimizations and a root solver. This fixes the
risk of falsely narrow limits from linearly interpolating a coarse likelihood
grid. Rate-only profiles now return conservative outer crossing brackets.

A research result requires local checks, a single closed macro profile, adequate
profile precision, no strong residual-structure flag, and explicit study-backed
input assumptions. Raw values remain for audit. The two previously exposed
v0.2 misses are regression cases, not fresh validation evidence:

'''
    for case in known:
        assessment=case['result']['research_screen']
        text+=f"- True k4={case['truth_k4']}: returned={assessment['passed']}; reasons: {', '.join(assessment['reasons'])}.\n"
    text+='''
A closed profile remains conditional on the selected model, input, noise and
preprocessing. It is not proof of biological accuracy or complete uncertainty
propagation. The public low-level fitting API remains a diagnostic numerical
fit; use the stricter research entry point for the defined screen.

## Synthetic registration verification

'''
    text+=f"Twelve fresh known rigid poses on three MRI-derived synthetic PET volumes yielded mean landmark error **{registration['mean_landmark_error_mm']:.3f} mm**, maximum **{registration['maximum_landmark_error_mm']:.3f} mm**. These satisfy the predefined numerical registration gate. The synthetic modality was produced by an intensity transform, blur and Poisson noise; it is deliberately not called measured human PET ground truth.\n"
    text+='''
The first registration generator attempt stopped before any case evaluation
because tiny negative blur values caused invalid fractional-power/count inputs.
A nonnegative synthetic-intensity correction was recorded before evaluating
any pose. No seed, optimizer or acceptance threshold was changed. The initial
freeze and the amendment are both retained.

## Why full research qualification remains blocked

**External kinetic evidence (EXT-01):** Existing PBR28/FDG cohorts were already
examined during development. They cannot become untouched external evidence
by renaming a split. A relevant independent four-tracer validation study says
its underlying participant data are institutionally restricted:
[study data availability](https://link.springer.com/article/10.1186/s13550-020-00664-8).
No access request or communication was sent. No outside reviewer was invented.

**Human PET/MRI truth (REG-02):** The public
[RIRE training archive](https://rire.insight-journal.org/download_training_data)
provides a potential source of independent transforms, but the PET, MRI and
transform URLs all returned HTTP 403 here. Alternate public gateways also
returned 403. No independently annotated dataset was available from the current sources. Existing fitted registrations and visual
overlays are not independent truth and were not used as such.

To complete these gates, obtain an independent cohort with documented input
processing and a frozen analysis reviewed independently, and at least ten
measured PET/MRI pairs with independent transforms/landmarks and an error
criterion agreed before evaluation. The adapters and gate machinery can record
that evidence, but cannot manufacture it. The mandatory gate states remain
BLOCKED and the overall research-qualified flag remains false.

## Reproduction

Use the frozen protocol and the included example generators. The holdout script
refuses to overwrite the saved dataset. Its seed and source/data hashes are
recorded. Tests run against both source and the built wheel. Only the recorded
macOS/Python environment is within the current software-qualification scope.
The source freeze separates numerical implementation from the later package
version metadata. No acceptance threshold was relaxed after observing results.
'''
    (folder/'QUALIFICATION_REPORT.md').write_text(text)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    labels=['Routine','Hard'];returned=[routine['released'],stress['released']];withheld=[routine['total']-returned[0],stress['total']-returned[1]]
    axes[0].bar(labels,returned,label='Estimate returned',color='#167d8d');axes[0].bar(labels,withheld,bottom=returned,label='Withheld for review',color='#c7d0db')
    axes[0].set(ylabel='Cases',title='Fresh tests: retention and rejection');axes[0].legend(fontsize=8)
    for kind,color in [('routine','#167d8d'),('stress','#bf6d2f')]:
        group=[r for r in results if r['kind']==kind and r['released']];values=sorted(abs(r['raw_relative_error'])*100 for r in group)
        axes[1].plot(np.arange(1,len(values)+1)/len(values),values,'.-',label=kind,color=color)
    axes[1].axhline(50,color='black',linestyle='--',label='50% large-error threshold')
    axes[1].set(xlabel='Fraction of returned estimates',ylabel='Absolute macro error (%)',title='Returned-estimate errors; rejected cases retained in JSON');axes[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(folder/'qualification_evidence.png',dpi=150);plt.close(fig)
    import base64
    picture=base64.b64encode((folder/'qualification_evidence.png').read_bytes()).decode()
    core_passed=all(g['status']=='passed' for g in gates if g['id'] not in ['EXT-01','REG-02'])
    core_message='Numerical and synthetic-registration criteria passed.' if core_passed else 'Numerical/software criteria are incomplete or failed.'
    rows=''.join(f"<tr><td>{html.escape(g['id'])}</td><td>{html.escape(g['status'].upper())}</td></tr>" for g in gates)
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research qualification: incomplete</title><style>body{{font:16px/1.6 system-ui;background:#edf2f6;color:#183047}}main{{max-width:1050px;margin:30px auto;background:white;padding:30px;border-radius:12px}}.warning{{background:#fff1d6;padding:18px}}table{{border-collapse:collapse;width:100%}}td{{padding:10px;border-bottom:1px solid #ddd}}img{{width:100%}}pre{{white-space:pre-wrap;font:14px/1.7 system-ui;overflow-wrap:anywhere}}</style><main><h1>Research qualification is incomplete</h1><p class="warning">{html.escape(core_message)} Mandatory independent kinetic and measured-human registration evidence is still blocked. This toolbox is not labeled fully validated.</p><table>{rows}</table><img alt="Retention and error results" src="data:image/png;base64,{picture}"><details><summary>Full evidence, assumptions and blockers</summary><pre>{html.escape(text)}</pre></details></main></html>'''
    (folder/'QUALIFICATION_REPORT.html').write_text(page)
    print(json.dumps(dict(overall=report['overall_status'],gates={g['id']:g['status'] for g in gates}),indent=2))
if __name__=='__main__':run(sys.argv[1])
