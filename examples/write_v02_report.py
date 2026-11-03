"""Generate release evidence and scientific plots from saved benchmark outputs."""
from pathlib import Path
import base64,html,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run(folder):
    folder=Path(folder)
    read=lambda name:json.loads((folder/name).read_text())
    real=read('real_v02_summary.json');comparison=read('cross_software_agreement.json')['summary']
    nifty=read('niftypad_agreement.json')['summary'];simulation=read('simulation_v02.json')
    challenge=read('identifiability_challenge.json');uncertainty=read('uncertainty_audit.json')
    paired=read('paired_subject_comparison.json');imaging=read('imaging_summary.json');fdg=read('fdg_v02.json')
    tests=(folder/'test_results.txt').read_text().strip().splitlines()[-1]
    fixed=real['fixed']['median_cross_validated_nrmse'];advanced=real['fitted_vascular_delay']['median_cross_validated_nrmse']
    text=f'''# petkin 0.2 release evidence

Reviewed and evaluated 2026-11-03. **An upgraded research candidate, not an
established “best available” or clinically validated toolkit.** The release
combines a focused literature review, independent-software comparisons,
uncertainty audits, simulation challenges and measured human PET processing.
See [the literature review](../../LITERATURE_REVIEW.md) for primary sources and
an explicit comparison-scope matrix.

## Main outcomes

- Automated checks: **{tests}**.
- Real PBR28 median five-fold prediction NRMSE: **{fixed:.2%}** with fixed vb/delay,
  versus **{advanced:.2%}** with fitted vascular fraction and delay. This is a
  descriptive **{1-advanced/fixed:.1%}** reduction in the cohort median, not a
  proven general improvement or calibrated parameter-accuracy result.
- Median absolute kinfitr point-estimate differences: Logan
  **{comparison['logan']['median_absolute_difference']:.4%}**, MA1
  **{comparison['ma1']['median_absolute_difference']:.4%}** across 20 real scans.
- NiftyPAD SRTM basis comparison on 12 synthetic ODE cases: median absolute BPnd
  difference **{nifty['median_absolute_BPnd_difference']:.3%}**, maximum
  **{nifty['maximum_absolute_BPnd_difference']:.3%}**.
- 72 routine simulation fits converged. The harder 24-case identifiability
  challenge still contained **{challenge['summary']['large_errors_not_flagged']} estimates with >50% error that passed the local QC screen**.

## What changed from v0.1

Log-rate deterministic multistarts; optional shared plasma/whole-blood delay
and vascular fitting; explicit known-SD versus relative weighting; local macro
uncertainty and residual diagnostics; failure reasons and null VT/Ki on blocking
checks; conditional bootstrap and nuisance-refitted profile likelihood; exact
piecewise-linear input integrals; MA1 and Patlak; SRTM and reference Logan;
strict BIDS timing/unit helpers; portable HTML QC and file-hash provenance;
geometry/moments rigid registration candidates, physical displacement/landmark
QC, and composed transforms for one final interpolation from input frames.

The default nonlinear model still fixes nuisance parameters unless requested.
Adding parameters is not automatically preferable. Every model/data assumption
must be chosen for the tracer and available measurements.

## Real PBR28 development benchmark

Twenty scans from ten paired test/retest participants, frontal cortex, from
[kinfitr's real PBR28 data](https://github.com/mathesong/kinfitr). The same cohort
had already been examined in v0.1 and is **not a new external held-out cohort**.
All nonlinear models use frame-averaged predictions and the source positive
frame weights, treated as relative. Full fits use eight deterministic starts;
five-fold fits use four. Every fifth eligible tissue frame is omitted in turn;
the full arterial input remains known. No settings were selected automatically
using the test-fold score.

| Model | Median five-fold NRMSE | Optimizer converged | Estimates withheld by local QC |
|---|---:|---:|---:|
| Fixed vb=.05, delay=0 | {fixed:.2%} | {real['fixed']['converged']}/20 | {len(real['fixed']['unusable_scans'])}/20 |
| Fit vb and delay | {advanced:.2%} | {real['fitted_vascular_delay']['converged']}/20 | {len(real['fitted_vascular_delay']['unusable_scans'])}/20 |

The fixed-model withheld scans are {', '.join(real['fixed']['unusable_scans'])}.
The expanded-model withheld scans are {', '.join(real['fitted_vascular_delay']['unusable_scans'])}.
A fitted nuisance parameter at a bound triggers review even if the curve score
improves. All raw estimates, fold results and reasons are retained.

A descriptive paired bootstrap resampled the ten people, keeping their repeat
scans together. The mean paired NRMSE improvement was
**{paired['mean_difference']*100:.2f} percentage points**, with a subject-bootstrap
95% interval **[{paired['subject_bootstrap_95_interval'][0]*100:.2f}, {paired['subject_bootstrap_95_interval'][1]*100:.2f}] percentage points**.
The interval includes zero. This small exposed cohort does not establish
statistically conclusive superiority.

## Previously unstable case: uncertainty evidence

| Case | k4 profile closed on tested grid | Unusable bootstrap draws | Raw VT bootstrap percentile interval |
|---|---|---:|---|
'''
    for r in uncertainty:
        b=r['bootstrap'];interval=b['raw_percentile_95_interval']
        text+=f"| {r['scan']}, {r['variant']} | {r['profile_k4']['closed']} | {b['unusable_fraction']:.0%} | [{interval[0]:.3f}, {interval[1]:.3f}] |\n"
    text+='''
Each bootstrap uses 100 Gaussian conditional draws with seed 20261103; intervals
retain invalid raw estimates. The fixed-vb/delay `rbqc_1` solution remains at
the k4 lower bound and is withheld. Its extreme VT is no longer presented as
usable. The expanded model produces a finite estimate and a closed k4 profile
under its assumptions. This supports model/input-timing sensitivity, not proof
that the new biological VT is correct. k4 profiles are not VT profiles, and the
local/conditional intervals exclude input, registration and model uncertainty.

## Independent software execution

The unmodified kinfitr source functions were retrieved at the commit recorded in
`kinfitr_source_manifest.json`. The retrieved DESCRIPTION reported 0.9.7.
Logan/MA1 use matched t*=40 min, vb=.05, unweighted regression and
frame-duration tissue integration. Small input/blood interpolation differences
remain. The full kinfitr application was not installed; its sourced modeling
functions ran with their needed R dependencies. Unavailable standard errors
from the upstream helper were not used. This compares point estimates only.

Three sourced kinfitr nonlinear 2TCM runs are also saved. Those use midpoint
observations rather than petkin's frame means, so their differences are
diagnostic comparisons, not an exact-equivalence test. The independent
kinfitr run also produces an extreme boundary VT for `rbqc_1` with fixed
nuisance parameters; the difficulty is not unique to petkin.

NiftyPAD 2.0.0's released SRTM basis/reference Logan functions ran without model
changes on twelve independent ODE cases. SRTM uses 1000 basis rates in
[5e-5,.02] s^-1. Discrete second-grid/basis predictions differ from petkin's
minute-grid frame-averaged nonlinear predictions. NiftyPAD's nonlinear SRTM
entry point was incompatible with heterogeneous curve-fit inputs in the
installed NumPy/SciPy; the working released basis method was used instead.
No human reference-region validation is claimed; PBR28 was not given an
inappropriate reference region.

PETPrep, PMOD and other packages were reviewed through primary studies/docs,
but were not executed in this environment. This report does not claim a
comprehensive head-to-head comparison with all leading tools.

## Routine independent-ODE simulations

Twenty-four parameter/input combinations, three model families, and independent
Gaussian noise at 1%, 5%, 15% produce 72 fits. Seed 20261103. Ground truth comes
from an adaptive ODE with an integrated observation state, independent of the
fitting convolution. The vascular fraction and delay are known and fixed;
SDs are known. The activity-noise floor is 5% of peak. These cases are generally
well-conditioned and are not a scanner/reconstruction simulation.

| Model | Noise scale | Median absolute macro error | Maximum absolute macro error | Usable / total |
|---|---:|---:|---:|---:|
'''
    for r in simulation['summary']:
        text+=f"| {r['model']} | {r['noise']:.0%} | {r['median_absolute_macro_error']:.2%} | {r['max_absolute_macro_error']:.2%} | {r['usable']}/{r['n']} |\n"
    text+='''
Macro outcomes are VT for reversible models and Ki for irreversible models.
Graphical-method results and approximate local interval coverage are also
saved, including when MA1 or Patlak performs worse. Eight cases per cell are
far too few to establish calibrated 95% uncertainty coverage.

## Hard challenge and important remaining failures

The challenge varies k4 over .0001, .001, .01, .05 min^-1; scan duration over
30,60,120 min; and two count sensitivities. Twenty-four cases use seed 1742,
C11-like decay, Poisson counts and background subtraction. WLS uses known
expected SDs. This is a count-observation stress test, not a full PET
reconstruction/noise simulation.

'''
    c=challenge['summary']
    caught=sum(abs(r['raw_relative_error'])>.5 and not r['fit']['usable_for_quantification'] for r in challenge['cases'])
    text+=f"There were **{c['raw_errors_over_50_percent']}/24 raw VT errors above 50%**. Local QC flagged **{caught}** of those; **{c['large_errors_not_flagged']} large errors passed**. In total {c['flagged']}/24 results were withheld. No QC threshold was changed to make this challenge pass. The median raw absolute error was {c['median_absolute_raw_error']:.1%}.\n"
    text+='''
This is a material limitation. Short/noisy scans may not identify equilibrium
binding. Optimizer convergence, favorable simulations and good curve prediction
cannot create missing kinetic information. Profiles/bootstraps add evidence
but are not guarantees under noise or misspecification. The passing local
screen must not be interpreted as a clinical acceptance certificate.

## Real FDG and imaging

The same three source FDG subjects were processed through the new composed
pipeline at >=4-mm resolution. Benchmark input frames are first downsampled;
one final interpolation follows from those inputs. This is not a claim of only
one interpolation from native raw resolution in the benchmark. The core API
can take original-resolution images directly. Mean images, composed transforms,
candidate metrics and overlays are saved under `imaging/`.

| Subject | Frames | Final interpolations from test input | Human landmark truth |
|---|---:|---:|---|
'''
    for r in imaging:text+=f"| {r['subject']} | {r['frames']} | {r['final_interpolations_per_frame']} | Not available |\n"
    text+='''
Overlays were visually inspected for gross brain alignment. True human
registration error remains unmeasured. Automated tests check transform order
using noncommuting transforms and known physical correspondences. They do not
replace independent anatomical landmarks or segmentations.

Exploratory irreversible 2T versus Patlak comparisons use the frozen regional
TACs from v0.1. No whole blood is available; plasma timing/endpoints and decay
conventions retain the original limitations.

| Subject | Irreversible 2T Ki | Patlak Ki, t*=20 | Weighted residual lag-1 correlation |
|---|---:|---:|---:|
'''
    for r in fdg:
        fit=r['fit'];pt=r['patlak']['results'][1]
        text+=f"| {r['subject']} | {fit['raw_Ki']:.5f} | {pt['raw_Ki']:.5f} | {fit['residual_diagnostics']['lag1_correlation']:.3f} |\n"
    text+='''
Ki units are mL plasma/mL tissue/min. Method disagreement and strong serial
residual structure indicate unresolved model/input/noise questions. No true Ki
is available. The earlier late-block prediction errors of 16–25% remain an
unresolved limitation; they were not tuned away using test data.

## Reproduction and remaining evidence gaps

Use the commands in `README.md` and `examples/REPRODUCE.md`. Source hashes,
R session, Python environment, raw fit results and settings accompany the report.
The optional comparator dependencies are separate from the core toolkit.
The paired bootstrap is descriptive; the real cohort was exposed during
development. No preregistered multicenter study, independent human-registration
accuracy, clinical threshold calibration, or full uncertainty propagation was
performed. Voxelwise performance, partial-volume correction and the broader
features of established pipelines are also outside this release.

The defensible conclusion is **substantially improved, auditable research
software with strong implementation-agreement evidence and clearly exposed
failure modes**. A best-in-class ranking remains unestablished.
'''
    (folder/'RELEASE_REPORT.md').write_text(text)
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    axes[0,0].bar(['Fixed vb/delay','Fitted vb/delay'],[100*fixed,100*advanced],color=['#64748b','#167d8d'])
    axes[0,0].set(ylabel='Median held-out NRMSE (%)',title='20 real PBR28 scans: descriptive CV comparison')
    axes[0,0].text(.03,.95,'Paired subject-bootstrap interval includes zero',transform=axes[0,0].transAxes,va='top',fontsize=9)
    refs={r['scan']:r for r in __import__('csv').DictReader((folder/'kinfitr_graphical.csv').open())}
    own=read('real_v02.json');x=[];y=[]
    for r in own:x.append(float(refs[r['scan']]['ma1_VT']));y.append(r['graphical']['ma1']['results'][1]['raw_VT'])
    axes[0,1].scatter(x,y,color='#167d8d');limits=[min(x),max(x)];axes[0,1].plot(limits,limits,'k--')
    axes[0,1].set(xlabel='kinfitr MA1 VT',ylabel='petkin MA1 VT',title='Independent source-function agreement')
    for r in uncertainty:
        grid=r['profile_k4']['grid'];label=f"{r['scan']}, {'fixed' if r['variant']=='fixed' else 'fitted nuisance'}"
        axes[1,0].semilogx([v['value'] for v in grid],[v['likelihood_ratio'] for v in grid],label=label)
    axes[1,0].axhline(uncertainty[0]['profile_k4']['threshold'],color='black',linestyle='--',label='Asymptotic 95% threshold')
    axes[1,0].set(xlabel='k4 (min⁻¹)',ylabel='Profile likelihood ratio',ylim=(0,10),title='Open and closed profiles (y-axis limited to 10)');axes[1,0].legend(fontsize=8)
    group=[r for r in simulation['summary'] if r['model']=='2t'];x=np.arange(3)
    for i,(key,label) in enumerate([('median_absolute_macro_error','2T fit'),('logan_median_absolute_error','Logan'),('ma1_median_absolute_error','MA1')]):
        axes[1,1].bar(x+(i-1)*.25,[100*r[key] for r in group],width=.25,label=label)
    axes[1,1].set(xticks=x,xticklabels=['1% noise','5% noise','15% noise'],ylabel='Median absolute VT error (%)',title='Routine simulations; hard cases reported separately');axes[1,1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(folder/'release_evidence.png',dpi=150);plt.close(fig)
    encoded=base64.b64encode((folder/'release_evidence.png').read_bytes()).decode()
    summary=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>petkin 0.2 evidence</title><style>body{{font:16px system-ui;line-height:1.6;color:#183047;background:#edf2f6;margin:0}}main{{max-width:1100px;margin:32px auto;padding:30px;background:white;border-radius:14px}}img{{width:100%;height:auto}}.note{{background:#fff1d6;padding:18px;border-radius:8px}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font:14px/1.7 system-ui}}a{{color:#086c82}}</style><main><h1>petkin 0.2 · research release</h1><p>Literature-informed methods, uncertainty and independent implementation comparisons.</p><p class="note"><b>Promising implementation evidence, with material limitations.</b> The difficult challenge contains two large errors that pass local screening. A global or clinical superiority claim is not established.</p><img alt="Release benchmark results and uncertainty profiles" src="data:image/png;base64,{encoded}"><details><summary>Complete report, assumptions and remaining gaps</summary><pre>{html.escape(text)}</pre></details></main></html>'''
    (folder/'RELEASE_REPORT.html').write_text(summary)

if __name__=='__main__':run(sys.argv[1])
