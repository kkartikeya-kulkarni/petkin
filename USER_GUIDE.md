# petkin 0.3.0rc1 — research qualification candidate

**Full research qualification is not established.** The fresh numerical and
synthetic-registration gates passed. A new independent PS13 cohort (20 scans,
10 people) now provides external prediction and software-comparison evidence,
but it does not complete external kinetic qualification. Measured-human PET/MRI
registration truth remains unavailable. See the
[qualification record](validation/v03/QUALIFICATION_REPORT.md),
[frozen protocol](validation/v03/PROTOCOL.json), and
[external evidence checklist](validation/v03/EXTERNAL_VALIDATION_CHECKLIST.md).
The v0.3 protocol covers regional arterial compartment estimates and rigid
registration verification. Other methods retain their earlier implementation
checks; this is not a blanket validation of every method for every tracer.

## Independent PS13 evaluation

See [external report](validation/external_ps13/REPORT.html) and
[frozen study protocol](validation/external_ps13/PROTOCOL.json). Primary held-frame
prediction error was 4.34% median and 37.59% worst case. Only 7/20 full fits passed
local QC; all 20 were withheld by the research screen. Logan and MA1 were compared
with unmodified kinfitr. All 20 automatic measured PET/MRI alignments are available
for review, with no independent accuracy claim. These subjects are now exposed
and cannot be reused as a fresh validation cohort after further development.

## Preferred research workflow

Use `analyze_research_tac` or `petkin research-fit` for the stronger conditional
screen. It requires a direct VT/Ki likelihood profile, root-refined confidence
crossings, acceptable precision, local/residual checks, and study-backed input
assumptions. Uncertain or unverified results are withheld while raw values remain
available for audit. Passing this per-result screen does not mean the full
toolbox has passed external research qualification.

```sh
petkin research-fit tac.csv input.csv result.json --model 2t --vb .05 \
  --assumptions study_assumptions.json --report quality.html
```

`examples/RESEARCH_ASSUMPTIONS_TEMPLATE.json` lists the required assumptions and
intentionally defaults to false. Record verification from study metadata and
processing records; do not assert a condition merely to obtain an estimate.
`absolute_sigma=True`/`--absolute-sigma` requires actual known frame SDs.
The lower-level `fit_compartment` / `petkin fit` remain diagnostic fits whose
local-QC flag does not replace the research screen.

Rate-only profiles now report conservative outer grid brackets. Direct
macroparameter profiles use refitted root crossings; they avoid falsely narrow
limits from coarse linear interpolation of the likelihood ratio.


A Python research toolkit for arterial compartment fitting, graphical analysis,
reference-tissue modeling and rigid PET/MRI preprocessing. This release focuses
on trustworthy failure handling, uncertainty and reproducible comparisons.
It has **not** established a global “best available” or clinical ranking.

Start with [the release evidence](validation/v02/RELEASE_REPORT.md) and
[the focused literature review](LITERATURE_REVIEW.md). Historical v0.1 results
remain under `validation/`; they are not retroactively presented as v0.2 tests.

## Install

Python 3.10+; the recorded validation environment was Python 3.13 on macOS.

```sh
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
pytest -q
```

Core dependencies are NumPy, SciPy and SimpleITK. Matplotlib/NiBabel are used by
examples. Optional `.[benchmarks]` adds the pinned NiftyPAD comparator and its
needed scikit-learn dependency. R is needed only for the external reference
comparison and original RData conversion, not for toolkit operation.

## Methods

| Family | Available methods | Reported outcomes |
|---|---|---|
| Arterial compartment models | 1T, reversible 2T, irreversible 2T | K1, k2, k3, k4 as appropriate; VT or Ki |
| Plasma graphical models | Logan, Ichise MA1, Patlak | VT or Ki, explicit t*, regression diagnostics |
| Reference tissue | SRTM, reference Logan with k2′ | R1, k2, BPnd and/or DVR |
| Uncertainty | Local covariance, conditional parametric bootstrap, profile likelihood | Standard errors, raw intervals, open/closed profile components |
| Imaging | Rigid motion correction, rigid PET/MRI registration, transform composition, ROI extraction | Images, transforms, motion/QC diagnostics |

Reference models require a biologically appropriate reference region. They are
not applied to the human PBR28 benchmark. There is no voxelwise kinetic mapping,
partial-volume correction, intraframe correction or image reconstruction here.

## Kinetic data and units

All kinetic times are **minutes**. K1 and Ki are mL plasma/mL tissue/min;
k2/k3/k4 and k2′ are min^-1; VT is mL plasma/mL tissue; BPnd/DVR/R1 are unitless.
PET and input concentrations must have matching units, clocks and decay origin.
Input CSV: `time,plasma` and optional `blood`. TAC CSV:
`start,end,activity` and optional positive `sigma`. Input must begin at zero and
cover requested times; frames must be ordered and nonoverlapping.

The supplied input must already have appropriate parent/metabolite, calibration
and dispersion processing. The library does not infer those corrections.
`petkin.bids` checks PET-BIDS timing, units and explicit image decay metadata,
converts seconds to minutes, and requires explicit missing-input endpoint
policies. It returns unresolved parent/decay verification requirements rather
than assuming blood metadata establish them.

## Fit with QC and uncertainty

```python
import numpy as np
from petkin import Frames, InputFunction, fit_compartment, profile_likelihood

frames = Frames(start_minutes, end_minutes)
input_function = InputFunction(sample_minutes, plasma_activity, whole_blood_activity)
fit = fit_compartment(frames, input_function, tissue_activity,
                      model="2t", sigma=frame_sd, absolute_sigma=True,
                      vb=0.05, n_starts=8, seed=0)

if fit["usable_for_quantification"]:
    estimate = fit["VT"]  # Passes local automated checks; further validation remains necessary.
else:
    print(fit["qc_reasons"])
    # VT is None. raw_VT is retained only for auditing an unsafe numerical estimate.

profile = profile_likelihood(frames, input_function, tissue_activity, fit,
                             "k4", np.geomspace(1e-6, 3, 51), sigma=frame_sd)
print(profile["closed"], profile["interval_components"])
```

Predictions are frame means of `(1-vb)*(C1+C2)+vb*Cb`, calculated by numerical
convolution/integration. The default 0.025-min step should be decreased to check
convergence for sharp inputs or short frames. Optimizers use log-rate
coordinates and deterministic multistarts. Default rate bounds are 1e-6–3
min^-1 (also applied numerically to K1); expose different bounds explicitly
when scientifically justified. Defaults are engineering choices, not universal
tracer priors.

`sigma` supplies frame SDs or relative weighting scales. **By default they are
relative**: residual variance is estimated. `absolute_sigma=True` requires
supplied calibrated SDs. This changes uncertainty, not the WLS optimum.

Local QC checks optimizer success, bounds, numerical rank, scaled condition
number and relative macroparameter SE. VT/Ki becomes null on a blocking flag;
`raw_VT`/`raw_Ki` preserves the numerical result. The default local relative-SE
flag is >50%, and condition flag is >1e8. These are documented engineering
screens, not validated clinical limits. The challenge benchmark includes large
errors that still pass this screen. Reported Wald intervals are approximate
and especially unreliable at bounds. Residual serial-correlation diagnostics
are descriptive, not calibrated acceptance tests.

Profiles refit nuisance parameters at every grid value. Missing confidence
crossings are open; grid endpoints are never reported as closed confidence
limits. Chi-square calibration is asymptotic. Bootstrap intervals retain
boundary/invalid raw estimates and report their frequency. Both procedures
condition on the model, input and preprocessing: they do not propagate blood,
reference-region, registration or model-selection uncertainty.

## Delay and vascular fraction

```python
fit = fit_compartment(frames, input_function, tissue_activity,
                      model="2t", fit_vb=True, fit_delay=True,
                      input_tail="hold")
```

A positive delay evaluates Cp(t-delay) and Cb(t-delay): the measured input
arrives later in tissue. A negative delay may require input samples beyond
scan end. `input_tail="error"` is the default; `"hold"` explicitly assumes a
constant terminal input. Fitting vb requires whole blood. Default vascular
bounds are 0–0.3 and delay bounds are -0.5–0.5 min. More parameters can increase
overfitting or boundary contact; this model remains optional. It improves some
PBR28 predictions and worsens others in the recorded development benchmark.

## Graphical analysis

```python
from petkin import logan, ma1, patlak, graphical_stability
result = ma1(frames, tissue_activity, input_function, tstar=40, vb=.05)
sensitivity = graphical_stability(frames, tissue_activity, input_function,
                                  [30, 40, 50], method="logan", vb=.05)
```

Graphical tissue integrals use piecewise-constant frame means; piecewise-linear
plasma/blood integrals are evaluated analytically. Tissue must be contiguous
from zero. Missing initial tissue is never filled implicitly. Logan excludes
nonpositive corrected late tissue; MA1 avoids dividing by tissue. Both can
remain biased by kinetics, noise and scan duration. MA1's denominator is checked
for instability. Patlak is for effectively irreversible kinetics after an
appropriate t*. t* sensitivity reports every requested setting without selecting
the most favorable result. Zero regression weights exclude frames from the fit,
not from the measured tissue integral.

```python
from petkin import ReferenceCurve, fit_srtm, logan_reference
reference = ReferenceCurve(reference_sample_minutes, reference_activity)
srtm = fit_srtm(frames, reference, tissue_activity)
ref_logan = logan_reference(frames, reference, tissue_activity,
                           tstar=30, k2prime=.15)
```

ReferenceCurve represents a sampled nonnegative reference concentration curve,
not an implicit interpolation of arbitrary frame means. Reference measurement
uncertainty is ignored. Validate the tracer/reference assumptions independently.

## Command line and portable reports

```sh
petkin fit tac.csv input.csv fit.json --model 2t --vb .05 --report qc.html
petkin fit tac.csv input.csv fit_delay.json --fit-vb --fit-delay --input-tail hold
petkin ma1 tac.csv input.csv ma1.json --tstar 40 --vb .05
petkin patlak tac.csv input.csv patlak.json --tstar 20
petkin profile tac.csv input.csv profile.json --fit-result fit.json \
  --parameter k4 --minimum .000001 --maximum 3 --points 51
petkin bootstrap tac.csv input.csv bootstrap.json --fit-result fit.json --samples 100
```

Use `--spacing linear` to profile delay or vb across zero. CLI JSON includes
input hashes, dependency versions, settings and timestamps. HTML QC reports
visibly withhold flagged estimates; raw results remain in a labeled audit
section. Source CSV time columns always use minutes.

## Images and transforms

```sh
petkin preprocess dynamic_pet.nii.gz T1w.nii.gz pet.json processed --reference 20
```

PET-BIDS frame start/duration arrays are required and checked against image
frame count. Motion estimation uses rigid correlation registration. PET/MRI
uses Mattes mutual information and deterministic geometry/moments
initializations; selecting the best objective is not anatomical accuracy
validation. Optional API registration masks must match image geometry.

Blank frames raise by default. `--preserve-blank` preserves and flags them with
identity transforms because motion is unmeasurable. Exclude them from fitting.
Motion displacement QC uses physical landmarks; it is not Power framewise
 displacement and does not establish accuracy.

The CLI uses `preprocess_dynamic`: temporary corrected images support estimation,
then original frames are resampled directly onto MRI with composed transforms.
**One final interpolation per frame.** The older `motion_correct` and
`pet_to_mri` APIs remain available; using them sequentially retains their two-stage
interpolation behavior. SimpleITK physical coordinates are **LPS millimeters**.
Transforms map fixed points to moving points:
`T_mri_to_original_frame = T_reference_to_frame(T_mri_to_reference(point))`.
The composition order is covered by a noncommuting-transform test.

`apply_composed_transforms` accepts externally estimated transforms after their
conventions have been explicitly converted. FSL matrices are not automatically
interpreted as ITK/LPS transforms. `landmark_registration_error` measures TRE
only when independent corresponding physical landmarks are supplied. Human
landmark accuracy remains unvalidated in this release.

`roi_tacs` requires integer labels on the exact PET spatial grid. Resample labels
with nearest-neighbor interpolation. Do not use linear interpolation for labels.

## Reproduce evidence

See `validation/v02/RELEASE_REPORT.md` for the frozen protocol and results.
The example scripts provide:

- 72 independent-ODE simulation fits across three models/noise levels.
- 24 harder short-scan, slow-dissociation, Poisson-count cases.
- Twenty real PBR28 scan comparisons, including five-fold frame prediction.
- Three real FDG method comparisons and historical MRI outputs.
- Profiles and 100-draw bootstraps for three illustrative model/scan cases.
- Unmodified kinfitr source-function and installed NiftyPAD comparisons.

RData/source and image fetch scripts pin or verify available hashes. Raw data
stay in the user's working data directories. Data/model choices and limitations
are described in the reports; these results are development evidence, not
preregistered multicenter validation.
