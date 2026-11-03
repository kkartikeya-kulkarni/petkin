# petkin 0.2 release evidence

Reviewed and evaluated 2026-11-03. **An upgraded research candidate, not an
established “best available” or clinically validated toolkit.** The release
combines a focused literature review, independent-software comparisons,
uncertainty audits, simulation challenges and measured human PET processing.
See [the literature review](../../LITERATURE_REVIEW.md) for primary sources and
an explicit comparison-scope matrix.

## Main outcomes

- Automated checks: **28 passed in 1.84s**.
- Real PBR28 median five-fold prediction NRMSE: **6.69%** with fixed vb/delay,
  versus **4.40%** with fitted vascular fraction and delay. This is a
  descriptive **34.3%** reduction in the cohort median, not a
  proven general improvement or calibrated parameter-accuracy result.
- Median absolute kinfitr point-estimate differences: Logan
  **0.0064%**, MA1
  **0.0050%** across 20 real scans.
- NiftyPAD SRTM basis comparison on 12 synthetic ODE cases: median absolute BPnd
  difference **0.096%**, maximum
  **0.241%**.
- 72 routine simulation fits converged. The harder 24-case identifiability
  challenge still contained **2 estimates with >50% error that passed the local QC screen**.

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
| Fixed vb=.05, delay=0 | 6.69% | 20/20 | 2/20 |
| Fit vb and delay | 4.40% | 20/20 | 4/20 |

The fixed-model withheld scans are rbqc_1, ytdh_1.
The expanded-model withheld scans are mhco_1, rtvg_2, rwrd_1, rwrd_2.
A fitted nuisance parameter at a bound triggers review even if the curve score
improves. All raw estimates, fold results and reasons are retained.

A descriptive paired bootstrap resampled the ten people, keeping their repeat
scans together. The mean paired NRMSE improvement was
**2.24 percentage points**, with a subject-bootstrap
95% interval **[-0.20, 4.85] percentage points**.
The interval includes zero. This small exposed cohort does not establish
statistically conclusive superiority.

## Previously unstable case: uncertainty evidence

| Case | k4 profile closed on tested grid | Unusable bootstrap draws | Raw VT bootstrap percentile interval |
|---|---|---:|---|
| cgyu_1, fixed | True | 0% | [2.123, 2.387] |
| rbqc_1, fixed | False | 88% | [2.054, 12029.940] |
| rbqc_1, fitted_vascular_delay | True | 3% | [2.396, 3.872] |

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
| 1t | 1% | 0.11% | 0.31% | 8/8 |
| 1t | 5% | 0.46% | 1.55% | 8/8 |
| 1t | 15% | 1.87% | 2.89% | 8/8 |
| 2t | 1% | 0.19% | 0.23% | 8/8 |
| 2t | 5% | 0.57% | 0.81% | 8/8 |
| 2t | 15% | 1.20% | 3.49% | 8/8 |
| 2t_irreversible | 1% | 0.14% | 0.31% | 8/8 |
| 2t_irreversible | 5% | 0.56% | 1.21% | 8/8 |
| 2t_irreversible | 15% | 1.59% | 3.82% | 8/8 |

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

There were **15/24 raw VT errors above 50%**. Local QC flagged **13** of those; **2 large errors passed**. In total 15/24 results were withheld. No QC threshold was changed to make this challenge pass. The median raw absolute error was 82.3%.

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
| s028 | 33 | 1 | Not available |
| s029 | 32 | 1 | Not available |
| s030 | 32 | 1 | Not available |

Overlays were visually inspected for gross brain alignment. True human
registration error remains unmeasured. Automated tests check transform order
using noncommuting transforms and known physical correspondences. They do not
replace independent anatomical landmarks or segmentations.

Exploratory irreversible 2T versus Patlak comparisons use the frozen regional
TACs from v0.1. No whole blood is available; plasma timing/endpoints and decay
conventions retain the original limitations.

| Subject | Irreversible 2T Ki | Patlak Ki, t*=20 | Weighted residual lag-1 correlation |
|---|---:|---:|---:|
| s028 | 0.04205 | 0.03346 | 0.831 |
| s029 | 0.04105 | 0.03513 | 0.822 |
| s030 | 0.05609 | 0.04706 | 0.855 |

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
