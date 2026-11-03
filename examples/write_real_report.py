"""Create a factual report from completed, saved human-data evaluations."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1])
fdg=json.loads((root/'fdg_stress/stress_results.json').read_text())
pbr=json.loads((root/'pbr28/pbr28_summary.json').read_text())
pbrrows=json.loads((root/'pbr28/pbr28_results.json').read_text())
source=json.loads((root/'pbr28/source.json').read_text())
text='''# Real PET data testing — 2026-11-03

The toolkit was exercised on **23 real PET scans across two public datasets**:
3 FDG scans with PET/MRI images, and 20 PBR28 scans from 10 people with regional
activity and blood curves. This extends the original two-subject demonstration.
Results are mixed: curve prediction is encouraging, but temporal extrapolation
is weaker and one reversible fit produces an unusable distribution volume.
A low prediction error does not establish reliable kinetic parameters.

## FDG PET/MRI: three subjects

Source: [OpenNeuro ds004513 v1.0.6](https://openneuro.org/datasets/ds004513/versions/1.0.6),
Castrillon et al., *The energetic costs of the human connectome*, CC0.
s028 and s029 were previously examined; s030 was a new evaluation subject.
No optimizer settings or rate bounds were tuned on these results.

Each subject underwent 4-mm resampling, rigid framewise motion correction,
PET-to-MRI registration and regional activity extraction. Blank frames were
preserved and flagged without inventing motion estimates, and excluded from
fitting. Overlays were inspected for gross brain alignment; they do not provide
landmark-based human registration accuracy.

Five interleaved frame folds were held out in turn. For each fold, a threshold
ROI was defined using only training PET frames (35% of the maximum
training-duration-weighted mean). The plasma input was supplied for the full
scan, and an irreversible two-tissue model was fitted to training tissue frames.
Scores are conditional on known plasma and frozen image preprocessing; the
motion reference and image transforms used all frames. This is not a fully
independent validation of the entire preprocessing pipeline.

NRMSE = prediction RMSE / mean held-out activity. A second challenge fitted
only the first 27 frames (17–20 minutes, depending on source timing) and
predicted the final 5–6 frames through minute 45, with the full plasma input known.

| Subject | Five-fold pooled NRMSE | Late-block NRMSE | Held-out late interval (min) | Fold Ki span / median Ki |
|---|---:|---:|---|---:|
'''
for r in fdg:
    text+=f"| {r['subject']} | {r['pooled_cross_validated_nrmse']:.2%} | {r['late_block_nrmse']:.2%} | {r['late_test_start_min']:g}–{r['late_test_end_min']:g} | {r['fold_Ki_span_relative_to_median']:.2%} |\n"
text+='''
The Ki span is a descriptive sensitivity measure, not a confidence interval.
Long-range predictions are substantially weaker than interleaved-frame
predictions. Omitted frames often have nearby training frames, making
interleaved validation an easier task than predicting the late scan.

Sensitivity checks refitted each curve at uptake thresholds 25%, 35%, 45%,
halved the integration step from 0.025 to 0.0125 min, and replaced the zero
initial plasma assumption with the first supplied plasma value. Changing the
ROI changes the tissue mixture and need not represent numerical instability.

| Subject | Ki, threshold 25% | Ki, 35% | Ki, 45% | Ki change, half step | Ki change, initial plasma hold |
|---|---:|---:|---:|---:|---:|
'''
for r in fdg:
    ks=[s['Ki'] for s in r['roi_sensitivity']]
    text+=f"| {r['subject']} | {ks[0]:.5f} | {ks[1]:.5f} | {ks[2]:.5f} | {r['step_halving_Ki_relative_change']:.3%} | {r['initial_plasma_hold_Ki_relative_change']:.2%} |\n"
text+='''
Ki units are mL plasma / mL tissue / min. These are assumption-dependent
estimates without known biological truth. FDG plasma was supplied at frame
midpoints, not as densely measured bolus samples. The example assumes zero
plasma at injection, linearly interpolates, and holds the terminal value to
scan end. Whole blood is unavailable (vb=0), and the blood metadata do not
independently establish the plasma decay convention. Threshold ROIs are not
anatomical brain segmentations. These limitations preclude absolute accuracy
claims. FDG is not used here to validate reversible Logan analysis.

## Reversible PBR28: 20 scans from 10 people

Source: the real [kinfitr PBR28 example dataset](https://github.com/mathesong/kinfitr),
with source processed metabolite-corrected arterial input and whole blood.
The [author's tutorial](https://granvillematheson.com/posts/2020/kinfitr-part-2-basics-iteration)
describes the input fields and modeling workflow. This dataset is distinct
from kinfitr's simulated `simref` dataset. Frontal cortex (FC) was selected
before inspecting these fits; scans are paired test/retest observations and
are not 20 independent people. Images/MRI were not supplied by this benchmark,
so it exercises the kinetic methods, not image registration.

All 20 scans used the existing reversible two-tissue fit, fixed vb=0.05,
source positive frame weights, and no input-delay estimation. Zero-duration
source rows were discarded and zero-weight frames were excluded from fitting
and scoring. Five interleaved folds generated out-of-fold predictions. Logan
was evaluated at t*=30,40,50 min with the same blood correction.

The source PET starts 4–49 seconds after time zero. For Logan integration only,
we explicitly inserted an unmeasured initial interval with zero activity.
This is an assumption, not an observation. Any input endpoint shorter than the
last frame was extended with its last value. Full choices and interval lengths
are recorded per scan in `pbr28_results.json`. Timing, vascular fraction and
initial-interval assumptions can affect both methods.

'''
text+=f"- Median five-fold prediction NRMSE: **{pbr['median_cross_validated_nrmse']:.2%}**; worst scan: **{pbr['max_cross_validated_nrmse']:.2%}**.\n"
text+=f"- Full fits reporting optimizer convergence: **{pbr['converged_full_fits']}/20**, with **1/20 at a rate bound**. Four of the 100 fold fits also hit bounds.\n"
text+=f"- Median absolute difference between Logan (t*=40) and compartment VT: **{pbr['median_absolute_logan_vs_compartment_difference']:.2%}**, including all scans. Agreement is not accuracy: methods share input data and assumptions.\n"
text+=f"- Median absolute test/retest VT difference, divided by the pair mean: **{pbr['median_absolute_test_retest_difference']:.2%}**, retaining the failed pair. Real repeat differences include biology and measurement error.\n"
text+='''
### A failure that matters

For `rbqc_1`, k4 reached the imposed lower bound (approximately 1e-6 min^-1).
The reversible formula VT=(K1/k2)*(1+k3/k4) consequently gives **VT≈10124**,
while Logan at 40 min gives **VT≈2.62**. Its out-of-fold curve NRMSE is still
about **7.94%**. Thus a visually plausible fit and optimizer convergence can
coexist with an unusable equilibrium estimate. The case is retained in the
JSON, summary statistics and report; it is omitted from the agreement plot's
axis range and explicitly annotated there. Do not use its compartment VT for
biological inference. No rate bounds were changed to hide the failure.

## Conclusion and reproducibility

The code runs on measured human data and generally predicts omitted frames
reasonably well. These results support a research prototype with explicit QC,
not a fully validated quantitative pipeline. Human registration accuracy,
absolute kinetic accuracy, robustness to delay/dispersion and clinically
acceptable error thresholds are not established. The instability case and
late-scan prediction errors are material limitations, not passing results.

Reproduce from the package root:

```sh
python examples/fetch_public.py public-pet
python examples/fetch_extended.py public-pet
python - <<'PYTHON'
import sys
sys.path.insert(0, 'examples')
from validate_public import run
run('public-pet', 'validation/real_extended/reproduction', subjects=['s028','s029','s030'])
PYTHON
python examples/stress_real_data.py validation/real_extended/reproduction validation/real_extended/fdg_stress
python examples/fetch_pbr28.py pbr28
Rscript examples/export_pbr28.R pbr28/pbr28.RData pbr28/export
python examples/validate_pbr28.py pbr28/export validation/real_extended/pbr28
```

R base is only needed to export the source RData into CSV. The kinetic analysis
is Python. Install package development dependencies as described in the README.
Raw PBR28 data remain in the working data directory; only results and fetch/
conversion scripts are distributed. A second additional FDG image download
(s031) was stopped before analysis because of slow transfer; it contributes
no results. Completed FDG scans are exactly s028, s029, s030.

Source checksums and the pinned kinfitr URL/commit are in the manifests.
All downloaded s030 objects were verified against hashes from the pinned
OpenNeuro metadata commit used by the original toolkit. Eleven core tests
still pass; see `test_results.txt`. Saved results, TACs, transforms, reduced
images and overlays accompany this report. No model selection or acceptance
threshold was tuned to make these real-data results pass.
'''
text+=f"\nPBR28 source commit: `{source['commit']}`; SHA256: `{source['sha256']}`.\n"
(root/'REAL_DATA_REPORT.md').write_text(text)
fig,ax=plt.subplots(figsize=(8,4))
x=np.arange(len(fdg));ax.bar(x-.18,[100*r['pooled_cross_validated_nrmse'] for r in fdg],width=.36,label='Five-fold omitted frames')
ax.bar(x+.18,[100*r['late_block_nrmse'] for r in fdg],width=.36,label='Final 5–6 frames predicted from earlier frames')
ax.set(xticks=x,xticklabels=[r['subject'] for r in fdg],ylabel='Normalized prediction RMSE (%)',title='Real FDG: the harder temporal test exposes larger errors')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.1),fontsize=8);fig.tight_layout();fig.savefig(root/'fdg_comparison.png',dpi=150)
