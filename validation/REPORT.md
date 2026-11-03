> Follow-up: see [extended testing on 23 real scans](real_extended/REAL_DATA_REPORT.md), including a failed equilibrium estimate and stronger prediction tests.

# Validation report

## Scope

This is a research implementation with numerical tests, an independently
simulated kinetic cohort, known-motion phantom checks and two public FDG
subjects reserved for out-of-sample evaluation. Numerical recovery in simulation
and real-data prediction are distinct forms of evidence. Neither establishes
clinical validity or absolute accuracy of real-world kinetic parameters.

## Automated tests

See `test_results.txt`: 11 tests passed. Coverage includes independent ODE
agreement for three compartment models, recovery of known parameters, analytic
one-tissue Logan VT, blood correction, step-size convergence, invalid inputs,
roundoff in real BIDS frame boundaries, blank-frame handling, rigid alignment
and ROI geometry.
The command-line fit and Logan workflows were also exercised successfully.

## Synthetic evaluation cohort

Twelve new synthetic parameter sets with fixed seed 20261102 were evaluated
separately from the deterministic unit-test cases. Ground truth came from an
adaptive ODE solver with an integrated observation state, independent of the
convolution implementation used to fit the data. Frames span 120 minutes;
independent Gaussian noise has SD 1% of activity with a 5%-peak activity floor.
This is not a scanner count/noise or reconstruction model.

| Parameter | Median absolute relative error | Maximum absolute relative error |
|---|---:|---:|
| K1 | 0.17% | 0.94% |
| k2 | 0.81% | 3.08% |
| k3 | 1.60% | 5.08% |
| k4 | 0.76% | 1.90% |

Median VT error: compartment fit **0.10%**; plasma Logan **0.09%**. All fits converged.

Simulation acceptance checks passed: all fits converge, median parameter errors
below 10%, and median VT errors below 5% for both methods. Thresholds are
engineering checks, not preregistered clinical endpoints. Logan t* = 50 min.
The plots, every TAC and its input curve, known parameters and full fit
residuals are saved alongside this report.

## Registration checks

An asymmetric 3D phantom was transformed with known rotations of
(0.025, -0.035, 0.04) rad and translation (3, -2, 2) mm. Three physical-space
landmarks were evaluated. A square-root intensity transform supplied a
multimodal contrast surrogate, not a biological MRI simulation.

- PET_motion: mean landmark error **0.237 mm**, maximum **0.423 mm**.
- PET_MRI_contrast_surrogate: mean landmark error **0.349 mm**, maximum **0.635 mm**.

All landmark errors were below the 1-mm acceptance threshold. The dynamic
PET-to-MRI workflow also produced the expected MRI spatial grid. These checks
do not measure registration accuracy on human images.

## Public held-out FDG evaluation

[OpenNeuro ds004513](https://openneuro.org/datasets/ds004513/versions/1.0.6),
*The energetic costs of the human connectome*, CC0. Source metadata commit:
`37590f49c0bf6cb11ea75ebdc36266277c87fe98`. Subjects s028 and s029 are from the
prospective replication cohort and were not used for simulation or optimizer
selection. Every fifth frame (indices 4, 9, 14, 19, 24, 29) was excluded from
fitting. Those six measured frames per subject were compared with predictions
from the remaining eligible nonblank frames. Blank frames are flagged
and preserved with identity transforms because motion is unmeasurable, then
excluded from fitting and prediction scoring. NRMSE is RMSE divided by the mean held-out activity.

| Subject | Held-out NRMSE | Ki (mL/mL/min) | Fit converged | At bound | Excluded blank frames |
|---|---:|---:|---|---|---|
| s028 | 6.61% | 0.04195 | True | False | [] |
| s029 | 6.20% | 0.04075 | True | False | [0] |

Processing includes 4-mm resampling, framewise rigid motion correction to the
last PET frame, duration-weighted mean PET registration to MRI, and a PET
uptake threshold ROI (35% of maximum mean activity). The irreversible two-tissue
model was fit with uniform frame weights and fixed vb=0. All transforms,
optimizer diagnostics, TACs, reduced images and overlays are provided.
The saved overlays were visually inspected for gross brain alignment; no
landmark or segmentation ground truth was available for human images.

These are exploratory prediction checks; no real-data accuracy pass threshold
was asserted. The curves show residual model mismatch. Input plasma is a
previously modeled curve supplied only at frame midpoints. We assume zero at
injection, linear interpolation, and a terminal constant extension. Whole blood
is unavailable. The available blood metadata did not independently establish
its decay-correction convention. These limitations, the threshold ROI and
resampling prevent interpreting the fitted rates as validated absolute values.
Human FDG data do not validate reversible Logan VT; Logan was validated here
against reversible simulated data only. Larger independent cohorts with
densely sampled calibrated plasma, MRI segmentation and expert registration
review remain necessary for real-world accuracy claims.

## Reproduction

Follow the commands in the package README. `public_manifest.json` records URLs,
file sizes and SHA256 digests of the downloaded source objects. Raw data are
kept outside the distribution; reduced demonstration images and source
metadata are included. `environment.txt` records installed dependency versions.
