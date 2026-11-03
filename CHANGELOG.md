# 0.4.0.dev1 — regional workflow

- Added `petkin extract` for supplied anatomical label images with physical-grid checks.
- Added `petkin analyze` for regional manifests, source weights, per-region error isolation,
  CSV exports, held-frame prediction, fit/residual plots and linked HTML reports.
- Preserved failed and flagged results; numerical models and thresholds are unchanged.
- Completed 120 PBR28 regional development fits across 20 scans: 113 passed local QC;
  median held-frame NRMSE 4.39%. These are exposed development data, not a new
  independent qualification. All scientific limitations remain explicit.
- Passed 48 tests, including independent-ODE image-to-region-to-fit recovery.
- Updated CI actions and avoided capacity-constrained hosted macOS ARM queues.

# GitHub packaging addendum — 2026-11-03

Numerical core unchanged from 0.3.0rc1 and documented pre-existing amendments.
Added MIT packaging metadata, GitHub installed-wheel CI, publication instructions,
current evidence summaries, release-integrity checks and an independent-landmark
evaluation tool. No research gate was relaxed or retrospectively marked passed.
New scorer tests use synthetic fixtures and are not measured-human validation.

# Changes

## 0.3.0rc1 — 2026-11-03

- Added direct VT/Ki profile likelihood with exact macro constraints and refitted
  confidence-crossing roots; rate profiles use conservative grid brackets.
- Added a stricter conditional research screen, explicit input assumptions,
  residual/precision checks and fail-closed qualification gates.
- Froze a local protocol before 130 fresh simulation outcomes and 12 synthetic
  registration poses; retained all rejected and failed results.
- Full declared-scope research qualification remains blocked by missing
  independent kinetic and measured-human registration evidence. This candidate
  must not be described as a fully validated toolbox.

## 0.2.0 — 2026-11-03

- Added deterministic log-rate multistart fitting with optional blood-volume
  and input-delay estimation, explicit extrapolation and noise-scale conventions.
- Withhold VT/Ki on blocking local QC; retain raw values, bound names, covariance,
  residual diagnostics, information criteria and settings for audit.
- Added conditional nuisance-refitted profile likelihood and parametric bootstrap.
- Added MA1, Patlak, SRTM and reference Logan, and exact linear-input integrals.
- Added strict PET-BIDS helpers, portable QC HTML and hashed input provenance.
- Compose transforms for one final resampling from original input frames; added
  registration masks/multiple initializations, displacement and landmark QC.
- Added independent kinfitr/NiftyPAD comparisons, 72 routine simulations,
  24 hard count-noise cases, real-data reanalysis and a focused literature review.

### Compatibility

`VT`/`Ki` can now be null when checks fail. Audit-only numerical values are in
`raw_VT`/`raw_Ki`. Supplied sigma is relative unless `absolute_sigma=True`.
Existing `motion_correct` and `pet_to_mri` retain their stage-specific behavior;
the CLI now uses `preprocess_dynamic` for one final interpolation.

### Known limitations

The hard challenge includes two >50%-error cases that pass local checks.
External clinical validation, independent human registration ground truth,
complete uncertainty propagation and a best-toolkit ranking are not established.
