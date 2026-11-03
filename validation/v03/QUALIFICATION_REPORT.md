# Research qualification record

## Overall decision: NOT QUALIFIED for the full declared scope

The requested use is **regional brain-PET research analysis**. Numerical and
synthetic-registration verification were evaluated against a locally frozen
protocol. Untouched external kinetic validation and measured-human PET/MRI
registration truth remain unavailable. Neither a passing simulation nor a
new release label substitutes for those mandatory gates.

Protocol: `PETKIN-RESEARCH-2026-11-03-A`. SHA256: `33335f0256ba2e1cba785b614feb3a1d98ffeda136960d577158a3349a5bc7a6`.
The criteria were frozen before the new holdout outcomes were generated.
This is a local prospective freeze, not a public preregistration or independent
certification. Implementation/generator hashes and a documented setup-only
registration-generator correction are included.

| Gate | Result | Requirement |
|---|---|---|
| SW-01 | PASSED | Automated regression and input contracts |
| NUM-01 | PASSED | Routine simulation released macro accuracy |
| NUM-02 | PASSED | Large-error risk on routine simulations |
| ROB-01 | PASSED | Fresh short/noisy/slow-dissociation challenge |
| REG-01 | PASSED | Known-transform MRI-derived synthetic PET |
| EXT-01 | BLOCKED | Untouched external kinetic cohort |
| REG-02 | BLOCKED | Independent measured human PET/MRI registration truth |

## Fresh kinetic tests

The new seed was 61820427. All 130 cases and
outcomes are retained in `holdout_data.json` and `holdout_results.json`.
Ninety routine cases span 1T, reversible 2T and irreversible 2T with known
processed inputs, timing, vascular fraction and Gaussian SDs. Forty hard
reversible cases vary scan duration, slow dissociation and count sensitivity,
with C11-like decay and Poisson/background-subtracted observations. Ground
truth comes from independent adaptive ODE integration of frame activity.
These are numerical/model tests, not independent human biology.

| Set | Returned estimates | Retention | Median absolute macro error among returned | 95th percentile error | Returned errors >50% |
|---|---:|---:|---:|---:|---:|
| Routine | 89/90 | 98.9% | 0.44% | 1.62% | 0 |
| Hard challenge | 8/40 | 20.0% | 4.19% | 24.97% | 0 |

The stricter screen deliberately withholds uncertain results. Rejection is not
successful recovery of the missing kinetic information. Retention floors in
the frozen criteria prevent passing merely by rejecting everything.

The one-sided 95% binomial upper summary for >50% error among returned routine
estimates is **3.31%**. In the hard
set it is **31.23%** because only
8 estimates were returned. Zero observed failures in that
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

- True k4=0.001: returned=False; reasons: macro_profile_open_or_incomplete.
- True k4=0.01: returned=False; reasons: macro_profile_open_or_incomplete.

A closed profile remains conditional on the selected model, input, noise and
preprocessing. It is not proof of biological accuracy or complete uncertainty
propagation. The public low-level fitting API remains a diagnostic numerical
fit; use the stricter research entry point for the defined screen.

## Synthetic registration verification

Twelve fresh known rigid poses on three MRI-derived synthetic PET volumes yielded mean landmark error **0.203 mm**, maximum **0.378 mm**. These satisfy the predefined numerical registration gate. The synthetic modality was produced by an intensity transform, blur and Poisson noise; it is deliberately not called measured human PET ground truth.

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
