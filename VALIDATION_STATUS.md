# Validation status

**Overall: not qualified for the full declared research scope.**
The original seven-gate protocol and its thresholds have not been relaxed.
Historical v03 records remain unchanged. The newer external addendum records
additional evidence without retrospectively passing its unresolved gates.

## Completed technical evidence

Software tests, numerical parameter recovery, fresh stress simulations,
synthetic geometry checks, an independent 20-scan PS13 evaluation, independent
kinfitr graphical comparisons, and 20 measured PET/MRI automatic alignments.
See the machine-readable records in `validation/v03/` and
`validation/external_ps13/`.

## PS13 extraction correction and anatomical follow-up

The historical intensity-threshold PS13 experiment is not anatomical brain
validation. In the first audited baseline scan (PSBB01 / S3 alias PS50), the old
mask has zero overlap with the FSL-estimated MRI brain mask after independent
FSL rigid alignment. Its original numerical records remain available, but those
results cannot establish regional brain kinetic performance. The reproduction
script now requires an explicit historical-experiment flag.

A separate exploratory workflow used MRI-defined Harvard–Oxford atlas regions
with FSL brain extraction, MRI normalization and PET alignment. All 96 regional
analyses completed: 2/48 baseline and 45/48 blocked regions passed local checks
with source timing unchanged. Timing sensitivity suggests a scan-specific
baseline alignment question; no fitted or hypothetical clock offset has been
adopted as verified timing. No motion correction was added in this exploratory
comparison, and these atlas-derived regions are not expert-labelled truth.
See `validation/ps13_anatomical_fsl/` for all reported hypotheses and scope.
PETPrep is an optional comparison pipeline, not a requirement to use PetKin.

## Outstanding EXT-01: external kinetic qualification

The v0.4.0.dev4 native preprocessing implementation adds N4 correction, candidate
tissue segmentation, template-mask propagation, anatomical normalization, motion
correction and optional regional GTM correction. These additions have their own
component and integration tests in `tests/test_preprocessing.py`; they do not
inherit full anatomical qualification from the older rigid-registration tests.
PETPrep was reviewed as a workflow reference, not copied or executed as a backend.
See `validation/native_preprocessing_v04/` for implementation and execution scope.

A reviewer independent of development must verify the source records supporting
PET–blood calibration, decay reference, input/tissue timing, parent correction,
vascular assumptions, and anatomical ROI geometry. Existing raw estimates and
prediction metrics do not supply these facts. The PS13 experiment used an
exploratory uptake mask and omitted motion correction. None of the four repeat
pairs passed local QC on both scans; all 20 research estimates were withheld.
Further development based on these outcomes needs a new untouched cohort for
confirmatory testing. These are unresolved requirements, not paperwork-only gaps.

## Outstanding REG-02: measured human accuracy

At least ten measured PET/MRI pairs require independent corresponding landmarks
or marker-derived reference transforms and a reviewer-approved error tolerance
fixed before evaluation. We currently have twenty image pairs and zero independent
references. RIRE's official project page says truth transforms remain sequestered
and its earlier submission service is not maintained:
https://rire.insight-journal.org/

An independent rater must mark native MRI and PET images without viewing this
algorithm's fitted alignment. Use physical LPS millimetres. The blank forms and
instructions are in `validation/external_ps13/registration_review/`.
`scripts/evaluate_registration_review.py` can calculate errors once a signed-off
protocol, independently supplied points, and transform files are available. It
cannot establish reviewer independence or anatomical correctness automatically.

## Claims suitable for a release

"An auditable PET kinetic-modeling research candidate with numerical verification,
independent software comparisons, external prediction experiments and explicit
quality-control withholding. Full external kinetic and measured-human registration
validation are pending."

No full-validation, state-of-the-art, clinical, or universal-tracer claim is
supported by the current evidence. Independent review must evaluate the failures,
not merely endorse successful plots.
