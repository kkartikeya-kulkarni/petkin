# Independent review handoff

No independent reviewer has been appointed or has approved this release.
This document records required work; it is not a certificate.

## Kinetic reviewer

Review all 20 PS13 outcomes, including every withheld estimate and the four
repeat-scan pairs. Verify PET–blood cross-calibration, common decay reference,
injection/input alignment, frame timing, metabolite processing, vascular model
and anatomical ROI placement from source records. Record the supporting source
for every assertion. The current exploratory uptake region and kinetic-only
processing do not complete regional anatomical or motion validation. Assess the
large test–retest discrepancies and specify development changes separately from
confirmatory evaluation. If development uses these outcomes, require an untouched
cohort for the new confirmation; these scans are already exposed.

## Registration reviewer / rater

Use native measured PET and MRI images from the full evidence archive. Do not
view automatic transforms, resampled PET, or overlays before placing reference
points. Confirm that landmarks are identifiable independently in both modalities,
and document localization uncertainty. Fill the landmark template with physical
LPS millimetres and at least three non-collinear corresponding points per scan.
The frozen gate requires at least ten measured pairs. Report participant clustering.

Specify acceptable mean and maximum errors in the protocol template, justify them
for the research task, name the reviewer, state independence, and lock the protocol
before evaluating errors. The blank templates intentionally fail evaluation.
Update CASES_TEMPLATE.json to bind it to the SHA256 of the approved protocol.

Run from the repository root, using completed evidence files:

```sh
python scripts/evaluate_registration_review.py reviewed-protocol.json independent-landmarks.csv reviewed-cases.json new-review-results.json
```

The case transform paths are relative to the cases JSON directory. Preserve the
recorded transform hashes. The evaluator reports conditional geometric scores;
it cannot authenticate raters, determine anatomy, or certify independent review.
A separate signed reviewer report must consider these scores and all assumptions.

## Author response

Record every reviewer finding and response. Do not edit historical failed results
or convert unavailable evidence into a passed gate. Publish reviewer names and
reports only with their permission. A user-approved source clarification and review-resource request was sent to the
PS13 corresponding author on 2026-11-03; no response or approval is recorded.
