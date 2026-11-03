# Evidence needed to complete research qualification

The numerical core is verified only within its frozen simulation/test scope.
These open items cannot be satisfied by relabeling development cases, changing
thresholds after results, or treating an existing registration as truth.

## EXT-01: external kinetic cohort

Obtain an independent research cohort not inspected during development, with:

- Deidentified participant, session, site and tracer identifiers.
- Dynamic regional tissue curves and exact frame boundaries.
- Calibrated processed plasma, appropriate parent/metabolite processing,
  blood where needed, documented clocks and decay-reference times.
- ROI provenance, scanner/reconstruction details and processing decisions.
- A reference analysis and independently reviewed choices of model, weights,
  bounds, vascular fraction, delay and graphical t*.
- Frozen intended-use and numerical-agreement/repeatability criteria before
  evaluating the cohort; uncertainty, withheld estimates and failures retained.

Published biological truth is generally unavailable; distinguish implementation
agreement/repeatability from biological accuracy. Do not use the 20 exposed
PBR28 scans or three exposed FDG scans as an untouched external cohort.

The independent multi-tracer kinfitr/PMOD validation paper's participant data
are institutionally restricted. No application or communication was made on
the user's behalf. A future institutional collaborator can request data under
its own approvals and sharing conditions.

## REG-02: measured-human registration accuracy

Obtain at least ten measured PET/MRI pairs with independent landmark
correspondences or reference transforms, and freeze an expert-reviewed error
tolerance before evaluation. Record patient-level independence and contrasts;
multiple contrast pairs from one person do not establish broad population
coverage. Document image geometry, anatomical coverage, coordinates and units.

Use `landmark_registration_error` with independent physical LPS points and the
fixed-to-moving convention. Examine failure rate and error distribution, not
just the mean. MRI-derived synthetic PET, visual overlays and transforms fitted
by the same method are not independent measured-human truth.

The RIRE training PET/MRI/reference archives returned HTTP 403 from the official
public links and alternate public gateways in this environment. Other RIRE
reference transforms are sequestered. Record this as unavailable evidence;
do not manufacture correspondence labels or assume identity alignment.

## Release decision

The frozen protocol, code/data hashes, source/wheel tests, full holdout outcomes
and external review records must be consistent. Mandatory failed, blocked or
not-run gates prevent a full declared-scope qualification. Changing the
acceptance criteria or scientific algorithm requires a new protocol and fresh
confirmatory evidence; old failures remain visible.
