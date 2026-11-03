# Research qualification workflow

The supplied v0.3 candidate is **not fully qualified** for its full declared
research scope. Numerical/synthetic gates passed; external measured-human and
kinetic evidence are blocked. Read `validation/v03/QUALIFICATION_REPORT.md`.

The original locally frozen protocol and its SHA256 are immutable audit records.
The holdout generator refuses to overwrite its dataset. Do not rerun with
selected seeds until a favorable result appears. A new policy or numerical
algorithm requires a new protocol, a new independent test set and retention of
old failures.

## Replay the frozen experiment

From the package root, copy `PROTOCOL.json` and `PROTOCOL.sha256` into a new
empty output directory, then run:

```sh
pytest -q
python examples/qualification_holdout.py fresh_protocol_run
python examples/qualification_registration.py validation/v02/imaging fresh_protocol_run
```

The included implementation freeze documents the code used for the original
run. The registration generator had one recorded setup correction before any
registration pose was evaluated; no criterion, optimizer or seed was changed.
Tests should also be run against an installed wheel. The candidate was tested
in the recorded macOS/Python environment; other platforms are not implicitly
qualified.

The original two missed errors are exposed regression cases:

```sh
python examples/qualification_known_failures.py \
  validation/v02/identifiability_challenge.json fresh_protocol_run
```

They cannot be counted as independent confirmatory evidence. All 130 new cases,
returned/withheld states and truth errors are saved. The larger conditional
risk interval in the small retained stress subset must remain visible.

## Qualification decision

`qualification_status` requires every mandatory gate to pass. Missing,
blocked, failed and not-run evidence prevents qualification. The included
`qualification_report.py` records the actual current external-data blockers;
it does not generate anatomical truth or independently approve a study.
External evidence and its criteria must be established and reviewed before
claiming the full intended-use qualification. See
`validation/v03/EXTERNAL_VALIDATION_CHECKLIST.md`.
