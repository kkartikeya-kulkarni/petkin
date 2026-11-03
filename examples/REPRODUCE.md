# Reproduce the 0.2 evidence

Run from the petkin package directory. Commands may download several GB of
public images. The supplied validation folder already contains frozen results;
use a new output folder to preserve the original audit.

```sh
pip install -e '.[dev,benchmarks]'
mkdir -p replay
pytest -q
python examples/validate_v02_simulation.py replay
python examples/challenge_v02.py replay
python examples/benchmark_niftypad.py replay
```

## Measured PBR28 data

R base converts the exact source RData into CSV. PET modeling then runs in Python.

```sh
python examples/fetch_pbr28.py pbr28
Rscript examples/export_pbr28.R pbr28/pbr28.RData pbr28/export
python examples/benchmark_v02_real.py pbr28/export replay
python examples/audit_v02_uncertainty.py pbr28/export replay
```

The full fits use eight starts and CV fits four, fixed seeds. No tuning or
selection of t* is performed using the held-out fold. The people/data were
already exposed during development; this is not a new external validation set.

## Source-level independent kinfitr comparisons

The modeling functions are downloaded unchanged at the pinned source commit.
The complete kinfitr application is not installed. The R `car` package is not
needed for these point-estimate comparisons; upstream unavailable SEs are not
used in the audit.

```sh
mkdir -p Rlib kinfitr_reference
Rscript -e 'install.packages(c("pracma","purrr","minpack.lm","tibble"), repos="https://cloud.r-project.org", lib="Rlib")'
python examples/fetch_kinfitr_reference.py kinfitr_reference
Rscript examples/reference_kinfitr.R Rlib kinfitr_reference pbr28/pbr28.RData replay
```

The release compares the saved CSV results against `real_v02.json` at t*=40 min.
Graphical regression is unweighted, vb=.05, with duration-based tissue integrals.
The three nonlinear examples are diagnostic because midpoint versus frame-mean
observations differ. The original source manifest and R session are included
with the release. Optional dependency versions may change; the recorded release
environment is in `validation/v02/environment.txt` and `kinfitr_R_session.txt`.

## FDG images and source metadata

```sh
python examples/fetch_public.py public-pet
python examples/fetch_extended.py public-pet
python examples/validate_v02_imaging.py public-pet replay
```

This processes s028/s029/s030 at >=4 mm. Inputs are downsampled before testing
the single-final-interpolation workflow. It does not establish human landmark
accuracy. For the original frozen FDG regional comparisons:

```sh
python examples/validate_v02_fdg.py validation/real_extended replay
```

Raw-input re-extraction is available through `validate_public.run`, as described
in the earlier real-data report. All frozen PET threshold ROI and input
boundary/decay limitations remain relevant.

## Reports

`write_v02_report.py` consumes the saved release evidence, including the
software-agreement and paired-subject summary files. Run against the provided
complete evidence folder to regenerate the portable report and plots:

```sh
python examples/write_v02_report.py validation/v02
```

No released threshold was tuned to remove the challenge failures. Compare
point estimates, withheld estimates, interval closure, failure rates and
prediction scores together. A good curve alone is insufficient.
