# PetKin: PET Preprocessing and Tracer Kinetic Modeling

[![Software verification](https://github.com/kkartikeya-kulkarni/petkin/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/kkartikeya-kulkarni/petkin/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://github.com/kkartikeya-kulkarni/petkin/blob/main/pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**PetKin** is a Python toolkit for preprocessing dynamic positron emission
tomography (PET) images and estimating tracer kinetics from regional
time–activity curves.

[Download v0.4.0.dev4](https://github.com/kkartikeya-kulkarni/petkin/releases/tag/v0.4.0-dev.4)
· [Usage guide](docs/USAGE.md)
· [Methods](docs/NATIVE_PREPROCESSING.md)
· [Report an issue](https://github.com/kkartikeya-kulkarni/petkin/issues)

## About

PetKin connects image preparation, regional curve extraction, kinetic modeling,
and result inspection in one research workflow. You can start with reconstructed
dynamic PET and a structural MRI, or analyze time–activity curves prepared by
another processing pipeline.

The native preprocessing workflow uses NumPy, SciPy, and SimpleITK. It provides
MRI bias correction, brain masking, candidate tissue segmentation, rigid motion
correction, PET–MRI registration, anatomical template mapping, and optional
regional partial-volume correction. The resulting curves can be passed directly
to PetKin's compartmental and graphical models.

PetKin supports:

- MRI bias correction with N4 and brain extraction using a supplied mask or
  template-mask propagation.
- Candidate CSF, gray-matter, and white-matter intensity classes and probability maps.
- Rigid interframe motion correction and alignment of PET to structural MRI.
- Affine and B-spline anatomical mappings, atlas-label propagation, and
  template-space PET output.
- Regional time–activity curve extraction and optional geometric transfer matrix
  (GTM) partial-volume correction with a supplied scanner point-spread function.
- One-tissue, reversible two-tissue, and irreversible two-tissue compartmental fits.
- Logan, MA1, Patlak, SRTM, and reference-tissue Logan analysis.
- Profile-likelihood and parametric-bootstrap uncertainty estimates.
- PET-BIDS blood-input preparation, PETPrep TAC import, and version-pinned
  OpenNeuro retrieval.
- HTML reports, motion and fit plots, residuals, quality flags, and source provenance.

**Current status:** research development release. See the
[validation record](VALIDATION_STATUS.md) for tested scope and outstanding
independent checks.

## Principles

1. **Reproducibility.** Record source hashes, processing settings, model choices,
   and random seeds so an analysis can be traced and repeated.
2. **Inspection.** Provide image overlays, regional coverage, motion summaries,
   fitted curves, and residuals alongside numerical results.
3. **Explicit assumptions.** Keep frame timing, concentration units, blood
   corrections, anatomical inputs, and scanner resolution visible to the analyst.
   Preserve failed or flagged outcomes for review.

## Installation

PetKin requires **Python 3.10 or later**.

```sh
git clone https://github.com/kkartikeya-kulkarni/petkin.git
cd petkin
python -m venv .venv
source .venv/bin/activate
python -m pip install .
petkin --help
```

On Windows, activate the environment with `.venv\Scripts\activate`.

The native workflow runs through the Python dependencies installed above. Supply
an MRI-grid brain mask or an appropriate anatomical template and its brain mask.
Templates and atlas assets are supplied separately under their own licenses.

## Getting started

### Try kinetic modeling with simulated data

Generate a dataset with known parameters, fit a two-tissue model, and run Logan
analysis:

```sh
python examples/make_demo_data.py demo
petkin fit demo/tac.csv demo/input.csv demo/fit.json \
  --model 2t --vb 0.05 --report demo/fit.html
petkin logan demo/tac.csv demo/input.csv demo/logan.json \
  --tstar 30 --vb 0.05 --report demo/logan.html
```

Open the HTML reports in a browser. The generating parameters are saved in
`demo/truth.json`. Modeling CSV times are in **minutes**.

### Preprocess PET and MRI

Supply reconstructed 4D PET, 3D T1-weighted MRI, the original PET-BIDS JSON,
a binary MRI-grid brain mask, and MRI-grid integer region labels:

```sh
petkin preprocess-native pet.nii.gz T1w.nii.gz pet.json preprocessed \
  --brain-mask brain_mask_in_T1w.nii.gz --labels regions_in_T1w.nii.gz
```

Open `preprocessed/index.html` to review alignment, candidate tissue classes,
regional coverage, and motion. Regional curves are written to
`preprocessed/regions/regions.json`. PET-BIDS frame timing remains in **seconds**
in the source JSON and is converted to minutes for modeling.

Template-based extraction, standard-space output, and optional GTM correction are
covered in the [native preprocessing guide](docs/NATIVE_PREPROCESSING.md).

### Analyze regional curves

Once a parent-plasma input has been prepared with matching units, time origin,
and decay reference:

```sh
petkin analyze preprocessed/regions/regions.json prepared-input/input.csv results \
  --model 2t --vb 0.05 --tstar 30
```

Open `results/index.html` for regional results, quality flags, fit and residual
plots, held-frame prediction checks, and CSV exports. Choose the model, vascular
fraction, and graphical start time for your tracer and acquisition protocol.
Use a fresh output directory for each preprocessing or regional-analysis run.

See the [usage guide](docs/USAGE.md) for blood-input preparation, CSV formats,
PETPrep imports, uncertainty estimates, and Python examples.

## Documentation and development

- [Command examples and input formats](docs/USAGE.md)
- [Native preprocessing methods and assumptions](docs/NATIVE_PREPROCESSING.md)
- [Detailed modeling guide](USER_GUIDE.md)
- [Validation status and evidence](VALIDATION_STATUS.md)
- [Methods and literature](LITERATURE_REVIEW.md)
- [Release history](CHANGELOG.md)

Bug reports, reproducible examples, and contributions are welcome through
[GitHub issues](https://github.com/kkartikeya-kulkarni/petkin/issues) and pull requests.
Include the PetKin version, command, relevant metadata, and error output when
reporting a problem.

To run the software tests:

```sh
python -m pip install ".[dev]"
python -m pytest tests -q
```

## Citation

Please cite the PetKin version used in your analysis, with its repository URL
and release identifier. Software citation metadata are available in
[CITATION.cff](CITATION.cff).

Also cite the methods, anatomical atlases, source datasets, and external software
used in your analysis. Method references and data provenance are linked in the
[methods guide](docs/NATIVE_PREPROCESSING.md) and
[third-party record](THIRD_PARTY.md).

## License

Original PetKin code is released under the [MIT license](LICENSE).
Dependencies, datasets, templates, and atlases retain their respective licenses.

## Acknowledgements

PetKin builds on the numerical and imaging capabilities of NumPy, SciPy,
SimpleITK, and Matplotlib, and on publicly shared PET datasets and methods.
PETPrep's workflows informed the architectural review for native preprocessing;
reviewed sources are recorded in the repository. See
[third-party provenance](THIRD_PARTY.md) for details.
