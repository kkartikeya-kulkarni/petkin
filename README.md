# petkin

A Python toolkit for dynamic PET tracer kinetic modeling, graphical analysis,
rigid motion correction, and PET-to-structural-MRI registration.

**Version:** 0.4.0.dev1 · Regional research workflow

## Installation

Requires Python 3.10 or later.

```sh
git clone https://github.com/kkartikeya-kulkarni/petkin.git
cd petkin
python -m venv .venv
source .venv/bin/activate
python -m pip install ".[dev]"
petkin --help
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Quick start

Generate a reproducible simulated dataset, fit a reversible two-tissue model,
and write a JSON result and HTML report:

```sh
python examples/make_demo_data.py demo
petkin fit demo/tac.csv demo/input.csv demo/fit.json \
  --model 2t --vb 0.05 --report demo/fit.html
petkin logan demo/tac.csv demo/input.csv demo/logan.json \
  --tstar 30 --vb 0.05 --report demo/logan.html
```

Open `demo/fit.html` or `demo/logan.html` in a browser. The generated
`demo/truth.json` contains the simulation parameters.

## Analyze a complete set of regions

Run all supplied regions in one command:

```sh
petkin analyze demo/regions.json demo/input.csv demo/results \
  --model 2t --vb 0.05 --tstar 30
```

Open `demo/results/index.html`. It links to an exportable `results.csv`, fit
plots, residuals, held-frame prediction checks and individual JSON/HTML reports.
Failed and flagged regions remain in the table; they are never silently dropped.
Use a new output directory for every run.

For your own extracted curves, create a region manifest with paths relative to
that manifest:

```json
{
  "regions": [
    {"id": "FC", "name": "Frontal cortex", "path": "frontal.csv"},
    {"id": "THA", "name": "Thalamus", "path": "thalamus.csv"}
  ]
}
```

Each regional CSV uses the TAC format below. Instead of `sigma`, a `weight`
column may supply nonnegative frame weights. Zero-weight frames are excluded
from compartment fitting and retained in the plotted observations. Do not supply
both `weight` and `sigma`.

## From images to regional results

Supply an anatomical label image from your segmentation workflow. After rigid
preprocessing, the labels and PET must share the same spatial grid:

```sh
petkin preprocess dynamic_pet.nii.gz T1w.nii.gz pet.json processed
petkin extract processed/pet_in_mri.nii.gz anatomical_labels.nii.gz pet.json regions
petkin analyze regions/regions.json input.csv regional-results \
  --model 2t --vb 0.05 --tstar 30
```

`extract` uses the supplied positive integer label IDs, checks physical geometry,
and reports region volumes. It performs no PET-intensity threshold segmentation
or implicit interpolation of label images. You can replace the generated
`region_ID` names in `regions.json` with anatomical names from your label atlas.
Choose the model, vascular fraction, input processing and graphical start time
for your study. Graphical methods that cannot operate on a region's timing are
reported as unavailable while compartmental results are retained.

## Prepare your data

The modeling commands use two CSV files. **CSV times are in minutes.** Tissue,
plasma and whole-blood activity must use matching concentration units and a
common time origin and decay reference. `plasma` is the processed parent-tracer
plasma input.

**Tissue activity — `tac.csv`:**

```csv
start,end,activity,sigma
0,0.5,250,25
0.5,1,1100,55
1,1.5,1800,70
```

`start` and `end` define each frame; `activity` is its mean concentration. The
optional `sigma` column provides positive standard deviations or relative uncertainty scales (weights are
proportional to `1 / sigma²`).
Use a complete scan with enough frames for the selected model; the rows above
only illustrate the format.

**Arterial input — `input.csv`:**

```csv
time,plasma,blood
0,0,0
0.1,400,320
0.2,1200,960
```

Provide increasing sample times beginning at zero and covering the scan. Include
`blood` when using a nonzero vascular fraction or fitting it. Input samples are
interpolated linearly. The rows above are a format example, not a complete input.

## Compartmental fitting

| Model flag | Model | Main outcome |
|---|---|---|
| `1t` | One tissue compartment | `VT` |
| `2t` | Reversible two tissue compartments | `VT` |
| `2t_irreversible` | Irreversible two tissue compartments | `Ki` |

```sh
petkin fit tac.csv input.csv fit.json --model 1t --vb 0.05
petkin fit tac.csv input.csv fit.json --model 2t --vb 0.05 --fit-delay
petkin fit tac.csv input.csv fit.json --model 2t_irreversible --vb 0.05
```

Useful options:

- `--starts 8` and `--seed 0`: reproducible multistart fitting.
- `--fit-vb`: estimate the vascular fraction; requires whole blood.
- `--fit-delay`: estimate the input delay within the configured bounds.
- `--dt 0.025`: prediction integration step, in minutes.
- `--absolute-sigma`: treat `sigma` as known measurement standard deviations.
- `--report result.html`: add a portable HTML report.

For the additional research screening workflow, copy the assumptions template
and complete it using your study records:

```sh
cp examples/RESEARCH_ASSUMPTIONS_TEMPLATE.json assumptions.json
petkin research-fit tac.csv input.csv research-fit.json \
  --model 2t --vb 0.05 --assumptions assumptions.json \
  --report research-fit.html
```

## Graphical analysis

Select `tstar` in minutes for your tracer and acquisition protocol.

```sh
petkin logan tac.csv input.csv logan.json --tstar 30 --vb 0.05
petkin ma1 tac.csv input.csv ma1.json --tstar 30 --vb 0.05
petkin patlak tac.csv input.csv patlak.json --tstar 30 --vb 0.05
```

Logan and MA1 estimate `VT`; Patlak estimates `Ki` for effectively irreversible
kinetics. Graphical commands require contiguous tissue frames from time zero.

## Motion correction and PET–MRI registration

Supply a 4D dynamic PET image, a 3D structural MRI, and the PET-BIDS JSON containing
`FrameTimesStart` and `FrameDuration` arrays. **PET-BIDS timing is in seconds**;
the preprocessing command converts it internally.

```sh
petkin preprocess dynamic_pet.nii.gz T1w.nii.gz pet.json processed
```

The output directory contains:

- `pet_in_mri.nii.gz`: motion-corrected PET in MRI space.
- `motion_corrected.nii.gz`: PET aligned to the selected PET reference frame.
- Frame transforms, the MRI-to-PET transform, and composed transforms.
- `qc.json`: processing settings, motion information and provenance.

Use `--reference INDEX` to select a zero-based reference frame, or
`--preserve-blank` to explicitly preserve blank frames with identity motion.
Transforms map fixed physical coordinates to moving physical coordinates in
**LPS millimetres**. Inspect the resulting alignment before extracting regions.

## Uncertainty estimates

Profile a parameter from an existing compartmental fit:

```sh
petkin profile demo/tac.csv demo/input.csv demo/vt-profile.json \
  --fit-result demo/fit.json --parameter VT \
  --minimum 0.5 --maximum 5 --points 41
```

Run a parametric bootstrap:

```sh
petkin bootstrap demo/tac.csv demo/input.csv demo/bootstrap.json \
  --fit-result demo/fit.json --samples 100 --seed 0
```

Choose profile ranges and bootstrap settings appropriate to your data.

## Python API

```python
import numpy as np
from petkin import Frames, InputFunction, fit_compartment, logan

tac = np.genfromtxt("demo/tac.csv", delimiter=",", names=True)
blood = np.genfromtxt("demo/input.csv", delimiter=",", names=True)
frames = Frames(tac["start"], tac["end"])
input_function = InputFunction(blood["time"], blood["plasma"], blood["blood"])

fit = fit_compartment(
    frames, input_function, tac["activity"],
    model="2t", sigma=tac["sigma"], vb=0.05,
)
graphical = logan(frames, tac["activity"], input_function, tstar=30, vb=0.05)
print(fit["parameters"])
print(graphical["VT"])
```

Reference-tissue methods are available through `ReferenceCurve`, `fit_srtm` and
`logan_reference`; see [USER_GUIDE.md](USER_GUIDE.md) for their assumptions and
examples.

## Read the outputs

JSON results include model parameters, predictions, residuals, uncertainty
information, quality flags and file provenance. Check `status`, `qc_reasons` and
`usable_for_quantification` before interpreting an outcome. A withheld `VT` or
`Ki` is stored as `null`; the corresponding `raw_VT` or `raw_Ki` remains available
for inspection. Optimization convergence alone does not establish suitability.

## Tests and documentation

```sh
python -m pytest tests -q
```

- [Full user guide](USER_GUIDE.md)
- [Validation status and evidence](VALIDATION_STATUS.md)
- [Methods and literature](LITERATURE_REVIEW.md)
- [Third-party provenance](THIRD_PARTY.md)

## License

Original code is available under the [MIT license](LICENSE). External datasets
and comparison packages retain their own licenses.
