# Scientific basis and comparison scope

Focused review checked on **2026-11-03**. This is a targeted methods/software
review, not an exhaustive systematic review of every PET method available.
A defensible “among the best” claim would need independent comparative studies
across additional tracers, sites and image-registration ground truth. The
present release does not make that ranking claim.

| Study or primary source | What it establishes | Implementation decision |
|---|---|---|
| [Logan et al., 1990](https://pubmed.ncbi.nlm.nih.gov/2384545/) | Graphical estimation for reversible tracer binding using plasma and tissue measurements | Retain plasma Logan with explicit t*, blood correction and integration conventions |
| [Ichise et al., 2002](https://pubmed.ncbi.nlm.nih.gov/12368666/) | Alternative linear formulations evaluated for noise effects | Add MA1 and report agreement/sensitivity, without claiming universal superiority |
| [Patlak et al., 1983](https://pubmed.ncbi.nlm.nih.gov/6822610/) | Linear graphical inference for irreversible transfer after appropriate equilibration | Add Patlak for the FDG/irreversible branch; do not substitute reversible VT for Ki |
| [Lammertsma & Hume, 1996](https://pubmed.ncbi.nlm.nih.gov/9345505/) | Simplified reference-tissue receptor modeling | Add SRTM with explicit frame means and reference-region assumptions |
| [Logan et al., 1996](https://journals.sagepub.com/doi/10.1097/00004647-199609000-00008) | Reference-region graphical DVR estimation | Include the reference Cr/k2′ term; require k2′ and a valid reference region |
| [Raue et al., 2009](https://doi.org/10.1093/bioinformatics/btp358) | Profile likelihood can expose practical identifiability and open confidence intervals | Add nuisance-refitted profile grids; never report a grid boundary as a closed interval |
| [Reducing modeling error of graphical methods, 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2900488/) | Finite-duration model error can affect Logan and MA1 in addition to noise | Report t* sensitivity and retain failure cases; a high R² is insufficient |
| [Tjerkaski et al., 2020, kinfitr validation](https://link.springer.com/article/10.1186/s13550-020-00664-8) | Agreement and repeatability were evaluated against established software across several tracers | Execute independent point-estimate comparisons and distinguish agreement from biological truth |
| [NiftyPAD, 2023](https://link.springer.com/article/10.1007/s12021-022-09616-0) | Reference-tissue methods were validated against other software; includes SRTM variants | Compare against installed NiftyPAD 2.0.0 SRTM basis on independent ODE cases |
| [PETPrep workflow documentation](https://petprep.readthedocs.io/en/latest/workflows.html) | Comprehensive preprocessing uses dedicated registration tools and transform concatenation | Compose transforms for final output, retain QC and allow external-transform resampling |
| [Reversible 2TCM identifiability preprint, 2025](https://arxiv.org/abs/2504.09959) | Mathematical identifiability depends on available measurements and assumptions | Treat local numerical screening as limited; do not equate optimizer convergence with identifiability |

## Competitor comparisons and limits

**kinfitr:** The retrieved DESCRIPTION reported version 0.9.7, dated 2026-11-01.
The audit executes the unmodified upstream Logan, MA1 and 2TCM modeling
functions, sourced at the commit pinned in
`validation/v02/kinfitr_source_manifest.json`. Required R dependencies were
installed in the working directory. This is a source-function comparison,
not an installation or certification of the full kinfitr application. The
upstream standard-error helper returns unavailable values without its optional
`car` dependency; this audit compares point estimates only. Graphical estimates
use matched t*=40 min, vb=.05, unweighted regression and frame-duration tissue
integration. Small differences remain from input interpolation and blood
sampling conventions. The three nonlinear comparisons are illustrative:
kinfitr's midpoint observations differ from petkin's frame-averaged predictions.
They are not used as an exact-equality acceptance test.

**NiftyPAD:** Version 2.0.0 was installed and its SRTM basis and reference Logan
functions executed. Twelve independent ODE cases test the reference-tissue
branch. NiftyPAD uses a second-grid convolution and a finite basis grid; these
numerical conventions differ from petkin's frame-averaged nonlinear fit. The
NiftyPAD nonlinear `srtm` entry point failed with the installed SciPy/NumPy
combination because its heterogeneous input tuple was coerced to a numeric
array. The released basis implementation worked and was used without modifying
its modeling code. This compatibility issue is not a claim that petkin is
scientifically superior to NiftyPAD.

**PETPrep:** Current workflow/usage documentation was reviewed, including
reference selection, robust registration, transforms, masks and QC. PETPrep,
ANTs, FreeSurfer and FSL were not installed or benchmarked here. Our SimpleITK
rigid pipeline is much smaller in scope. Supporting externally estimated
transforms and a single final resampling does not reproduce PETPrep's full
segmentation, partial-volume correction or anatomical validation.

**PMOD and other commercial tools:** Their role in published validation was
reviewed. No licensed executable was available or run. Results from another
study are not treated as our own head-to-head benchmark.

## What this release improves

The release adds reproducible log-rate multistarts, optional vascular/delay
fitting, clear relative-versus-absolute noise weights, local covariance and
macroparameter uncertainty, explicit failure reasons and withheld unsafe VT/Ki,
conditional bootstrap, and profile-likelihood diagnostics. MA1, Patlak, SRTM
and reference Logan extend the original model set. PET-BIDS parsing rejects
unknown units/correction assumptions, and final PET-to-MRI output is resampled
from original frames using composed transforms.

The inference policies are conservative engineering rules, not clinically
validated acceptance thresholds. The difficult simulation challenge still
contains large errors that pass the local screen. Profiles and bootstrap are
conditional on the selected model, input and preprocessing and cannot resolve
all missing information or model misspecification.

## Scope still missing for a leading-toolkit ranking

Independent multicenter/tracer cohorts not used during development; human
landmark/segmentation registration accuracy; validated partial-volume and
attenuation handling; propagation of blood/reference/image uncertainty;
voxelwise throughput and memory benchmarks; broader external-software and
clinical reproducibility studies; and a release/install matrix across Python,
operating systems and dependency versions. These remain evidence gaps, not
features or validations implied by the current results.
