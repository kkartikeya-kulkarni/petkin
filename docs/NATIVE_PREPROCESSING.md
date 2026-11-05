# Native preprocessing

`petkin preprocess-native` is an original PetKin workflow using NumPy, SciPy and
SimpleITK. It does not invoke PETPrep, FSL, ANTs executables, or FreeSurfer. The
earlier `preprocess` command remains available for the smaller rigid-only workflow.

## Processing stages

1. Validate image dimensions, finite observations, timing, physical grids and
   explicitly supplied masks/labels. Inputs are never overwritten.
2. Use an MRI-grid binary brain mask, or propagate a supplied template mask into
   the MRI using affine and optional B-spline registration. This is an estimated
   brain extraction, requiring inspection, rather than an uptake threshold.
3. Estimate an N4 MRI bias field within positive masked voxels and correct the
   full-resolution MRI. Zero voxels remain zero and are excluded from field fitting;
   at least 95% of masked voxels must be positive. Negative masked intensities fail.
4. Fit a three-component Gaussian intensity mixture using expectation maximization.
   Save candidate CSF/GM/WM classes and posterior maps, ordered by T1 intensity.
   These are intensity-based candidates, not reviewed tissue anatomy or parcellations.
5. Estimate rigid frame motion, then register a duration-weighted PET mean to the
   brain-extracted MRI using mutual information. The default reference frame is the
   frame with highest total activity; `--reference` sets an explicit zero-based index.
6. Compose transforms and sample each original frame directly onto the MRI grid.
   Save individual and composed transforms, displacement summaries, and images.
7. Optionally propagate template atlas labels using nearest-neighbour sampling.
   Keep integer IDs and report missing IDs and regional brain-mask overlaps.
8. Extract regional mean curves without changing clocks, decay correction or units.
9. Optionally estimate regional GTM-corrected activities from a user-supplied scanner
   point-spread function. Corrected and uncorrected curves are saved separately.
10. Optionally map each original frame directly to a supplied template grid through
    a separately estimated template-to-MRI pull transform, avoiding a second PET
    image interpolation. Save positive-Jacobian screening statistics and transforms.
11. Write overlays, motion curves, all settings, package versions, source hashes,
    and a processing status. Partial failure cannot produce a completed status.

All transform files map output-grid physical LPS millimetres to source-image LPS
millimetres. A file named `template_to_mri_pull.h5` therefore maps template points
to MRI sampling coordinates. Do not interpret filenames as forward voxel transforms.
No inverse nonlinear map is assumed or silently approximated.

## Requirements and interpretation

Inputs are reconstructed scalar 4D PET, scalar 3D T1w MRI, and the original PET
JSON. A reviewed MRI-grid binary brain mask or a template with its brain mask is
required. MRI-grid regional labels and template-grid atlas labels are alternatives.
The template and atlas are user-supplied; PetKin does not download atlas assets or
grant rights to them. MRI contrasts other than conventional T1w are not supported
by the named tissue-class interpretation.

Deformable registration uses mutual information, an affine initialization and a
smooth B-spline warp, with a nominal 40 mm control-point spacing. Both stages use
multiresolution estimation. The sampled displacement Jacobian must be finite and
positive. This screen does not prove anatomical correctness or global invertibility;
independent human registration accuracy remains unmeasured.

Motion summaries measure displacement of seven physical landmarks, not Power
framewise displacement. Motion within a frame, attenuation-map mismatch, scatter,
reconstruction, and blood calibration are outside this pipeline. Blank frames fail
by default. `--preserve-blank` keeps them and explicitly records unmeasurable motion.
Tracer uptake changes can impair frame alignment; inspect early and low-count frames.

GTM estimates a dense region-to-region spillover matrix by Gaussian-blurring label
indicators in physical space. A nonnegative least-squares solve estimates each
frame's regional activities. Label zero, if present, is fitted as a single uniform
background region over the remaining field of view. Heterogeneous unlabelled tissue
or extracranial uptake can violate that assumption: provide an appropriate complete
regional partition or omit PVC. The image outside the field of view is assumed zero.
At most 512 regions including background are supported; matrices with condition
number greater than 1e6 are rejected. Residuals and condition number are recorded.
Noise amplification, PSF uncertainty, anatomical uncertainty and induced regional
correlations are not propagated into the downstream kinetic confidence intervals.
This method provides regional estimates, not voxelwise corrected PET images.

## Sources and validation scope

PETPrep's parser and workflows for motion, registration, TAC extraction and PVC
were inspected at commit `cf9c8cc9c7d5ddb662c18b2d9653b3ea92c1e58d` as an
architectural reference. No PETPrep implementation was copied into this module.
The algorithms, defaults and scope differ; equivalence or superiority is not claimed.

- [PETPrep source](https://github.com/nipreps/petprep/tree/cf9c8cc9c7d5ddb662c18b2d9653b3ea92c1e58d)
- [N4ITK, Tustison et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3071855/)
- [SimpleITK registration framework](https://simpleitk.readthedocs.io/en/master/registrationOverview.html)
- [Rousset et al., regional partial-volume correction](https://jnm.snmjournals.org/content/39/5/904)

Tests use synthetic data with known tissue classes, bias fields, transforms and
regional activities. End-to-end tests include template-mask propagation, label
sampling, standard-space PET generation, failure recording and unchanged source
timing. These establish implementation behavior for the tested examples, not broad
clinical or tracer-independent performance. See the separate execution records for
the measured-data development checks; independent qualification remains incomplete.
