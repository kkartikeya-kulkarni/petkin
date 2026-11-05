# Third-party provenance

The MIT license covers original petkin code. It does not relicense datasets,
articles, or external packages. No kinfitr or NiftyPAD implementation source is
bundled in this repository; the comparator scripts execute separately obtained
upstream functions/packages. Version pins and source hashes are preserved in
the example manifests and validation records.

- PS13: OpenNeuro ds004868, source dataset declares CC0. The source manifest
  records URLs, file sizes and SHA256 hashes; all 20 PET files matched original
  annex hashes despite renamed participants. This repository includes deidentified
  derived numerical evidence and source metadata. Dense image evidence is separate.
- Earlier FDG/PBR28 studies: provenance and limitations are recorded in the
  original user guide, example source manifests and full evidence archive.
- RIRE: reference only; official data are CC BY 3.0 US. No inaccessible RIRE data
  or reference transforms are represented as included or tested.
- NumPy, SciPy, SimpleITK, matplotlib, nibabel, kinfitr and NiftyPAD are separate
  dependencies or comparison tools governed by their upstream licenses.
- PETPrep (Apache-2.0) was inspected as an architectural reference for native
  preprocessing. No PETPrep source implementation is bundled or called by the
  native pipeline. Reviewed commit and file hashes are recorded under
  `validation/native_preprocessing_v04/UPSTREAM_REVIEW.json`.

There is no endorsement or independent certification by the dataset authors.
Do not assign collaborator identities, institutional affiliations or reviewers
without their agreement.
