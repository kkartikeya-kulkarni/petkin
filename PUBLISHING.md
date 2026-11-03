# Repository and releases

Repository: https://github.com/kkartikeya-kulkarni/petkin

The initial upload is private. Version 0.3.0rc1 is a research prerelease.
The root README describes installation and usage. Detailed evidence and scope
remain in VALIDATION_STATUS.md and the validation directory.

GitHub Actions builds and tests the installed wheel on Linux and macOS, using
Python 3.10, 3.12 and 3.13. Check the actual run results before making claims about
supported platforms. Local verification is recorded in validation/release/.

Build release files with `python -m build`. Keep source image data and dense
validation inputs outside Git. Source manifests and reproduction instructions
are included; the full imaging evidence bundle remains separately available.

Use prerelease tags while research qualification is incomplete. Preserve all
historical results when adding new evidence. CI success verifies software and
does not replace study-specific validation.
