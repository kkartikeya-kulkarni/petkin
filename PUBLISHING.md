# Repository and releases

Repository: https://github.com/kkartikeya-kulkarni/petkin

The repository is private. Version 0.4.0.dev1 is the active regional workflow on
`main` and the only downloadable GitHub release. The old 0.3.0rc1 release and tag
were deleted at the owner's request. Historical code and validation evidence remain
in Git history and the validation directory.
The root README describes installation and usage. Detailed evidence and scope
remain in VALIDATION_STATUS.md and the validation directory.

GitHub Actions builds and tests the installed wheel on Linux with Python 3.10,
3.12 and 3.13, and macOS Intel with Python 3.12. The initial macOS ARM runner
queue failed to acquire two runners; this was infrastructure failure, not test
failure. Local ARM macOS installed-wheel testing is also recorded. Check the actual run results before making claims about
supported platforms. Local verification is recorded in validation/release/.

Build release files with `python -m build`. Keep source image data and dense
validation inputs outside Git. Source manifests and reproduction instructions
are included; the full imaging evidence bundle remains separately available.

Use prerelease tags while research qualification is incomplete. Preserve all
historical results when adding new evidence. CI success verifies software and
does not replace study-specific validation.
