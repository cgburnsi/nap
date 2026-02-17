# Conversion Worklog

## 2026-02-17

### Completed
- Added missing conversion/design source-of-truth docs:
  - `docs/conversion/spec.md`
  - `docs/design/decisions.md`
- OCR cleanup pass applied to Fortran listings used by `docs/NAP.tex`:
  - `legacy/fortran/map.f` (replaced with cleaner in-repo version)
  - `legacy/fortran/inlet.f`
  - `legacy/fortran/onedim.f`
  - `legacy/fortran/wall.f`
  - `legacy/fortran/fortran_main.f`
  - `legacy/fortran/exitt.f` (provisional reconstruction from PDF text blocks)
  - `legacy/fortran/shock.f` (provisional reconstruction from PDF text blocks)
  - `legacy/fortran/plot.f` (provisional reconstruction from PDF text blocks)
  - `legacy/fortran/geomcb.f` (provisional reconstruction from PDF text blocks)
  - `legacy/fortran/mtlup.f` (provisional reconstruction from PDF text blocks)
  - `legacy/fortran/diff.f` (provisional reconstruction from PDF text blocks)

### Remaining OCR-heavy Files
- None currently tracked as unresolved OCR-heavy in Appendix C Fortran set.

### Notes
- OCR cleanup required line-by-line source interpretation from the reference PDF for multiple Appendix C listings.
- `exitt.f` and `shock.f` were upgraded from raw OCR fragments to readable reconstructions and aligned to Appendix C line blocks in their degraded equation regions.
- `shock.f` pressure-viscosity energy term was re-aligned to Appendix C line blocks and its prior in-file ambiguity marker was removed.
- `plot.f` was upgraded from partial OCR transcription to a readable reconstruction with restored contour-cell intersection logic from Appendix C line blocks.
- `geomcb.f` was upgraded from raw OCR text to a readable reconstruction; minor formula and format literal normalization was required where scan quality was poor.
- `mtlup.f` and `diff.f` were upgraded from raw OCR fragments to readable reconstructions based on Appendix C; comments/diagnostic literals were normalized where scan artifacts obscured text.
- `mtlup.f` second-order `N=2` fallback was normalized to first-order interpolation flow; diagnostic format text was normalized for readability.
- Legacy 2024 milestone snapshot files were annotated as historical; active status now points to conversion trackers.
- For milestone 1, preference is to preserve meaning and mark unresolved ambiguities explicitly rather than invent uncertain code tokens.

### Milestone 1 Exit Gate
- Decision: `PASS` for Milestone 1 readability gate.
- Rationale: no unresolved critical OCR ambiguities remain in Appendix C listings.
- Remaining non-critical normalization choices:
  - `geomcb.f`: some long FORMAT literals were wording-normalized from degraded scan regions.
  - `mtlup.f`: diagnostic FORMAT wording normalized from degraded scan regions.
  - Multiple listings retain conservative token normalization for readability (spacing/punctuation/case) where source glyphs were corrupted.
