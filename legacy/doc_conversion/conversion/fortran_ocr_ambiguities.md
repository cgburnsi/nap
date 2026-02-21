# Fortran OCR Ambiguities

This file tracks Appendix C Fortran listings where the scanned source is too degraded for high-confidence transcription.

## High Confidence (normalized)
- `legacy/fortran/map.f`
- `legacy/fortran/fortran_main.f` (partial normalization)
- `legacy/fortran/inlet.f` (partial normalization)
- `legacy/fortran/wall.f` (partial normalization)
- `legacy/fortran/onedim.f` (partial normalization)

## Medium Confidence (reconstructed from PDF blocks)
- `legacy/fortran/exitt.f` (reconstructed and line-aligned to Appendix C compatibility-equation block)
- `legacy/fortran/shock.f` (reconstructed and line-aligned to Appendix C viscosity-energy block)
- `legacy/fortran/plot.f` (reconstructed contour-cell intersection logic from Appendix C line blocks)
- `legacy/fortran/geomcb.f` (reconstructed from Appendix C blocks; minor geometric/formatted text normalization)
- `legacy/fortran/mtlup.f` (reconstructed interpolation logic; diagnostic format text normalized)
- `legacy/fortran/diff.f` (reconstructed finite-difference logic; comments normalized)

## Low Confidence (heavy OCR corruption)
- None currently tracked.

## Gate Summary
- Milestone 1 critical-ambiguity gate: `PASS`
- Remaining notes are non-critical wording/format normalizations in degraded scan regions.

## Ambiguity Types Observed
- Lost or substituted tokens (`0/O`, `1/l`, punctuation and operators).
- Broken continuation records and line labels.
- Truncated lines and mixed line-number gutters in OCR text.
- Common block name/array shape corruption.

## Conversion Policy for These Files
- Do not silently invent uncertain code tokens.
- Keep current OCR transcription visible for traceability.
- Add explicit OCR notes near listings in `docs/NAP.tex`.
- Reconstruct only with page-level verification against the original PDF.
