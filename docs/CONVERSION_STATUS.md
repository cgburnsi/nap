# NAP Document Conversion

## 2026-02-17 Update

- Appendix C listings are now split into two confidence tiers:
  - normalized OCR listings with conservative fixes (`map.f`, `fortran_main.f`, `inlet.f`, `wall.f`, `onedim.f`)
  - reconstructed-from-PDF listings with explicit ambiguity/normalization notes (`exitt.f`, `shock.f`, `plot.f`, `geomcb.f`, `mtlup.f`, `diff.f`)
- OCR ambiguity notes are now explicit in `docs/NAP.tex` at each affected listing section.
- Detailed ambiguity tracker added: `docs/conversion/fortran_ocr_ambiguities.md`.
- Conversion/design source-of-truth stubs added:
  - `docs/conversion/spec.md`
  - `docs/design/decisions.md`
  - `docs/conversion/worklog.md`
- Milestone 1 exit-gate assessment recorded as `PASS` in `docs/conversion/worklog.md`:
  - no unresolved critical OCR ambiguities remain in Appendix C listings
  - remaining notes are non-critical wording/format normalizations in degraded diagnostic FORMAT text

## Current Status: Milestone 1 (Draft LaTeX Template)

This directory contains the ongoing conversion of Michael C. Cline's 1977 NAP (Nozzle Analysis Program) document from PDF to LaTeX and Markdown formats.

## Files

### Main Conversion Files
- **NAP.tex** — Complete LaTeX template (29 KB, 23 pages when compiled)
- **NAP-Milestone1.pdf** — Current compiled PDF output (259 KB)
- **NAP.md** — Markdown conversion (in progress)

### References
- **../references/NAP - Two Dimensional Invisid Flow.pdf** — Original source document (1977)

### Figures
- **figures/** — Placeholder images for nozzle geometries and Mach contours
  - All figures are currently placeholders marked with "[Original image required]"
  - These will be replaced with actual extracted/scanned images in Milestone 2

### Source Code
- **../nap/fortran code/** — Original Fortran source files referenced in Appendix C
  - All Fortran files are cataloged in NAP.tex appendices
  - Files: fortran_main.f, geom.f, inlet.f, inter.f, masflo.f, onedim.f, wall.f, etc.

## Conversion Progress (Current Snapshot)

### Done in current branch
- Appendix C listings are now present in `docs/NAP.tex` with per-listing OCR notes.
- OCR cleanup completed for normalized files:
  - `map.f`, `fortran_main.f`, `inlet.f`, `wall.f`, `onedim.f`
- OCR reconstruction completed for degraded files:
  - `exitt.f`, `shock.f`, `plot.f`, `geomcb.f`, `mtlup.f`, `diff.f`
- Ongoing decisions and ambiguity tracking now live in:
  - `docs/conversion/worklog.md`
  - `docs/conversion/fortran_ocr_ambiguities.md`
  - `docs/design/decisions.md`

### Remaining milestones for conversion
- Validate all medium-confidence reconstructed listings against the PDF at equation/token level.
- Synchronize Markdown output with canonical LaTeX content.
- Replace placeholder figures for milestone-2 fidelity.

## What's Needed From You

1. **Original Images** — Do you have scanned/extracted versions of:
   - Nozzle geometry diagrams (Figures 1a, 1b, 1c)
   - Mach number contour plots (Figures 2, 3, 4)
   - Any other figures in the original PDF?

2. **Fortran Code** — The source files are in `../nap/fortran code/`. Should we:
   - Pretty-print them directly in the appendices?
   - Or include them as-is with minimal formatting?

3. **Section II Content** — This requires careful extraction from the PDF. Do you have:
   - Typed notes on subroutine descriptions?
   - Input/output format documentation?
   - Or should we OCR/extract from the PDF?

4. **References** — Should we look up the complete citation information for all 18 references, or just include what's available?

## Building the PDF

To recompile the LaTeX document:

```bash
cd /home/cgburnsi/code/nap/docs
pdflatex -interaction=nonstopmode -output-directory=.build NAP.tex
cp .build/NAP.pdf NAP-Milestone1.pdf
```

Requires: `pdflatex`, `graphicx` package (standard TeX Live installation)

## Milestone 2 (Production-Grade)

Once all content is complete, Milestone 2 will focus on:
- High-fidelity equation formatting (verify all symbols, indices, spacing)
- Replace placeholder figures with actual images
- Complete all references with full citations
- Professional PDF layout with proper bookmarks and cross-references
- Synchronized Markdown version for web/GitHub readability

## Naming Convention

- **NAP.tex** — Source (always maintained)
- **NAP-Milestone1.pdf** — Draft output
- **NAP-Final.pdf** — Milestone 2 output (when ready)

## Conversion Notes

- Original document: 21,174 lines (extracted by pdftotext)
- LaTeX template: 23 pages
- Font: 11pt Computer Modern (standard)
- Language: English, original 1977 technical documentation
- All equations preserved in LaTeX format
- All table data preserved accurately
