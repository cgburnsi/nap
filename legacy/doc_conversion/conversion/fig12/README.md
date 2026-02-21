# Figure 12 Extraction Notes

Source: `docs/references/NAP - Two Dimensional Invisid Flow.pdf`

Identified Figure 12 page range for Case No. 3 output:
- Page 36: input deck / run setup summary (Case No. 3)
- Page 37: run progression lines (`N=90` ... `N=340`)
- Page 38: solution surface table (No. 345, first block; `L=1..`)
- Page 39: solution surface table (No. 345 continuation; `L=10..18`)
- Page 40: solution surface table continuation (`L=19..27`)
- Page 41: solution surface table continuation (`L=28..31`) + mass/thrust summary

Interpretation note:
- OCR reads pages 40-41 headers as `NO. 395`, but this is likely a scan/OCR artifact.
- Evidence: `L` continues monotonically from page 39 to 41, and the reported times are effectively the same (`.00237634` vs `.00237639`), consistent with one continued printout rather than a new surface block.

Raw extracted text files:
- `fig12_p36.txt`
- `fig12_p37.txt`
- `fig12_p38.txt`
- `fig12_p39.txt`
- `fig12_p40.txt`
- `fig12_p41.txt`

These are `pdftotext -layout` outputs and preserve table spacing reasonably well for follow-on parsing into CSV/Markdown tables.

Structured table artifact:
- `fig12_surface345_longtable.tex` (OCR-normalized full 31 x 6 grid)
- Current row coverage: 186/186 expected rows (31 x 6 grid)
- Some rows on pages 40--41 required OCR interpretation for ambiguous digits/symbols; keep the raw page extracts for auditability.
