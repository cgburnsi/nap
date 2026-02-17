# Decision Log

## 2026-02-17
- Keep `docs/NAP.tex` as canonical conversion output for publication-oriented PDF builds.
- Keep Fortran Appendix listings sourced from `legacy/fortran/*.f` so OCR cleanup is centralized in one place.
- Apply conservative OCR fixes first (high-confidence token/identifier corrections), and record unresolved ambiguity instead of speculative rewrites.
