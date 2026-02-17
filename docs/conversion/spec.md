# NAP Conversion Spec

## Scope
- Convert `docs/references/NAP - Two Dimensional Invisid Flow.pdf` into readable Markdown and LaTeX.
- During conversion milestones, prioritize document fidelity over solver implementation.

## Milestone 1 (Readable Conversion)
- Preserve section order and technical meaning.
- Preserve equations, tables, figure captions/references, and Fortran listings.
- Allow formatting simplifications that improve readability.

## Milestone 2 (Faithful Conversion)
- LaTeX is canonical output.
- Preserve numbering hierarchy, symbols/indices, and table structure exactly.
- Track all ambiguities explicitly instead of silently normalizing.

## Current Working Rule
- Source-of-truth for content disputes: `docs/references/NAP - Two Dimensional Invisid Flow.pdf`.
