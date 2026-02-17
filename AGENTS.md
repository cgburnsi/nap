# AGENTS.md — How to work in this repo

This file is a **map + operating rules** for agents working on this document conversion project.
It is not the full design spec. When in doubt, follow the docs listed below.

## Project intent 
Convert the original NAP document and Fortran source code into human-readable Markdown and production-grade LaTeX formats.

### Milestones
1. Quick, readable conversion of the original document including the Fortran source code
2. Production-grade faithful conversion of the original document including the Fortran source code

## Execution priority (must follow)
- Work milestone-by-milestone in order; do not start coding until document conversion is complete and approved.
- During milestone 1, optimize for readability and speed of conversion over perfect layout fidelity.
- During milestone 2, optimize for faithful structure/content parity with the original source, with LaTeX as canonical output for publishable PDFs.
- After document conversion is complete, solver implementation milestones will be defined in a separate roadmap.

## Document conversion contract

### Global conversion defaults
- Convert in reading order from abstract to end of document.
- Use major-section granularity by default.
- When conflicts exist, prioritize `docs/references/NAP - Two Dimensional Invisid Flow.pdf` over derived/converted artifacts.

### Milestone 1 (quick readable conversion)
- Goal: human-readable Markdown that preserves technical meaning.
- Required preservation:
  - section headings and section order
  - equations (can be simplified formatting if mathematically equivalent)
  - tables (can be simplified if values/labels remain accurate)
  - figure captions and references
  - Fortran listings with line integrity preserved
- Allowed simplifications:
  - minor typography/layout differences
  - normalized whitespace and line wrapping
  - modernized Markdown structure for readability

### Milestone 2 (production-grade faithful conversion)
- Goal: high-fidelity LaTeX conversion as canonical format, with synchronized Markdown.
- Required preservation:
  - heading hierarchy and numbering parity
  - equation symbols and indices exactly preserved
  - table structure/labels/units preserved exactly
  - references/footnotes/cross-references traceable and complete
  - Fortran code blocks aligned with original structure and comments
- Any ambiguity in source text should be called out explicitly in notes rather than silently normalized.
- Milestone-1 exit gate: no unresolved critical ambiguities in source interpretation.

## Approval policy
- Loose mode: the agent may make small and medium scoped edits autonomously, including adding files and restructuring docs when it improves milestone progress.
- Ask before:
  - deleting files
  - major restructuring of conversion outputs
  - scope/goal changes for the document conversion

## Commit message convention
- Use Conventional Commits: `type(scope): subject`
- Preferred types: `docs`, `refactor`, `fix`.

## Source of truth (read these first)
- README.md — project overview and conversion goals
- docs/conversion/spec.md — rules for milestones 1-2 NAP document conversion
- docs/design/decisions.md — decision log (add entries when changing direction)

If any of these files are missing, create them as small stubs.

## Definition of Done (for conversion work)
A documentation change is "done" only if:
1. It preserves the original meaning and structure (per milestone contract).
2. Technical content (equations, tables, Fortran code) is accurate.
3. Cross-references and citations are traceable to the original source.
4. Formatting is clean and consistent with the rest of the conversion.
5. No unrelated refactors are mixed in.

## Work style rules

### Keep changes small and reviewable
- Prefer multiple small PRs over one large PR.
- Each PR should have a clear goal and short summary.

### Preserve source fidelity
- When in doubt, match the original structure and layout.
- Call out deviations explicitly in notes or comments.
- Track all changes to section order or content omissions.

## Repo conventions (expected layout for document conversion phase)

- docs/references/       — original reference documents
- docs/conversion/       — conversion rules, specs, and progress tracking
- docs/design/           — design decisions and architectural notes
- legacy/                — early drafts and reference implementations
