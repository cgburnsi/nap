# AGENTS.md — How to work in this repo (Python)

This file is a **map + operating rules** for agents working on this codebase.
It is not the full design spec. When in doubt, follow the docs listed below.

## Project intent 
Build a small, correct, legible, modular, Python library that implements the steady-state numerical models described in A.S. Kesten's 1968 paper.

### Milestones
1. Quick, readable conversion of the original document including the Fortran source code
2. Production-grade faithful conversion of the original document including the Fortran source code
3. minimal, testable 1D solver that matches the results from the original paper
4. Minimal, testable 2D solver that matches the results from the original paper

## Execution priority (must follow)
- Work milestone-by-milestone in order; do not start later milestones early unless explicitly approved.
- During milestone 1, optimize for readability and speed of conversion over perfect layout fidelity.
- During milestone 2, optimize for faithful structure/content parity with the original source, with LaTeX as canonical output for publishable PDFs.

## Document conversion contract
### Global conversion defaults
- Convert in reading order from abstract to end of document.
- Use major-section granularity by default.
- When conflicts exist, prioritize `docs/reference/Kesten_1968_original.pdf` and `docs/fortran/*.f` over derived/converted artifacts.

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
- Goal: high-fidelity dual conversion, with LaTeX as canonical project documentation and Markdown kept synchronized.
- Required preservation:
  - heading hierarchy and numbering parity
  - equation symbols and indices exactly preserved
  - table structure/labels/units preserved exactly
  - references/footnotes/cross-references traceable and complete
  - Fortran code blocks aligned with original structure and comments
- Any ambiguity in source text should be called out explicitly in notes rather than silently normalized.
- Milestone-1 exit gate: no unresolved critical ambiguities in source interpretation.

## Numeric tolerance policy
- Use combined threshold: `abs(err) <= max(abs_tol, rel_tol * max(1.0, abs(reference)))`.
- Starter tolerances (intentionally pragmatic for reconstruction):
  - liquid region: `abs_tol=1e-6`, `rel_tol=5e-4`
  - vapor region: `abs_tol=1e-6`, `rel_tol=1e-3`
  - liquid-vapor region: `abs_tol=1e-6`, `rel_tol=2e-3`
- For tabulated outputs, require at least 95% of points to pass per table; log failing rows for review.
- Tightening rule: after baseline parity is reached, reduce relative tolerances by 2x per milestone if tests remain stable.

## Approval policy
- Loose mode: the agent may make small and medium scoped edits autonomously, including adding files and restructuring docs when it improves milestone progress.
- Ask before:
  - deleting files
  - major architecture rewrites
  - dependency or tooling policy changes

## Commit message convention
- Use Conventional Commits: `type(scope): subject`
- Preferred types: `docs`, `feat`, `fix`, `test`, `refactor`, `chore`.

## Source of truth (read these first)
- README.md — how to install/run and what the repo is for
- ARCHITECTURE.md — module boundaries and data ownership rules
- docs/conversion/spec.md — rules for milestones 1-2 conversion outputs
- docs/implementation/spec.md — scope and acceptance for milestones 3-4 solver implementation
- docs/verification/plan.md — verification strategy + acceptance gates
- docs/design/decisions.md — decision log (add entries when changing direction)

If any of these files are missing, create them as small stubs rather than
inventing undocumented architecture.

## Definition of Done (for any change)
A change is “done” only if:
1. It runs without errors for the intended use case.
2. Tests pass (unit + any applicable regression tests).
3. New behavior has tests (or an explicit note why it cannot).
4. Docs are updated when behavior/CLI/file formats change.
5. The change is deterministic (within stated tolerances).
6. No unrelated refactors are mixed in.
7. Dependency rules (below) are respected.
8. Results match the values provided in the Kesten report from 1968.

## Work style rules
### Keep changes small and reviewable
- Prefer multiple small PRs over one large PR.
- Each PR should have a clear goal and short summary.

### Be deterministic
- Avoid nondeterminism. If randomness is required, seed it and document why.
- Numerical outputs should be stable across machines within tolerances.
- If results differ by CPU/BLAS/platform, encode tolerance-based checks.

### Respect module boundaries
- No “god modules.” Keep files focused.
- Avoid circular imports.
- Modules should follow the code layout in the Kesten paper.
- Keep data ownership clear: who constructs, mutates, and consumes arrays.

### Avoid premature complexity
- No GPU/MPI/multiprocessing until verification harness is solid.
- Prefer simple numerics and clean interfaces over cleverness.

## Dependencies (hard rule)
This project must run with only:
- Python standard library
- NumPy
- Matplotlib

No other packages may be added (runtime or development). This includes (but is
not limited to): scipy, numba, pandas, sympy, meshio, shapely, pyvista, tqdm,
rich, click/typer, pydantic, pytest, ruff, black, etc.

### Dependency gate (must follow)
If any task appears to require a new dependency:
1. **Stop immediately.**
2. Write a short note (issue/PR comment/docs snippet) that includes:
   - What feature is blocked
   - Why stdlib + NumPy + Matplotlib are insufficient
   - Two alternative approaches that keep the dependency set unchanged
   - If still justified: the proposed dependency, its scope, and why it’s worth it
3. Do **not** add the dependency or refactor towards it until explicit approval.

### Guidance under this constraint
- Use `argparse` for CLI (stdlib), not click/typer.
- Use stdlib `unittest` for tests (no pytest).
- Use stdlib `json`, `csv`, `pathlib`, `logging` for I/O and instrumentation.
- Keep plotting in `examples/` or `tools/` when possible (avoid importing
  matplotlib in the core solver path unless necessary).

## Repo conventions (expected layout)
(If the repo differs, update this section.)

- src/                       core library code
- app/                       command line entry points (thin wrappers)
- tests/                     unittest tests
- examples/                  runnable examples (kept small)
- docs/                      design + verification notes
- tools/                     helper scripts (goldens/plots/etc.)

## Setup (local)
From the repo root:

```sh
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install numpy matplotlib
```
