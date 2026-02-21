# Implementation Roadmap

## Milestones
1. Establish stable orchestration skeleton (config -> grid -> state -> solver -> verifier -> plots).
2. Rebuild baseline grid/geometry and 1D initialization path.
3. Restore transient marching and boundary treatment.
4. Reproduce baseline reference metrics and plots.
5. Refactor internals while holding parity.

## Current Focus
- Keep geometry+grid coupled until baseline parity is stable.

## Current Status (2026-02-21)
- `main.py` orchestration is active and runnable:
  - config -> grid -> state -> solver -> verifier -> plots -> summary
- Implemented:
  - CD nozzle wall geometry build with computed `XT` and `LT`
  - ONEDIM-style initialization (subsonic-sonic-supersonic branch control)
  - verification against baseline CSV
  - regression no-worsening gates with margin
  - plotting set in `artifacts/verification/`
- Open work:
  - pressure-ratio convention alignment to reference plotting method
  - transient marching/convergence implementation
  - post-ONEDIM solver physics improvements

## Next Baby Steps
1. Align wall pressure ratio convention with reference method and re-check overlay.
2. Add convergence/residual gates once transient marching exists.
3. Decide when to lock absolute thresholds vs relative no-regression gates.
