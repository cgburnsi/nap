# Current State Snapshot (2026-02-21)

## Branch / Repo State
- Active branch: `feat/regression-rate-limits`
- Main branch contains merged work through commit: `0343089`
- This branch includes additional local updates not yet merged.

## What is Working
- `main.py` runs end-to-end.
- Grid geometry for baseline CD nozzle is computed from parameters.
- ONEDIM initialization runs with branch-safe area-Mach solve.
- Verification checks run against baseline CSV and enforce hard gates.
- Plots are generated in `artifacts/verification/`:
  - `cd_baseline_grid.png`
  - `cd_baseline_solution.png`
  - `cd_baseline_onedim_2d.png`
  - `cd_baseline_mach_contours.png`
  - `cd_baseline_wall_pressure_ratio.png`
  - `cd_baseline_residual.png`

## Active Verification/Gating Behavior
- Hard `XT` check.
- Hard inlet Mach mode check for subsonic-sonic-supersonic mode.
- Hard CSV domain validity check.
- Hard finite-metric check.
- Hard no-regression metric gates for:
  - `mach_mae`
  - `p_mae_psia`
  - `wall_p_ratio_mae`

## Known Gaps
- No transient marching implementation yet.
- Residual history is still placeholder-level.
- Pressure ratio convention alignment may still need refinement.

## Next Recommended Baby Steps
1. Pressure ratio alignment with reference convention and re-check overlay.
2. Add residual/convergence gates after transient marching exists.
3. Decide when to convert relative no-regression gates into absolute target gates.

## Resume Checklist
1. `git checkout feat/regression-rate-limits`
2. `python3 main.py`
3. Inspect `artifacts/verification/*.png`
4. Confirm verification summary and gate lines in console output
