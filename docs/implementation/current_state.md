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
- Hard residual sanity check (non-empty, finite, non-negative residual history).
- Residual shape criteria are configured and auto-activate when history length is sufficient:
  - `residual_trend_min_points = 10`
  - `residual_end_over_start_max = 0.50`
  - `residual_nonincreasing_frac_min = 0.70`
- Hard CSV domain validity check.
- Hard finite-metric check.
- Hard no-regression metric gates for:
  - `mach_mae`
  - `p_mae_psia`
  - `wall_p_ratio_mae`

## Known Gaps
- No transient marching implementation yet.
- Residual history is still placeholder-level (shape criteria configured, but currently not active with one-point history).
- Pressure ratio convention alignment to legacy/reference plotting method may still need refinement if new reference evidence appears.
- Current runtime tracks two wall pressure-ratio diagnostics:
  - canonical/gated: `P/Pt`
  - diagnostic-only: `P/P_wall,inlet`
- 2026-02-21 method check result:
  - `P/Pt_MAE = 0.030674061764350037`
  - `P/P_wall,inlet_MAE = 0.030011552085192968`
  - relative MAE improvement of alternate vs canonical: `2.16%`
  - decision threshold for switching was `>= 5%`; recommendation: keep canonical `P/Pt`

## Next Recommended Baby Steps
1. Keep `P/Pt` as canonical unless stronger reference-method evidence suggests a change.
2. Add residual/convergence gates after transient marching exists.
3. Decide when to convert relative no-regression gates into absolute target gates.

## Resume Checklist
1. `git checkout feat/regression-rate-limits`
2. `python3 main.py`
3. Inspect `artifacts/verification/*.png`
4. Confirm verification summary and gate lines in console output
