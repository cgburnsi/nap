# Verification Plan

## Required Checks
- Grid sanity: dimensions, throat index/location.
- Initialization sanity: inlet/throat/exit Mach behavior.
- March sanity: residual trend and convergence status.
- Reference check: comparison against published/table data when available.

## Pass/Fail Gate
A change is acceptable only if baseline behavior is preserved or deviations are explicitly documented.

## Active Hard Gates (Current)
- Baseline CSV file must exist and have required columns.
- Geometry gate:
  - `abs(XT - XT_expected) <= XT_tol`
- ONEDIM mode gate (for current baseline mode):
  - inlet Mach must remain subsonic
- Residual sanity gate:
  - residual history must be non-empty, finite, and non-negative
- Residual shape scaffold (activates when history has more than one point):
  - if residual history length is less than `residual_trend_min_points`, report as insufficient for convergence-shape gating
  - once length is at least `residual_trend_min_points`, require:
    - `residual_end / residual_start <= residual_end_over_start_max`
    - fraction of non-increasing residual steps `>= residual_nonincreasing_frac_min`
- CSV domain gate:
  - raw reference points must satisfy `0 <= y <= y_wall(x)` within tolerance.
- Metric sanity gate:
  - computed verification metrics must be finite.
- Regression gates ("no regression + margin"):
  - `mach_mae <= baseline_mach_mae * (1 + margin)`
  - `p_mae_psia <= baseline_p_mae_psia * (1 + margin)`
  - `wall_p_ratio_mae <= baseline_wall_p_ratio_mae * (1 + margin)` (`P/Pt` convention)

## Notes
- Current margin is configured in code (`regression_margin_frac`) and defaults to 20%.
- Default residual-shape thresholds are configured in code:
  - `residual_trend_min_points = 10`
  - `residual_end_over_start_max = 0.50`
  - `residual_nonincreasing_frac_min = 0.70`
- `wall_p_ratio_mae` is currently defined with `P/Pt` as the canonical baseline convention.
- `P/P_wall,inlet` is tracked as diagnostic-only until reference convention alignment is finalized.
