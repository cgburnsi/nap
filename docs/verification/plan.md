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
- CSV domain gate:
  - raw reference points must satisfy `0 <= y <= y_wall(x)` within tolerance.
- Metric sanity gate:
  - computed verification metrics must be finite.
- Regression gates ("no regression + margin"):
  - `mach_mae <= baseline_mach_mae * (1 + margin)`
  - `p_mae_psia <= baseline_p_mae_psia * (1 + margin)`
  - `wall_p_ratio_mae <= baseline_wall_p_ratio_mae * (1 + margin)`

## Notes
- Current margin is configured in code (`regression_margin_frac`) and defaults to 20%.
- Residual/convergence hard gates will be added after transient marching is implemented.
