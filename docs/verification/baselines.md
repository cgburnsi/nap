# Baseline Cases

## Case: Baseline CD nozzle
- Grid: `21 x 8`
- Purpose: baseline regression for architecture and solver changes
- Track:
  - `LT`, `XT`
  - residual history
  - key MAE/MAPE values against reference table

### Reference Assets
- `docs/references/NAP_original_scan.pdf`
- `docs/references/NAP_latex_conversion.pdf`
- `docs/verification/baselines/cd_nozzle_converged_solution.csv`
- `docs/verification/baselines/cd_nozzle_wall_from_scan.png`
- `docs/verification/baselines/cd_nozzle_mach_pressure_reference.png`

### Current Baseline Snapshot (for no-regression gates)
- `XT_expected = 2.554`
- `XT_tol = 0.005`
- `baseline_mach_mae = 0.07736578378450953`
- `baseline_p_mae_psia = 1.8779507955979222`
- `baseline_wall_p_ratio_mae = 0.030674061764350037`
- `regression_margin_frac = 0.20`
