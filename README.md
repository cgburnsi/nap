# nap

A faithful document conversion project: reproducing Michael C. Cline's NAP (Nozzle Analysis Program, January 1977) in human-readable Markdown and LaTeX.

## Current phase: Document Conversion
We are focused **exclusively** on converting the original NAP document and its associated Fortran source code into accessible, readable formats.

**Milestone 1:** Quick, readable Markdown conversion  
**Milestone 2:** Production-grade faithful conversion with LaTeX as canonical format

After document conversion is complete and verified, solver implementation will begin in a separate roadmap.

## Current phase workflow
- Conversion rules: `docs/conversion/spec.md`
- Conversion style guide: `docs/conversion/style-guide.md`
- Conversion source mapping: `docs/conversion/source-map.md`
- Conversion acceptance criteria: `docs/conversion/acceptance.md`
- Conversion progress log: `docs/conversion/worklog.md`
- Conversion QA checklist: `docs/conversion/qa_checklist.md`
- Implementation rules (milestones 3-4): `docs/implementation/spec.md`
- Implementation module plan: `docs/implementation/module-plan.md`
- Implementation test matrix: `docs/implementation/test-matrix.md`
- Verification tolerances: `docs/verification/tolerances.md`
- Verification golden index: `docs/verification/golden-index.md`
- Verification regression protocol: `docs/verification/regression-protocol.md`
- Vapor port blockers: `docs/verification/fortran-vapor-port-blockers.md`
- Decision history: `docs/design/decisions.md`

## Repository layout
- `docs/references/` - original reference documents
- `docs/conversion/` - conversion rules, specs, and progress tracking
- `docs/design/` - design decisions and architectural notes
- `legacy/` - early drafts and reference implementations

## Local setup
```sh
python -m venv .venv
source .venv/bin/activate
# Clone the repo
git clone <repo-url>
cd nap

# No dependencies required for document conversion
# (Python not needed unless code is generated)
## Run baseline regression example
```sh
python examples/ex02_baseline_regression.py
```

## Run liquid physics slice
```sh
python app/run_solver.py --mode physics --region liquid
python app/run_solver.py --mode physics --region vapor
python app/run_solver.py --mode physics --region liquid_vapor
```

## Plot outputs
```sh
# Interactive plot window
python app/run_solver.py --mode physics --region liquid --plot

# Save PNG
python app/run_solver.py --mode physics --region liquid --plot-output artifacts/plots/liquid_physics.png
python app/run_solver.py --mode iterate --initial 0 --target 1 --max-iterations 50 --plot-output artifacts/plots/iterate_residual.png

# Compare physics vs baseline in one plot (especially useful for vapor single-point case)
python app/run_solver.py --mode physics --region vapor --plot --plot-compare
python app/run_solver.py --mode physics --region vapor --plot-output artifacts/plots/vapor_compare.png --plot-compare

# Temperature vs catalyst bed length (full-bed profile)
python app/run_solver.py --mode physics --region liquid --plot-temp-bed --plot
python app/run_solver.py --mode physics --region liquid --plot-temp-bed --plot-output artifacts/plots/temp_vs_bed_physics.png
python app/run_solver.py --mode physics --region liquid --plot-temp-bed --plot-compare --plot-output artifacts/plots/temp_vs_bed_compare.png
python app/run_solver.py --mode physics --region liquid --plot-temp-bed --temp-bed-vapor-model fortran_inspired --plot-output artifacts/plots/temp_vs_bed_fortran_inspired.png

# Use your reference curve file directly (model/reference/both)
python app/run_solver.py --mode physics --region liquid --plot-temp-bed --temp-curve-source reference --plot-output artifacts/plots/temp_vs_bed_reference.png
python app/run_solver.py --mode physics --region liquid --plot-temp-bed --temp-curve-source both --plot-output artifacts/plots/temp_vs_bed_model_vs_reference.png
```

Note: `--plot-temp-bed` now uses a dedicated full-bed equation integration model for `physics` source, rather than stitching sparse region outputs.
WARNING: current `--mode physics` outputs are experimental equation-model results with simplifying assumptions and are not yet fully validated. Use `--temp-bed-vapor-model fortran_inspired` to run the new coupled vapor-balance prototype.

## Run regression gate
```sh
python app/run_solver.py --mode regress --region liquid --source baseline
python app/run_solver.py --mode regress --region liquid --source physics
python app/run_solver.py --mode regress --region vapor --source physics
python app/run_solver.py --mode regress --region liquid_vapor --source physics
```

## Build LaTeX PDF
```sh
tools/build_pdf.sh
```
- Build details and options: `docs/latex/BUILD.md`

## Notes
- Dependencies are intentionally restricted to Python stdlib + NumPy + Matplotlib.
- Keep architectural decisions in `docs/design/decisions.md`.
- Keep verification gates and acceptance criteria in `docs/verification/plan.md`.
