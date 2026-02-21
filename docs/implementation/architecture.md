# Implementation Architecture

## Runtime Flow
`RunConfig -> Grid -> State -> Solver -> RunResult -> Verifier/Plots`

## Responsibilities
- `Grid`: geometry, mesh, metrics, boundary index maps
- `State`: flow arrays and time levels
- `Solver`: initialization + marching + convergence checks
- `Verifier`: reference comparisons
- `Plots`: diagnostics and report figures
