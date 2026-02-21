# Decision Log

Record significant architecture and numerical choices here with date and rationale.

## 2026-02-21
- Keep `main.py` as the active implementation entrypoint to support Spyder workflow and variable explorer visibility.
- Keep geometry and grid coupled in `Grid` during early parity work; defer separation until baseline behavior is stable.
- Use a modern regression harness layer (in addition to NAP-style convergence logic) for development safety:
  - hard `XT` gate against baseline expectation
  - hard CSV domain validity check (`0 <= y <= y_wall(x)` with tolerance)
  - hard "no-regression" metric gates with configurable margin
- Adopt stepwise implementation policy:
  - establish plots/verification first
  - then improve numerical methods incrementally
  - keep each change as a small baby step
