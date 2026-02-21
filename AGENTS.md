# AGENTS.md — NAP Python Agent Map

This file is intentionally short. Use it as a navigation map for implementation work.

## Mission
Build and verify a Python implementation of NAP with a stable baseline and incremental refactors.

## Read These First
- `docs/README.md` - docs index
- `docs/implementation/roadmap.md` - milestone order and current focus
- `docs/implementation/architecture.md` - target runtime structure
- `docs/implementation/workflow.md` - coding/review workflow
- `docs/verification/plan.md` - required checks before merge
- `docs/verification/baselines.md` - baseline case definitions and tracked metrics
- `docs/design/decisions.md` - architecture and numerical decisions log

## Working Rules
- Preserve baseline numerical behavior unless intentional and documented.
- Keep changes small and scoped.
- Run verification checks for solver/grid changes.
- Log meaningful direction changes in `docs/design/decisions.md`.

## Code Style Preferences
- Prefer minimal typing; do not add type hints unless explicitly requested or required.
- Dataclass annotations are acceptable when dataclasses are used.
- Prefer straightforward, readable code over abstraction-heavy patterns.
- Keep functions focused and easy to follow.
- Avoid broad rewrites when a local edit is sufficient.

## Work Modes
Use one of these labels in prompts to control execution style.

- `Mode: Baby Steps`
  - Make one small change at a time.
  - Explain the immediate next step before editing.
  - Run minimal checks and stop for confirmation before continuing.

- `Mode: Normal`
  - Complete the requested task end-to-end with standard checks.
  - Batch related edits when they are low risk and easy to review.

- `Mode: Sprint`
  - Move quickly across multiple related edits in one pass.
  - Return a concise summary and key verification results at the end.

## Repo Focus
- `main.py` active implementation entrypoint
- `legacy/code_archive/ver1/` archived code reference
- `legacy/fortran/` historical reference source
