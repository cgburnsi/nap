from dataclasses import dataclass, field


@dataclass
class RunConfig:
    """Top-level runtime settings for one NAP run."""

    case_name: str = "baseline_cd_nozzle"
    lmax: int = 21
    mmax: int = 8
    nmax: int = 400
    tconv: float = 5.0e-3
    fdt: float = 0.25


@dataclass
class Grid:
    """
    Grid container.
    For now, geometry+grid are intentionally coupled for rapid NAP parity work.
    """

    lmax: int
    mmax: int
    built: bool = False

    def build(self) -> "Grid":
        # Placeholder: geometry contour + mesh generation goes here.
        self.built = True
        return self


@dataclass
class State:
    """Flow/state storage container."""

    lmax: int
    mmax: int
    initialized: bool = False

    def initialize(self, grid: Grid) -> "State":
        if not grid.built:
            raise RuntimeError("Grid must be built before state initialization.")
        # Placeholder: allocate/initialize u,v,p,rho,pt,tt,theta arrays here.
        self.initialized = True
        return self


@dataclass
class RunResult:
    """Output bundle from a solver run."""

    converged: bool
    iterations: int
    residual_history: list[float] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


class Solver:
    def __init__(self, config: RunConfig):
        self.config = config

    def solve(self, grid: Grid, state: State) -> RunResult:
        if not grid.built or not state.initialized:
            raise RuntimeError("Grid and state must be prepared before solve().")
        # Placeholder: interior march + boundary routines + convergence check.
        return RunResult(
            converged=False,
            iterations=0,
            residual_history=[],
            metadata={"case_name": self.config.case_name},
        )


class Verifier:
    def verify(self, result: RunResult) -> dict:
        # Placeholder: Table 2.2 / reference comparisons go here.
        return {
            "checked": True,
            "count": 0,
            "u_mae_fps": float("nan"),
            "v_mae_fps": float("nan"),
            "p_mae_psia": float("nan"),
            "mach_mae": float("nan"),
        }


class Plots:
    def __init__(self, result: RunResult, verify_report: dict):
        self.result = result
        self.verify_report = verify_report

    def plot_grid(self) -> None:
        # Placeholder: grid/boundary plot.
        return

    def plot_solution(self) -> None:
        # Placeholder: solution/diagnostic plot.
        return

    def plot_residuals(self) -> None:
        # Placeholder: residual history plot.
        return





if __name__ == "__main__":
    ''' NAP: Converging-Diverging Nozzle Calculations '''
    
    
    # 1) Build config.
    config = RunConfig()

    # 2) Build Grid.
    grid = Grid(lmax=config.lmax, mmax=config.mmax).build()

    # 3) Build State.
    state = State(lmax=config.lmax, mmax=config.mmax).initialize(grid)

    # 4) Run Solver.
    solver = Solver(config)
    result = solver.solve(grid, state)

    # 5) Run Verifier.
    verifier = Verifier()
    verify_report = verifier.verify(result)

    # 6) Run Plots.
    plots = Plots(result, verify_report)
    plots.plot_grid()
    plots.plot_solution()
    plots.plot_residuals()

    # 7) Print concise run summary.
    print("NAP run summary")
    print(f"case={config.case_name}")
    print(f"grid={config.lmax}x{config.mmax}, converged={result.converged}, iterations={result.iterations}")
    print(
        "verify: "
        f"count={verify_report['count']}, "
        f"U_MAE={verify_report['u_mae_fps']}, "
        f"V_MAE={verify_report['v_mae_fps']}, "
        f"P_MAE={verify_report['p_mae_psia']}, "
        f"Mach_MAE={verify_report['mach_mae']}"
    )
