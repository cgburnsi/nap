import os
import numpy as np


class RunConfig:
    """Top-level runtime settings for one NAP run."""

    def __init__(
        self,
        case_name="baseline_cd_nozzle",
        lmax=21,
        mmax=8,
        nmax=400,
        tconv=5.0e-3,
        fdt=0.25,
        xi=0.31,
        ri=2.50,
        rt=0.80,
        xe=4.05,
        rci=0.80,
        rct=0.50,
        ani_deg=44.88,
        ane_deg=15.0,
        baseline_csv_path="docs/verification/baselines/cd_nozzle_converged_solution.csv",
        xt_expected=2.554,
        xt_tol=0.005,
    ):
        self.case_name = case_name
        self.lmax = lmax
        self.mmax = mmax
        self.nmax = nmax
        self.tconv = tconv
        self.fdt = fdt
        self.xi = xi
        self.ri = ri
        self.rt = rt
        self.xe = xe
        self.rci = rci
        self.rct = rct
        self.ani_deg = ani_deg
        self.ane_deg = ane_deg
        self.baseline_csv_path = baseline_csv_path
        self.xt_expected = xt_expected
        self.xt_tol = xt_tol


class Grid:
    """
    Grid container.
    For now, geometry+grid are intentionally coupled for rapid NAP parity work.
    """

    def __init__(
        self,
        lmax,
        mmax,
        xi=0.31,
        ri=2.50,
        rt=0.80,
        xe=4.05,
        rci=0.80,
        rct=0.50,
        ani_deg=44.88,
        ane_deg=15.0,
    ):
        self.lmax = lmax
        self.mmax = mmax
        self.xi = xi
        self.ri = ri
        self.rt = rt
        self.xe = xe
        self.rci = rci
        self.rct = rct
        self.ani_deg = ani_deg
        self.ane_deg = ane_deg
        self.built = False

    def build(self):
        if self.lmax < 2:
            raise RuntimeError("lmax must be at least 2 for geometry build.")

        ani = np.deg2rad(self.ani_deg)
        ane = np.deg2rad(self.ane_deg)

        self.xw = np.linspace(self.xi, self.xe, self.lmax)
        self.yw = np.zeros(self.lmax)
        self.nxny = np.zeros(self.lmax)

        xtan = self.xi + self.rci * np.sin(ani)
        rtan = self.ri + self.rci * (np.cos(ani) - 1.0)
        rt1 = self.rt - self.rct * (np.cos(ani) - 1.0)
        xt1 = xtan + (rtan - rt1) / np.tan(ani)
        if xt1 < xtan:
            xt1 = xtan
            rt1 = rtan

        self.xt = xt1 + self.rct * np.sin(ani)
        xt2 = self.xt + self.rct * np.sin(ane)

        rt2 = self.rt + self.rct * (1.0 - np.cos(ane))
        self.re = rt2 + (self.xe - xt2) * np.tan(ane)

        self.lt = 0
        for l in range(self.lmax):
            x = self.xw[l]
            if self.xi <= x <= xtan:
                self.yw[l] = self.ri + self.rci * (np.cos(np.arcsin((x - self.xi) / self.rci)) - 1.0)
                self.nxny[l] = (x - self.xi) / (self.yw[l] - self.ri + self.rci)
            elif xtan < x <= xt1:
                self.yw[l] = rt1 + (xt1 - x) * np.tan(ani)
                self.nxny[l] = np.tan(ani)
            elif xt1 < x <= self.xt:
                self.yw[l] = self.rt + self.rct * (1.0 - np.cos(np.arcsin((self.xt - x) / self.rct)))
                self.nxny[l] = (self.xt - x) / (self.rct + self.rt - self.yw[l])
            elif self.xt < x <= xt2:
                self.yw[l] = self.rt + self.rct * (1.0 - np.cos(np.arcsin((x - self.xt) / self.rct)))
                self.nxny[l] = (self.xt - x) / (self.rct + self.rt - self.yw[l])
            elif xt2 < x <= self.xe:
                self.yw[l] = rt2 + (x - xt2) * np.tan(ane)
                self.nxny[l] = -np.tan(ane)

            if l > 0 and self.yw[l] < self.yw[self.lt]:
                self.lt = l

        # Centerbody placeholders for current CD nozzle baseline (none).
        self.xcb = self.xw.copy()
        self.ycb = np.zeros(self.lmax)
        self.nxnycb = np.zeros(self.lmax)

        # Build simple structured coordinates for plotting/debug.
        eta = np.linspace(0.0, 1.0, self.mmax)
        self.x = np.repeat(self.xw[:, None], self.mmax, axis=1)
        self.y = np.zeros((self.lmax, self.mmax))
        for l in range(self.lmax):
            self.y[l, :] = self.ycb[l] + eta * (self.yw[l] - self.ycb[l])

        self.built = True
        return self


class State:
    """Flow/state storage container."""

    def __init__(self, lmax, mmax):
        self.lmax = lmax
        self.mmax = mmax
        self.initialized = False

    def initialize(self, grid):
        if not grid.built:
            raise RuntimeError("Grid must be built before state initialization.")

        if grid.lmax != self.lmax or grid.mmax != self.mmax:
            raise RuntimeError(
                f"State/Grid shape mismatch: state=({self.lmax},{self.mmax}), "
                f"grid=({grid.lmax},{grid.mmax})"
            )

        # NAP-style storage: (L, M, level) with two time levels.
        nlevels = 2
        shape_3d = (self.lmax, self.mmax, nlevels)
        self.u = np.zeros(shape_3d, dtype=float, order="F")
        self.v = np.zeros(shape_3d, dtype=float, order="F")
        self.p = np.zeros(shape_3d, dtype=float, order="F")
        self.rho = np.zeros(shape_3d, dtype=float, order="F")

        # Boundary condition vectors across radial index M.
        self.pt = np.zeros(self.mmax, dtype=float)
        self.tt = np.zeros(self.mmax, dtype=float)
        self.theta = np.zeros(self.mmax, dtype=float)

        self.initialized = True
        return self


class RunResult:
    """Output bundle from a solver run."""

    def __init__(self, converged, iterations, residual_history=None, metadata=None):
        self.converged = converged
        self.iterations = iterations
        self.residual_history = [] if residual_history is None else residual_history
        self.metadata = {} if metadata is None else metadata


class Solver:
    def __init__(self, config: RunConfig):
        self.config = config

    def solve(self, grid, state):
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
    def verify(self, config, grid, result):
        if not os.path.exists(config.baseline_csv_path):
            raise RuntimeError(f"Baseline CSV not found: {config.baseline_csv_path}")
        if not hasattr(grid, "xt"):
            raise RuntimeError("Grid does not define XT; cannot run baseline verification.")

        xt_err = abs(float(grid.xt) - float(config.xt_expected))
        if xt_err > float(config.xt_tol):
            raise RuntimeError(
                f"XT check failed: XT={grid.xt:.6f}, expected={config.xt_expected:.6f}, tol={config.xt_tol:.6f}"
            )

        lt_value = getattr(grid, "lt", None)
        return {
            "checked": True,
            "xt": float(grid.xt),
            "xt_expected": float(config.xt_expected),
            "xt_error": float(xt_err),
            "xt_tol": float(config.xt_tol),
            "lt_report": lt_value,
            "count": 0,
            "u_mae_fps": float("nan"),
            "v_mae_fps": float("nan"),
            "p_mae_psia": float("nan"),
            "mach_mae": float("nan"),
        }


class Plots:
    def __init__(self, grid, result: RunResult, verify_report: dict):
        self.grid = grid
        self.result = result
        self.verify_report = verify_report

    def plot_grid(self) -> None:
        import matplotlib.pyplot as plt

        if not self.grid.built:
            raise RuntimeError("Grid must be built before plotting.")

        os.makedirs("artifacts/verification", exist_ok=True)
        output_path = "artifacts/verification/cd_baseline_grid.png"

        fig, ax = plt.subplots(figsize=(10, 4.8), dpi=140)

        for l in range(self.grid.lmax):
            ax.plot(self.grid.x[l, :], self.grid.y[l, :], color="0.75", lw=0.8, zorder=1)
        for m in range(self.grid.mmax):
            ax.plot(self.grid.x[:, m], self.grid.y[:, m], color="0.75", lw=0.8, zorder=1)

        ax.plot(self.grid.xw, self.grid.yw, color="#d62728", lw=2.0, label="wall", zorder=3)
        ax.plot(self.grid.xcb, self.grid.ycb, color="#1f77b4", lw=2.0, label="centerline", zorder=3)
        ax.plot(self.grid.x[self.grid.lt, :], self.grid.y[self.grid.lt, :], color="#9467bd", ls="--", lw=1.6, label="throat index", zorder=2)

        ax.set_title("CD Nozzle Structured Grid")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.20)
        ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved grid plot: {output_path}")

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
    grid = Grid(
        lmax=config.lmax,
        mmax=config.mmax,
        xi=config.xi,
        ri=config.ri,
        rt=config.rt,
        xe=config.xe,
        rci=config.rci,
        rct=config.rct,
        ani_deg=config.ani_deg,
        ane_deg=config.ane_deg,
    ).build()

    # 3) Build State.
    state = State(lmax=config.lmax, mmax=config.mmax).initialize(grid)

    # 4) Run Solver.
    solver = Solver(config)
    result = solver.solve(grid, state)

    # 5) Run Verifier.
    verifier = Verifier()
    verify_report = verifier.verify(config, grid, result)

    # 6) Run Plots.
    plots = Plots(grid, result, verify_report)
    plots.plot_grid()
    plots.plot_solution()
    plots.plot_residuals()

    # 7) Print concise run summary.
    print("NAP run summary")
    print(f"case={config.case_name}")
    print(f"grid={config.lmax}x{config.mmax}, converged={result.converged}, iterations={result.iterations}")
    print(
        "verify: "
        f"XT={verify_report['xt']:.3f} (target {verify_report['xt_expected']:.3f}, err {verify_report['xt_error']:.3e}), "
        f"LT={verify_report['lt_report']}, "
        f"count={verify_report['count']}, "
        f"U_MAE={verify_report['u_mae_fps']}, "
        f"V_MAE={verify_report['v_mae_fps']}, "
        f"P_MAE={verify_report['p_mae_psia']}, "
        f"Mach_MAE={verify_report['mach_mae']}"
    )
