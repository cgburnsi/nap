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
        verify_level=1,
        xt_expected=2.554,
        xt_tol=0.005,
        gamma=1.4,
        rgas=53.35,
        pt_in_psia=70.0,
        tt_in_f=80.0,
        onedim_mode="subsonic_sonic_supersonic",
        csv_domain_tol=1.0e-3,
        regression_margin_frac=0.20,
        baseline_mach_mae=0.07736578378450953,
        baseline_p_mae_psia=1.8779507955979222,
        baseline_wall_p_ratio_mae=0.030674061764350037,
        residual_trend_min_points=10,
        residual_end_over_start_max=0.50,
        residual_nonincreasing_frac_min=0.70,
        use_synthetic_residual_history=False,
        synthetic_residual_points=10,
        enable_transient_stub=True,
        enable_inlet_bc=False,
        enable_wall_bc=False,
        enable_exit_bc=False,
        exit_back_pressure_psia=None,
        transient_stub_dt_base=1.0e-4,
        transient_cfl_target=0.005,
        enable_transient_logging=True,
        transient_log_every=1,
        transient_u_rel_jump_max=5.0,
        transient_v_rel_jump_max=5.0,
        transient_p_rel_jump_max=5.0,
        enable_plots=True,
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
        self.verify_level = verify_level
        self.xt_expected = xt_expected
        self.xt_tol = xt_tol
        self.gamma = gamma
        self.rgas = rgas
        self.pt_in_psia = pt_in_psia
        self.tt_in_f = tt_in_f
        self.onedim_mode = onedim_mode
        self.csv_domain_tol = csv_domain_tol
        self.regression_margin_frac = regression_margin_frac
        self.baseline_mach_mae = baseline_mach_mae
        self.baseline_p_mae_psia = baseline_p_mae_psia
        self.baseline_wall_p_ratio_mae = baseline_wall_p_ratio_mae
        self.residual_trend_min_points = residual_trend_min_points
        self.residual_end_over_start_max = residual_end_over_start_max
        self.residual_nonincreasing_frac_min = residual_nonincreasing_frac_min
        self.use_synthetic_residual_history = use_synthetic_residual_history
        self.synthetic_residual_points = synthetic_residual_points
        self.enable_transient_stub = enable_transient_stub
        self.enable_inlet_bc = enable_inlet_bc
        self.enable_wall_bc = enable_wall_bc
        self.enable_exit_bc = enable_exit_bc
        self.exit_back_pressure_psia = exit_back_pressure_psia
        self.transient_stub_dt_base = transient_stub_dt_base
        self.transient_cfl_target = transient_cfl_target
        self.enable_transient_logging = enable_transient_logging
        self.transient_log_every = transient_log_every
        self.transient_u_rel_jump_max = transient_u_rel_jump_max
        self.transient_v_rel_jump_max = transient_v_rel_jump_max
        self.transient_p_rel_jump_max = transient_p_rel_jump_max
        self.enable_plots = enable_plots


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

    def _area_ratio_from_mach(self, mach):
        """Isentropic A/A* relation for a perfect gas."""
        gamma = self.config.gamma
        acoef = 2.0 / (gamma + 1.0)
        bcoef = (gamma - 1.0) / (gamma + 1.0)
        ccoef = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        abm = acoef + bcoef * mach * mach
        return (abm**ccoef) / max(mach, 1.0e-12)

    def _solve_mach_from_area_ratio(self, area_ratio, branch):
        """
        Branch-safe bisection solve for Mach from A/A*.
        branch: 'subsonic' or 'supersonic'
        """
        area_target = max(float(area_ratio), 1.0 + 1.0e-10)

        if branch == "subsonic":
            lo = 1.0e-6
            hi = 0.999
        elif branch == "supersonic":
            lo = 1.001
            hi = 8.0
        else:
            raise RuntimeError(f"Unknown branch '{branch}'")

        f_lo = self._area_ratio_from_mach(lo) - area_target
        f_hi = self._area_ratio_from_mach(hi) - area_target
        if f_lo * f_hi > 0.0:
            raise RuntimeError(
                f"Area-Mach bracket failure on {branch} branch: "
                f"A/A*={area_target:.6f}, f(lo)={f_lo:.6e}, f(hi)={f_hi:.6e}"
            )

        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f_mid = self._area_ratio_from_mach(mid) - area_target
            if abs(f_mid) < 1.0e-10 or abs(hi - lo) < 1.0e-10:
                return mid
            if f_lo * f_mid <= 0.0:
                hi = mid
                f_hi = f_mid
            else:
                lo = mid
                f_lo = f_mid
        return 0.5 * (lo + hi)

    def _enforce_onedim_branch(self, l, lt, mach):
        """Enforce subsonic/sonic/supersonic branch by location for baseline mode."""
        if self.config.onedim_mode != "subsonic_sonic_supersonic":
            return mach

        if l < lt:
            return min(mach, 0.999)
        if l == lt:
            return 1.0
        return max(mach, 1.001)

    def solve(self, grid, state):
        if not grid.built or not state.initialized:
            raise RuntimeError("Grid and state must be prepared before solve().")

        # ONEDIM-style initialization from geometry area ratio.
        pt_psf = self.config.pt_in_psia * 144.0
        tt_r = self.config.tt_in_f + 460.0
        gamma = self.config.gamma
        rgas = self.config.rgas
        gam1 = gamma / (gamma - 1.0)
        gam2 = (gamma - 1.0) / 2.0

        state.pt[:] = self.config.pt_in_psia
        state.tt[:] = self.config.tt_in_f
        state.theta[:] = 0.0

        rstar = max(grid.yw[grid.lt], 1.0e-12)
        mach_l = np.zeros(grid.lmax)
        for l in range(grid.lmax):
            area_ratio = max((grid.yw[l] ** 2) / (rstar**2), 1.0 + 1.0e-10)
            if l < grid.lt:
                mach = self._solve_mach_from_area_ratio(area_ratio, branch="subsonic")
            elif l == grid.lt:
                mach = 1.0
            else:
                mach = self._solve_mach_from_area_ratio(area_ratio, branch="supersonic")
            mach = self._enforce_onedim_branch(l, grid.lt, mach)
            mach_l[l] = mach

            dem = 1.0 + gam2 * mach * mach
            p = pt_psf / (dem**gam1)
            t = tt_r / dem
            rho = p / (rgas * t)
            a = np.sqrt(gamma * p / rho)
            u = mach * a

            state.p[l, :, 0] = p
            state.rho[l, :, 0] = rho
            state.u[l, :, 0] = u
            state.v[l, :, 0] = 0.0

        # Keep both levels in sync until predictor/corrector marching is added.
        state.p[:, :, 1] = state.p[:, :, 0]
        state.rho[:, :, 1] = state.rho[:, :, 0]
        state.u[:, :, 1] = state.u[:, :, 0]
        state.v[:, :, 1] = state.v[:, :, 0]

        iterations = 0
        residual_history = [1.0]
        residual_source = "solver"
        converged = False
        transient_cfl_u_max = float("nan")
        transient_dt_min = float("nan")
        transient_dt_max = float("nan")
        if bool(self.config.enable_transient_stub):
            # Pseudo-time march kernel (interior-only for now).
            interior = (slice(1, grid.lmax - 1), slice(1, grid.mmax - 1))
            cfl_u_max = 0.0
            converged = False
            n0 = 0
            n1 = 1
            dt_base = float(self.config.transient_stub_dt_base)
            min_sweeps = max(1, int(self.config.residual_trend_min_points) - 1)
            inlet_u_ref = state.u[0, :, 0].copy()
            inlet_p_ref = state.p[0, :, 0].copy()
            inlet_rho_ref = state.rho[0, :, 0].copy()
            inlet_v_ref = state.v[0, :, 0].copy()
            if self.config.exit_back_pressure_psia is None:
                p_back_psf = float(np.mean(state.p[-1, :, 0]))
            else:
                p_back_psf = float(self.config.exit_back_pressure_psia) * 144.0

            for it in range(1, int(self.config.nmax) + 1):
                # Freeze/copy boundaries and non-marched fields first.
                state.u[:, :, n1] = state.u[:, :, n0]
                state.p[:, :, n1] = state.p[:, :, n0]
                state.rho[:, :, n1] = state.rho[:, :, n0]
                state.v[:, :, n1] = state.v[:, :, n0]

                if grid.lmax >= 3 and grid.mmax >= 3:
                    # Adaptive pseudo-time step based on local convective-acoustic speed.
                    char_max = 0.0
                    for l in range(1, grid.lmax - 1):
                        dx = max(grid.x[l + 1, 1] - grid.x[l - 1, 1], 1.0e-12)
                        for m in range(1, grid.mmax - 1):
                            dy = max(grid.y[l, m + 1] - grid.y[l, m - 1], 1.0e-12)
                            u0 = state.u[l, m, n0]
                            v0 = state.v[l, m, n0]
                            p0 = max(state.p[l, m, n0], 1.0e-10)
                            rho0 = max(state.rho[l, m, n0], 1.0e-12)
                            a0 = np.sqrt(max(self.config.gamma * p0 / rho0, 1.0e-12))
                            char_local = (abs(u0) + a0) / dx + (abs(v0) + a0) / dy
                            char_max = max(char_max, char_local)
                    if char_max > 0.0:
                        dt = min(dt_base, float(self.config.transient_cfl_target) / char_max)
                    else:
                        dt = dt_base
                    transient_dt_min = dt if not np.isfinite(transient_dt_min) else min(transient_dt_min, dt)
                    transient_dt_max = dt if not np.isfinite(transient_dt_max) else max(transient_dt_max, dt)

                    for l in range(1, grid.lmax - 1):
                        dx = max(grid.x[l + 1, 1] - grid.x[l - 1, 1], 1.0e-12)
                        cfl_row = np.max(np.abs(state.u[l, 1 : grid.mmax - 1, n0])) * dt / dx
                        cfl_u_max = max(cfl_u_max, float(cfl_row))
                        for m in range(1, grid.mmax - 1):
                            dy = max(grid.y[l, m + 1] - grid.y[l, m - 1], 1.0e-12)
                            u0 = state.u[l, m, n0]
                            v0 = state.v[l, m, n0]
                            p0 = state.p[l, m, n0]
                            rho0 = max(state.rho[l, m, n0], 1.0e-12)

                            du_dx = (state.u[l + 1, m, n0] - state.u[l - 1, m, n0]) / dx
                            du_dy = (state.u[l, m + 1, n0] - state.u[l, m - 1, n0]) / dy
                            dv_dx = (state.v[l + 1, m, n0] - state.v[l - 1, m, n0]) / dx
                            dv_dy = (state.v[l, m + 1, n0] - state.v[l, m - 1, n0]) / dy
                            dp_dx = (state.p[l + 1, m, n0] - state.p[l - 1, m, n0]) / dx
                            dp_dy = (state.p[l, m + 1, n0] - state.p[l, m - 1, n0]) / dy
                            drho_dx = (state.rho[l + 1, m, n0] - state.rho[l - 1, m, n0]) / dx
                            drho_dy = (state.rho[l, m + 1, n0] - state.rho[l, m - 1, n0]) / dy

                            div_u = du_dx + dv_dy
                            du_dt = -(u0 * du_dx + v0 * du_dy) - dp_dx / rho0
                            dv_dt = -(u0 * dv_dx + v0 * dv_dy) - dp_dy / rho0
                            drho_dt = -(u0 * drho_dx + v0 * drho_dy + rho0 * div_u)
                            dp_dt = -(u0 * dp_dx + v0 * dp_dy + self.config.gamma * p0 * div_u)

                            state.u[l, m, n1] = u0 + dt * du_dt
                            state.v[l, m, n1] = v0 + dt * dv_dt
                            state.rho[l, m, n1] = max(rho0 + dt * drho_dt, 1.0e-10)
                            state.p[l, m, n1] = max(p0 + dt * dp_dt, 1.0e-10)

                    # Boundary-condition updates (phase 2) after interior update.
                    if bool(self.config.enable_inlet_bc):
                        # Characteristic-style subsonic inlet using prescribed totals.
                        l_in = 0
                        l_im = 1
                        gam = self.config.gamma
                        gm1 = gam - 1.0
                        cp = gam * self.config.rgas / gm1
                        for m in range(grid.mmax):
                            ui = state.u[l_im, m, n1]
                            vi = state.v[l_im, m, n1]
                            pi = max(state.p[l_im, m, n1], 1.0e-10)
                            rhoi = max(state.rho[l_im, m, n1], 1.0e-12)
                            ai = np.sqrt(max(gam * pi / rhoi, 1.0e-12))
                            rminus = ui - 2.0 * ai / gm1

                            tt0 = float(state.tt[m]) + 460.0
                            pt0 = float(state.pt[m]) * 144.0

                            # Solve scalar inlet a from 1D characteristic+energy relation.
                            a = max(ai, 1.0)
                            for _ in range(20):
                                u = rminus + 2.0 * a / gm1
                                f = a * a + 0.5 * gm1 * u * u - gam * self.config.rgas * tt0
                                df = 2.0 * a + 2.0 * u
                                if abs(df) < 1.0e-12:
                                    break
                                da = -f / df
                                a = max(a + da, 1.0e-6)
                                if abs(da) < 1.0e-8:
                                    break
                            u_in = rminus + 2.0 * a / gm1
                            t_in = max(a * a / (gam * self.config.rgas), 1.0e-6)
                            mach_in = abs(u_in) / max(a, 1.0e-12)
                            dem = 1.0 + 0.5 * gm1 * mach_in * mach_in
                            p_in = pt0 / (dem ** (gam / gm1))
                            rho_in = p_in / (self.config.rgas * t_in)

                            theta = np.deg2rad(float(state.theta[m]))
                            v_in = u_in * np.tan(theta)

                            state.u[l_in, m, n1] = u_in
                            state.v[l_in, m, n1] = v_in
                            state.p[l_in, m, n1] = max(p_in, 1.0e-10)
                            state.rho[l_in, m, n1] = max(rho_in, 1.0e-12)
                    else:
                        # Preserve original inlet profile when inlet BC is inactive.
                        state.u[0, :, n1] = inlet_u_ref
                        state.p[0, :, n1] = inlet_p_ref
                        state.rho[0, :, n1] = inlet_rho_ref
                        state.v[0, :, n1] = inlet_v_ref

                    if bool(self.config.enable_wall_bc):
                        # Wall tangency BC + normal-gradient compatibility proxies.
                        j_wall = grid.mmax - 1
                        j_im1 = grid.mmax - 2
                        state.v[:, j_wall, n1] = -state.u[:, j_wall, n1] * grid.nxny
                        state.p[:, j_wall, n1] = state.p[:, j_im1, n1]
                        state.rho[:, j_wall, n1] = state.rho[:, j_im1, n1]

                    if bool(self.config.enable_exit_bc):
                        # Regime-based exit BC: supersonic extrapolation, subsonic pressure-corrected.
                        l_exit = grid.lmax - 1
                        l_im1 = grid.lmax - 2
                        l_im2 = grid.lmax - 3 if grid.lmax >= 3 else l_im1
                        gam = self.config.gamma
                        for m in range(grid.mmax):
                            u1 = state.u[l_im1, m, n1]
                            v1 = state.v[l_im1, m, n1]
                            p1 = max(state.p[l_im1, m, n1], 1.0e-10)
                            rho1 = max(state.rho[l_im1, m, n1], 1.0e-12)
                            a1 = np.sqrt(max(gam * p1 / rho1, 1.0e-12))
                            mach1 = np.sqrt(u1 * u1 + v1 * v1) / max(a1, 1.0e-12)
                            if mach1 >= 1.0:
                                state.u[l_exit, m, n1] = state.u[l_im1, m, n1] + (
                                    state.u[l_im1, m, n1] - state.u[l_im2, m, n1]
                                )
                                state.v[l_exit, m, n1] = state.v[l_im1, m, n1] + (
                                    state.v[l_im1, m, n1] - state.v[l_im2, m, n1]
                                )
                                state.p[l_exit, m, n1] = state.p[l_im1, m, n1] + (
                                    state.p[l_im1, m, n1] - state.p[l_im2, m, n1]
                                )
                                state.rho[l_exit, m, n1] = state.rho[l_im1, m, n1] + (
                                    state.rho[l_im1, m, n1] - state.rho[l_im2, m, n1]
                                )
                            else:
                                p_exit = p_back_psf
                                t1 = max(p1 / (rho1 * self.config.rgas), 1.0e-6)
                                rho_exit = p_exit / (self.config.rgas * t1)
                                u_exit = u1 + (p1 - p_exit) / max(rho1 * a1, 1.0e-12)
                                state.u[l_exit, m, n1] = u_exit
                                state.v[l_exit, m, n1] = v1
                                state.p[l_exit, m, n1] = p_exit
                                state.rho[l_exit, m, n1] = rho_exit
                        state.p[l_exit, :, n1] = np.maximum(state.p[l_exit, :, n1], 1.0e-10)
                        state.rho[l_exit, :, n1] = np.maximum(state.rho[l_exit, :, n1], 1.0e-12)

                    u_new = state.u[interior + (n1,)]
                    u_old = state.u[interior + (n0,)]
                    v_new = state.v[interior + (n1,)]
                    v_old = state.v[interior + (n0,)]
                    p_new = state.p[interior + (n1,)]
                    p_old = state.p[interior + (n0,)]
                    rho_new = state.rho[interior + (n1,)]
                    rho_old = state.rho[interior + (n0,)]

                    if not (
                        np.all(np.isfinite(u_new))
                        and np.all(np.isfinite(v_new))
                        and np.all(np.isfinite(p_new))
                        and np.all(np.isfinite(rho_new))
                    ):
                        raise RuntimeError("Transient march produced non-finite interior state.")
                    if np.any(p_new <= 0.0) or np.any(rho_new <= 0.0):
                        raise RuntimeError("Transient march produced non-positive interior pressure/density.")

                    rel_u = np.abs(u_new - u_old) / np.maximum(np.abs(u_old), 1.0e-12)
                    v_scale = max(float(np.max(np.abs(u_old))) * 0.1, 1.0)
                    rel_v = np.abs(v_new - v_old) / np.maximum(np.abs(v_old), v_scale)
                    rel_p = np.abs(p_new - p_old) / np.maximum(np.abs(p_old), 1.0e-12)
                    rel_rho = np.abs(rho_new - rho_old) / np.maximum(np.abs(rho_old), 1.0e-12)
                    res_u = float(np.max(rel_u)) if rel_u.size > 0 else 0.0
                    res_v = float(np.max(rel_v)) if rel_v.size > 0 else 0.0
                    res_p = float(np.max(rel_p)) if rel_p.size > 0 else 0.0
                    res_rho = float(np.max(rel_rho)) if rel_rho.size > 0 else 0.0
                    if res_u > float(self.config.transient_u_rel_jump_max):
                        raise RuntimeError(
                            "Transient march exceeded relative u-jump guardrail: "
                            f"res={res_u:.6g} > max={self.config.transient_u_rel_jump_max:.6g}"
                        )
                    if res_v > float(self.config.transient_v_rel_jump_max):
                        raise RuntimeError(
                            "Transient march exceeded relative v-jump guardrail: "
                            f"res={res_v:.6g} > max={self.config.transient_v_rel_jump_max:.6g}"
                        )
                    if res_p > float(self.config.transient_p_rel_jump_max):
                        raise RuntimeError(
                            "Transient march exceeded relative p-jump guardrail: "
                            f"res={res_p:.6g} > max={self.config.transient_p_rel_jump_max:.6g}"
                        )
                    res = max(res_u, res_v, res_p, res_rho)
                else:
                    res = 0.0

                residual_history.append(max(res, 1.0e-14))
                iterations = it
                if bool(self.config.enable_transient_logging):
                    every = max(1, int(self.config.transient_log_every))
                    if it == 1 or it % every == 0 or it >= min_sweeps:
                        print(
                            "march step "
                            f"{it:04d}: res={res:.6e}, dt={dt:.6e}, cfl_u_max={cfl_u_max:.6e}"
                        )
                if it >= min_sweeps and res <= float(self.config.tconv):
                    converged = True
                    n0 = n1
                    break
                n0, n1 = n1, n0

            # Preserve level 0 as the initialization snapshot and store
            # the marched final state on level 1 for post-processing.
            state.u[:, :, 1] = state.u[:, :, n0]
            state.p[:, :, 1] = state.p[:, :, n0]
            state.rho[:, :, 1] = state.rho[:, :, n0]
            state.v[:, :, 1] = state.v[:, :, n0]

            residual_source = "solver"
            transient_cfl_u_max = float(cfl_u_max)
        elif bool(self.config.use_synthetic_residual_history) and iterations == 0:
            nres = max(2, int(self.config.synthetic_residual_points))
            # Deterministic placeholder sequence to exercise residual-shape gates.
            residual_history = np.geomspace(1.0, 0.25, nres).tolist()
            residual_source = "synthetic"

        return RunResult(
            converged=converged,
            iterations=iterations,
            residual_history=residual_history,
            metadata={
                "case_name": self.config.case_name,
                "state": state,
                "gamma": gamma,
                "residual_source": residual_source,
                "transient_cfl_u_max": transient_cfl_u_max,
                "transient_dt_min": transient_dt_min,
                "transient_dt_max": transient_dt_max,
                "inlet_bc_active": bool(self.config.enable_inlet_bc),
                "wall_bc_active": bool(self.config.enable_wall_bc),
                "exit_bc_active": bool(self.config.enable_exit_bc),
                "mach_inlet": float(mach_l[0]),
                "mach_throat": float(mach_l[grid.lt]),
                "mach_exit": float(mach_l[-1]),
            },
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

        mach_inlet = result.metadata.get("mach_inlet")
        if config.onedim_mode == "subsonic_sonic_supersonic":
            if mach_inlet is None:
                raise RuntimeError("Missing inlet Mach in RunResult metadata.")
            if float(mach_inlet) >= 1.0:
                raise RuntimeError(
                    f"Inlet Mach check failed for subsonic-sonic-supersonic mode: M_inlet={mach_inlet:.6f}"
                )

        state = result.metadata.get("state")
        if state is None:
            raise RuntimeError("RunResult metadata does not include state for verification.")
        residual_source = result.metadata.get("residual_source")
        valid_residual_sources = {"synthetic", "solver"}
        if residual_source not in valid_residual_sources:
            raise RuntimeError(
                "Residual source check failed: expected one of "
                f"{sorted(valid_residual_sources)}, got {residual_source!r}"
            )
        if residual_source == "solver":
            cfl_u_max = result.metadata.get("transient_cfl_u_max")
            if cfl_u_max is None or not np.isfinite(float(cfl_u_max)) or float(cfl_u_max) < 0.0:
                raise RuntimeError(
                    "Transient CFL check failed for solver residual source: "
                    f"transient_cfl_u_max={cfl_u_max!r}"
                )
            dt_min = result.metadata.get("transient_dt_min")
            dt_max = result.metadata.get("transient_dt_max")
            if (
                dt_min is None
                or dt_max is None
                or not np.isfinite(float(dt_min))
                or not np.isfinite(float(dt_max))
                or float(dt_min) <= 0.0
                or float(dt_max) <= 0.0
                or float(dt_min) > float(dt_max)
            ):
                raise RuntimeError(
                    "Transient dt range check failed for solver residual source: "
                    f"transient_dt_min={dt_min!r}, transient_dt_max={dt_max!r}"
                )

        residuals = np.asarray(result.residual_history, dtype=float)
        if residuals.size == 0:
            raise RuntimeError("Residual sanity check failed: residual history is empty.")
        if not np.all(np.isfinite(residuals)):
            bad_idx = int(np.where(~np.isfinite(residuals))[0][0])
            raise RuntimeError(
                "Residual sanity check failed: residual history contains non-finite values. "
                f"First bad index={bad_idx}, value={residuals[bad_idx]}"
            )
        if np.any(residuals < 0.0):
            bad_idx = int(np.where(residuals < 0.0)[0][0])
            raise RuntimeError(
                "Residual sanity check failed: residual history contains negative values. "
                f"First bad index={bad_idx}, value={residuals[bad_idx]}"
            )
        residual_shape_check = {
            "status": "not_applicable_yet",
            "n_points": int(residuals.size),
            "start": float(residuals[0]),
            "end": float(residuals[-1]),
            "end_over_start": float("nan"),
            "nonincreasing_frac": float("nan"),
            "min_points_required": int(config.residual_trend_min_points),
            "end_over_start_max": float(config.residual_end_over_start_max),
            "nonincreasing_frac_min": float(config.residual_nonincreasing_frac_min),
            "passed": True,
        }
        if 1 < residuals.size < int(config.residual_trend_min_points):
            residual_shape_check["status"] = "insufficient_points"
        elif residuals.size >= int(config.residual_trend_min_points):
            diffs = np.diff(residuals)
            nonincreasing_frac = float(np.mean(diffs <= 0.0))
            end_over_start = (
                float(residuals[-1] / residuals[0]) if abs(float(residuals[0])) > 1.0e-12 else float("inf")
            )
            # Convergence-shape criteria for marching runs.
            passed_end_ratio = bool(end_over_start <= float(config.residual_end_over_start_max))
            passed_nonincreasing = bool(
                nonincreasing_frac >= float(config.residual_nonincreasing_frac_min)
            )
            passed = bool(passed_end_ratio and passed_nonincreasing)
            residual_shape_check = {
                "status": "evaluated",
                "n_points": int(residuals.size),
                "start": float(residuals[0]),
                "end": float(residuals[-1]),
                "end_over_start": end_over_start,
                "nonincreasing_frac": nonincreasing_frac,
                "min_points_required": int(config.residual_trend_min_points),
                "end_over_start_max": float(config.residual_end_over_start_max),
                "nonincreasing_frac_min": float(config.residual_nonincreasing_frac_min),
                "passed_end_ratio": passed_end_ratio,
                "passed_nonincreasing": passed_nonincreasing,
                "passed": passed,
            }
            if not passed:
                raise RuntimeError(
                    "Residual shape check failed: "
                    f"end_over_start={end_over_start:.6g} "
                    f"(max={config.residual_end_over_start_max:.6g}), "
                    f"nonincreasing_frac={nonincreasing_frac:.6g} "
                    f"(min={config.residual_nonincreasing_frac_min:.6g})."
                )

        data = np.genfromtxt(config.baseline_csv_path, delimiter=",", names=True)
        rows = np.atleast_1d(data)
        if rows.size == 0:
            raise RuntimeError(f"No rows found in baseline CSV: {config.baseline_csv_path}")

        required_cols = ["L", "M", "U", "V", "P", "RHO", "Mach"]
        for col in required_cols:
            if col not in data.dtype.names:
                raise RuntimeError(
                    f"Baseline CSV missing required column '{col}'. "
                    f"Found columns: {data.dtype.names}"
                )
        if "X" in data.dtype.names and "Y" in data.dtype.names:
            x_csv = np.asarray(rows["X"], dtype=float)
            y_csv = np.asarray(rows["Y"], dtype=float)
            finite_xy = np.isfinite(x_csv) & np.isfinite(y_csv)
            if np.any(finite_xy):
                x_use = x_csv[finite_xy]
                y_use = y_csv[finite_xy]
                y_wall_csv = np.interp(x_use, grid.xw, grid.yw)
                tol = float(config.csv_domain_tol)
                outside = (y_use < -tol) | (y_use > y_wall_csv + tol)
                if np.any(outside):
                    n_bad = int(np.sum(outside))
                    i0 = int(np.where(outside)[0][0])
                    raise RuntimeError(
                        "CSV domain check failed: "
                        f"{n_bad} points are outside nozzle bounds (tol={tol}). "
                        f"First bad point: x={x_use[i0]:.6f}, y={y_use[i0]:.6f}, "
                        f"y_wall={y_wall_csv[i0]:.6f}"
                    )

        gamma = float(result.metadata.get("gamma", config.gamma))
        nlevels = state.u.shape[2]
        level = int(config.verify_level)
        if level < 0 or level >= nlevels:
            raise RuntimeError(
                f"verify_level={level} is out of range for available levels [0, {nlevels - 1}]"
            )

        u_model = state.u[:, :, level]
        v_model = state.v[:, :, level]
        p_model_psia = state.p[:, :, level] / 144.0
        rho_model = state.rho[:, :, level]
        a_model = np.sqrt(np.maximum(1.0e-12, gamma * state.p[:, :, level] / rho_model))
        mach_model = np.sqrt(u_model * u_model + v_model * v_model) / a_model

        u_ref = np.full((grid.lmax, grid.mmax), np.nan)
        v_ref = np.full((grid.lmax, grid.mmax), np.nan)
        p_ref_psia = np.full((grid.lmax, grid.mmax), np.nan)
        mach_ref = np.full((grid.lmax, grid.mmax), np.nan)
        wall_ref_mask = np.zeros((grid.lmax, grid.mmax), dtype=bool)
        m_wall_ref = int(np.nanmax(np.asarray(rows["M"], dtype=float)))

        for row in rows:
            l_idx = int(round(float(row["L"]))) - 1
            m_idx = int(round(float(row["M"]))) - 1
            if 0 <= l_idx < grid.lmax and 0 <= m_idx < grid.mmax:
                u_ref[l_idx, m_idx] = float(row["U"])
                v_ref[l_idx, m_idx] = float(row["V"])
                p_ref_psia[l_idx, m_idx] = float(row["P"])
                mach_ref[l_idx, m_idx] = float(row["Mach"])
                if int(round(float(row["M"]))) == m_wall_ref:
                    wall_ref_mask[l_idx, m_idx] = True

        ref_mask = np.isfinite(mach_ref)
        if not np.any(ref_mask):
            raise RuntimeError("No valid reference points mapped from baseline CSV.")

        def _mae(model_vals, ref_vals, mask):
            vals = np.abs(model_vals[mask] - ref_vals[mask])
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                return float("nan")
            return float(np.mean(vals))

        def _mape(model_vals, ref_vals, mask, floor=1.0e-8):
            mv = model_vals[mask]
            rv = ref_vals[mask]
            finite = np.isfinite(mv) & np.isfinite(rv)
            if np.sum(finite) == 0:
                return float("nan")
            denom = np.maximum(np.abs(rv[finite]), floor)
            return float(np.mean(np.abs((mv[finite] - rv[finite]) / denom)) * 100.0)

        pt_psia = float(np.mean(state.pt)) if np.any(state.pt > 0.0) else float("nan")
        wall_mask = wall_ref_mask & np.isfinite(p_ref_psia)
        p_ratio_mae = float("nan")
        p_ratio_mape = float("nan")
        p_ratio_inlet_mae = float("nan")
        p_ratio_inlet_mape = float("nan")
        if np.any(wall_mask) and np.isfinite(pt_psia) and pt_psia > 0.0:
            p_ratio_model = p_model_psia / pt_psia
            p_ratio_ref = p_ref_psia / pt_psia
            p_ratio_mae = _mae(p_ratio_model, p_ratio_ref, wall_mask)
            p_ratio_mape = _mape(p_ratio_model, p_ratio_ref, wall_mask)

        # Diagnostic-only alternate convention:
        # normalize wall static pressure by inlet wall static pressure.
        j_wall = int(m_wall_ref) - 1
        if 0 <= j_wall < grid.mmax:
            wall_ref_line = p_ref_psia[:, j_wall]
            wall_model_line = p_model_psia[:, j_wall]
            valid_line = np.isfinite(wall_ref_line) & np.isfinite(wall_model_line)
            if np.any(valid_line):
                l_inlet = int(np.where(valid_line)[0][0])
                p_ref_inlet = float(wall_ref_line[l_inlet])
                p_model_inlet = float(wall_model_line[l_inlet])
                if abs(p_ref_inlet) > 1.0e-12 and abs(p_model_inlet) > 1.0e-12:
                    p_ratio_inlet_ref = p_ref_psia / p_ref_inlet
                    p_ratio_inlet_model = p_model_psia / p_model_inlet
                    p_ratio_inlet_mae = _mae(p_ratio_inlet_model, p_ratio_inlet_ref, wall_mask)
                    p_ratio_inlet_mape = _mape(p_ratio_inlet_model, p_ratio_inlet_ref, wall_mask)

        metrics = {
            "count": int(np.sum(ref_mask)),
            "u_mae_fps": _mae(u_model, u_ref, ref_mask),
            "u_mape_pct": _mape(u_model, u_ref, ref_mask),
            "v_mae_fps": _mae(v_model, v_ref, ref_mask),
            "v_mape_pct": _mape(v_model, v_ref, ref_mask),
            "p_mae_psia": _mae(p_model_psia, p_ref_psia, ref_mask),
            "p_mape_pct": _mape(p_model_psia, p_ref_psia, ref_mask),
            "mach_mae": _mae(mach_model, mach_ref, ref_mask),
            "mach_mape_pct": _mape(mach_model, mach_ref, ref_mask),
            "wall_p_ratio_mae": p_ratio_mae,
            "wall_p_ratio_mape_pct": p_ratio_mape,
            "wall_p_ratio_inlet_norm_mae": p_ratio_inlet_mae,
            "wall_p_ratio_inlet_norm_mape_pct": p_ratio_inlet_mape,
        }

        # Reference-method check for wall pressure ratio convention.
        # Use a conservative threshold so we only recommend switching when
        # there is a clear improvement against the current baseline reference.
        method_check = {
            "canonical": "P/Pt",
            "alternate": "P/P_wall,inlet",
            "canonical_mae": float(metrics["wall_p_ratio_mae"]),
            "alternate_mae": float(metrics["wall_p_ratio_inlet_norm_mae"]),
            "improvement_frac": float("nan"),
            "recommended": "keep_canonical",
        }
        can_mae = method_check["canonical_mae"]
        alt_mae = method_check["alternate_mae"]
        if np.isfinite(can_mae) and np.isfinite(alt_mae) and can_mae > 0.0:
            improvement = (can_mae - alt_mae) / can_mae
            method_check["improvement_frac"] = float(improvement)
            # 5% threshold to avoid churn from noise-level differences.
            if improvement >= 0.05:
                method_check["recommended"] = "switch_to_alternate"

        for key, value in metrics.items():
            if key == "count":
                continue
            if not np.isfinite(value):
                raise RuntimeError(f"Verification metric '{key}' is not finite.")

        # Modern regression gates: fail if key errors regress beyond baseline + margin.
        margin = float(config.regression_margin_frac)
        gates = {
            "mach_mae": {
                "baseline": float(config.baseline_mach_mae),
                "current": float(metrics["mach_mae"]),
            },
            "p_mae_psia": {
                "baseline": float(config.baseline_p_mae_psia),
                "current": float(metrics["p_mae_psia"]),
            },
            "wall_p_ratio_mae": {
                "baseline": float(config.baseline_wall_p_ratio_mae),
                "current": float(metrics["wall_p_ratio_mae"]),
            },
        }
        for gate_name, gate in gates.items():
            baseline = gate["baseline"]
            current = gate["current"]
            allowed = baseline * (1.0 + margin)
            gate["allowed"] = allowed
            gate["passed"] = bool(current <= allowed)
            if not gate["passed"]:
                raise RuntimeError(
                    f"Regression gate failed for {gate_name}: "
                    f"current={current:.6g} > allowed={allowed:.6g} "
                    f"(baseline={baseline:.6g}, margin={margin:.1%})"
                )

        lt_value = getattr(grid, "lt", None)
        return {
            "checked": True,
            "xt": float(grid.xt),
            "xt_expected": float(config.xt_expected),
            "xt_error": float(xt_err),
            "xt_tol": float(config.xt_tol),
            "verify_level": level,
            "lt_report": lt_value,
            "count": metrics["count"],
            "u_mae_fps": metrics["u_mae_fps"],
            "u_mape_pct": metrics["u_mape_pct"],
            "v_mae_fps": metrics["v_mae_fps"],
            "v_mape_pct": metrics["v_mape_pct"],
            "p_mae_psia": metrics["p_mae_psia"],
            "p_mape_pct": metrics["p_mape_pct"],
            "mach_mae": metrics["mach_mae"],
            "mach_mape_pct": metrics["mach_mape_pct"],
            "wall_p_ratio_mae": metrics["wall_p_ratio_mae"],
            "wall_p_ratio_mape_pct": metrics["wall_p_ratio_mape_pct"],
            "wall_p_ratio_inlet_norm_mae": metrics["wall_p_ratio_inlet_norm_mae"],
            "wall_p_ratio_inlet_norm_mape_pct": metrics["wall_p_ratio_inlet_norm_mape_pct"],
            "wall_ratio_method_check": method_check,
            "residual_points": int(residuals.size),
            "residual_start": float(residuals[0]),
            "residual_end": float(residuals[-1]),
            "residual_min": float(np.min(residuals)),
            "residual_shape_check": residual_shape_check,
            "residual_source": residual_source,
            "regression_margin_frac": margin,
            "regression_gates": gates,
        }


class Plots:
    def __init__(self, config, grid, result: RunResult, verify_report: dict):
        self.config = config
        self.grid = grid
        self.result = result
        self.verify_report = verify_report

    def _plot_level(self, state):
        nlevels = state.u.shape[2]
        level = int(getattr(self.config, "verify_level", 0))
        return min(max(level, 0), nlevels - 1)

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
        import matplotlib.pyplot as plt

        state = self.result.metadata.get("state")
        if state is None:
            raise RuntimeError("RunResult metadata does not include state for solution plotting.")

        os.makedirs("artifacts/verification", exist_ok=True)
        output_path = "artifacts/verification/cd_baseline_solution.png"
        level = self._plot_level(state)
        output_path_2d = f"artifacts/verification/cd_baseline_onedim_2d_level{level}.png"

        p = state.p[:, :, level] / 144.0  # psia
        u = state.u[:, :, level]
        x = self.grid.x[:, 0]
        j_center = 0
        j_wall = self.grid.mmax - 1

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=140)
        ax0, ax1 = axes

        ax0.plot(x, p[:, j_center], lw=2.0, color="#1f77b4", label="centerline")
        ax0.plot(x, p[:, j_wall], lw=2.0, color="#d62728", label="wall")
        ax0.set_title(f"Static Pressure (Level {level})")
        ax0.set_xlabel("x")
        ax0.set_ylabel("Pressure (psia)")
        ax0.grid(alpha=0.2)
        ax0.legend(loc="best")

        ax1.plot(x, u[:, j_center], lw=2.0, color="#2ca02c", label="centerline")
        ax1.plot(x, u[:, j_wall], lw=2.0, color="#9467bd", label="wall")
        ax1.set_title(f"Axial Velocity (Level {level})")
        ax1.set_xlabel("x")
        ax1.set_ylabel("u")
        ax1.grid(alpha=0.2)
        ax1.legend(loc="best")

        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved solution plot: {output_path}")

        # 2D view with grid overlay from selected marched level.
        mach = np.sqrt(state.u[:, :, level] ** 2 + state.v[:, :, level] ** 2) / np.sqrt(
            np.maximum(1.0e-12, self.result.metadata.get("gamma", 1.4) * state.p[:, :, level] / state.rho[:, :, level])
        )

        fig2, ax2 = plt.subplots(figsize=(10, 4.8), dpi=140)
        cf = ax2.contourf(self.grid.x, self.grid.y, mach, levels=20, cmap="viridis")
        for l in range(self.grid.lmax):
            ax2.plot(self.grid.x[l, :], self.grid.y[l, :], color="white", lw=0.35, alpha=0.6)
        for m in range(self.grid.mmax):
            ax2.plot(self.grid.x[:, m], self.grid.y[:, m], color="white", lw=0.35, alpha=0.6)
        ax2.plot(self.grid.xw, self.grid.yw, color="k", lw=1.4, label="wall")
        ax2.plot(self.grid.xcb, self.grid.ycb, color="k", lw=1.2, ls="--", label="centerline")
        ax2.set_title(f"Mach Field with Structured Grid Overlay (Level {level})")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_aspect("equal", adjustable="box")
        ax2.legend(loc="best")
        cbar = fig2.colorbar(cf, ax=ax2)
        cbar.set_label("Mach")
        fig2.tight_layout()
        fig2.savefig(output_path_2d, bbox_inches="tight")
        plt.show()
        plt.close(fig2)
        print(f"Saved 2D marched solution plot: {output_path_2d}")

    def plot_residuals(self) -> None:
        import matplotlib.pyplot as plt

        os.makedirs("artifacts/verification", exist_ok=True)
        output_path = "artifacts/verification/cd_baseline_residual.png"

        res = np.asarray(self.result.residual_history, dtype=float)
        it = np.arange(1, len(res) + 1)

        fig, ax = plt.subplots(figsize=(7, 4), dpi=140)
        if len(res) > 0:
            ax.plot(it, res, marker="o", lw=1.8, color="#1f77b4")
        ax.set_title("Residual History")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Residual")
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved residual plot: {output_path}")

    def plot_fig2_style(self) -> None:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri

        state = self.result.metadata.get("state")
        if state is None:
            raise RuntimeError("RunResult metadata does not include state for Fig. 2 style plot.")

        os.makedirs("artifacts/verification", exist_ok=True)
        output_path = "artifacts/verification/cd_baseline_mach_contours.png"

        gamma = self.result.metadata.get("gamma", 1.4)
        mach = np.sqrt(state.u[:, :, 0] ** 2 + state.v[:, :, 0] ** 2) / np.sqrt(
            np.maximum(1.0e-12, gamma * state.p[:, :, 0] / state.rho[:, :, 0])
        )

        fig, ax_top = plt.subplots(figsize=(8.6, 4.8), dpi=150)

        # Top panel: Mach contours and nozzle wall.
        levels = [0.6, 1.0, 1.6]
        cs = ax_top.contour(
            self.grid.x,
            self.grid.y,
            mach,
            levels=levels,
            colors="k",
            linewidths=[1.4, 1.8, 1.4],
        )
        ax_top.clabel(
            cs,
            cs.levels,
            inline=True,
            fmt={0.6: "M = 0.6", 1.0: "1.0", 1.6: "1.6"},
            fontsize=9,
        )

        # Overlay verification Mach contours from baseline CSV, if available.
        # Interpolate scattered CSV data to a dense structured grid first so
        # reference contour lines are smooth and confined to the nozzle domain.
        if os.path.exists(self.config.baseline_csv_path):
            data = np.genfromtxt(self.config.baseline_csv_path, delimiter=",", names=True)
            rows = np.atleast_1d(data)
            if rows.size > 0 and all(c in data.dtype.names for c in ("X", "Y", "Mach")):
                x_ref = np.asarray(rows["X"], dtype=float)
                y_ref = np.asarray(rows["Y"], dtype=float)
                mach_ref = np.asarray(rows["Mach"], dtype=float)
                finite = np.isfinite(x_ref) & np.isfinite(y_ref) & np.isfinite(mach_ref)
                if np.sum(finite) >= 6:
                    tri = mtri.Triangulation(x_ref[finite], y_ref[finite])
                    interp = mtri.LinearTriInterpolator(tri, mach_ref[finite])

                    # Dense sampling grid for smooth contours.
                    xg = np.linspace(np.min(x_ref[finite]), np.max(x_ref[finite]), 220)
                    yg = np.linspace(0.0, np.max(y_ref[finite]), 180)
                    Xg, Yg = np.meshgrid(xg, yg)
                    Mg = interp(Xg, Yg)

                    # Mask outside nozzle wall/centerline bounds.
                    ywall = np.interp(xg, self.grid.xw, self.grid.yw)
                    Ywall = np.tile(ywall, (yg.size, 1))
                    outside = (Yg < 0.0) | (Yg > Ywall)
                    Mg = np.ma.masked_where(outside, Mg)

                    cs_ref = ax_top.contour(
                        Xg,
                        Yg,
                        Mg,
                        levels=levels,
                        colors="k",
                        linewidths=[1.0, 1.3, 1.0],
                        linestyles="--",
                    )
                    if len(cs_ref.allsegs) > 0:
                        ax_top.clabel(
                            cs_ref,
                            cs_ref.levels,
                            inline=True,
                            fmt={0.6: "ref 0.6", 1.0: "ref 1.0", 1.6: "ref 1.6"},
                            fontsize=8,
                        )
                        ax_top.plot([], [], color="k", ls="--", lw=1.2, label="CSV Mach contours")

        ax_top.plot(self.grid.xw, self.grid.yw, color="k", lw=2.0)
        ax_top.plot(self.grid.xcb, self.grid.ycb, color="k", lw=1.3)
        ax_top.set_ylabel("Radius")
        ax_top.set_xlabel("Axial Distance")
        ax_top.set_ylim(bottom=0.0)
        ax_top.grid(alpha=0.15)
        ax_top.set_title("Mach Number Contours (ONEDIM Field)")
        ax_top.legend(loc="best")

        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved Mach contour plot: {output_path}")

    def plot_wall_pressure_comparison(self) -> None:
        import matplotlib.pyplot as plt

        state = self.result.metadata.get("state")
        if state is None:
            raise RuntimeError("RunResult metadata does not include state for wall pressure plotting.")
        if not os.path.exists(self.config.baseline_csv_path):
            raise RuntimeError(f"Baseline CSV not found: {self.config.baseline_csv_path}")

        os.makedirs("artifacts/verification", exist_ok=True)
        output_path = "artifacts/verification/cd_baseline_wall_pressure_ratio.png"

        # Model wall pressure ratio from current state.
        j_wall = self.grid.mmax - 1
        x_model = self.grid.x[:, j_wall]
        pt_psf = float(np.mean(state.pt)) * 144.0 if np.any(state.pt > 0.0) else 1.0
        p_ratio_model = state.p[:, j_wall, 0] / max(pt_psf, 1.0e-12)

        # Reference wall pressure ratio from baseline CSV.
        data = np.genfromtxt(self.config.baseline_csv_path, delimiter=",", names=True)
        rows = np.atleast_1d(data)
        if rows.size == 0:
            raise RuntimeError(f"No rows found in baseline CSV: {self.config.baseline_csv_path}")
        if "M" not in data.dtype.names or "X" not in data.dtype.names or "P" not in data.dtype.names:
            raise RuntimeError(f"CSV missing one of required columns: M, X, P. Found {data.dtype.names}")

        m_wall = int(np.nanmax(np.asarray(rows["M"], dtype=float)))
        wall_rows = [r for r in rows if int(round(float(r["M"]))) == m_wall]
        if len(wall_rows) == 0:
            raise RuntimeError("No wall rows found in baseline CSV.")
        wall_rows = sorted(wall_rows, key=lambda r: float(r["X"]))
        x_ref = np.array([float(r["X"]) for r in wall_rows], dtype=float)
        p_ref_psia = np.array([float(r["P"]) for r in wall_rows], dtype=float)
        pt_psia = float(np.mean(state.pt)) if np.any(state.pt > 0.0) else 1.0
        p_ratio_ref = p_ref_psia / max(pt_psia, 1.0e-12)
        p_model_psia_wall = state.p[:, j_wall, 0] / 144.0
        p_ratio_model_inlet = p_model_psia_wall / max(float(p_model_psia_wall[0]), 1.0e-12)
        p_ratio_ref_inlet = p_ref_psia / max(float(p_ref_psia[0]), 1.0e-12)

        p_ratio_ref_on_model = np.interp(x_model, x_ref, p_ratio_ref)
        p_ratio_ref_inlet_on_model = np.interp(x_model, x_ref, p_ratio_ref_inlet)
        mae_pt = float(np.mean(np.abs(p_ratio_model - p_ratio_ref_on_model)))
        mae_inlet = float(np.mean(np.abs(p_ratio_model_inlet - p_ratio_ref_inlet_on_model)))

        fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), dpi=150, sharex=True)
        ax0, ax1 = axes

        ax0.plot(x_model, p_ratio_model, color="k", lw=2.0, label="Model")
        ax0.scatter(x_ref, p_ratio_ref, s=20, facecolors="none", edgecolors="k", label="CSV")
        ax0.set_title(f"P/Pt (MAE={mae_pt:.4f})")
        ax0.set_xlabel("Axial Distance")
        ax0.set_ylabel("Wall Pressure Ratio")
        ax0.set_ylim(bottom=0.0)
        ax0.grid(alpha=0.15)
        ax0.legend(loc="best")

        ax1.plot(x_model, p_ratio_model_inlet, color="k", lw=2.0, label="Model")
        ax1.scatter(x_ref, p_ratio_ref_inlet, s=20, facecolors="none", edgecolors="k", label="CSV")
        ax1.set_title(f"P/P_wall,inlet (MAE={mae_inlet:.4f})")
        ax1.set_xlabel("Axial Distance")
        ax1.set_ylabel("Wall Pressure Ratio")
        ax1.set_ylim(bottom=0.0)
        ax1.grid(alpha=0.15)
        ax1.legend(loc="best")
        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved wall pressure comparison plot: {output_path}")





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
    if config.enable_plots:
        plots = Plots(config, grid, result, verify_report)
        plots.plot_grid()
        plots.plot_solution()
        plots.plot_residuals()
        plots.plot_fig2_style()
        plots.plot_wall_pressure_comparison()
    else:
        print("plots: skipped (enable_plots=False)")

    # 7) Print concise run summary.
    print("NAP run summary")
    print(f"case={config.case_name}")
    print(f"grid={config.lmax}x{config.mmax}, converged={result.converged}, iterations={result.iterations}")
    print(
        "verify: "
        f"XT={verify_report['xt']:.3f} (target {verify_report['xt_expected']:.3f}, err {verify_report['xt_error']:.3e}), "
        f"verify_level={verify_report['verify_level']}, "
        f"LT={verify_report['lt_report']}, "
        f"count={verify_report['count']}, "
        f"U_MAE={verify_report['u_mae_fps']}, "
        f"V_MAE={verify_report['v_mae_fps']}, "
        f"P_MAE={verify_report['p_mae_psia']}, "
        f"Mach_MAE={verify_report['mach_mae']}"
    )
    print(
        "verify wall ratio: "
        f"P/Pt_MAE={verify_report['wall_p_ratio_mae']}, "
        f"P/Pt_MAPE%={verify_report['wall_p_ratio_mape_pct']}, "
        f"P/Pwall,inlet_MAE={verify_report['wall_p_ratio_inlet_norm_mae']}, "
        f"P/Pwall,inlet_MAPE%={verify_report['wall_p_ratio_inlet_norm_mape_pct']}"
    )
    method_check = verify_report["wall_ratio_method_check"]
    print(
        "wall ratio method check: "
        f"canonical={method_check['canonical']}, "
        f"alternate={method_check['alternate']}, "
        f"improvement_frac={method_check['improvement_frac']}, "
        f"recommendation={method_check['recommended']}"
    )
    print(
        "residual sanity: "
        f"points={verify_report['residual_points']}, "
        f"start={verify_report['residual_start']}, "
        f"end={verify_report['residual_end']}, "
        f"min={verify_report['residual_min']}"
    )
    print(
        "time march: "
        f"enabled={config.enable_transient_stub}, "
        f"source={verify_report['residual_source']}, "
        f"converged={result.converged}, "
        f"iterations={result.iterations}, "
        f"tconv={config.tconv}"
    )
    print(f"residual source: {verify_report['residual_source']}")
    print(f"transient cfl_u max: {result.metadata.get('transient_cfl_u_max')}")
    print(
        "transient dt range: "
        f"min={result.metadata.get('transient_dt_min')}, "
        f"max={result.metadata.get('transient_dt_max')}"
    )
    residual_shape = verify_report["residual_shape_check"]
    print(
        "residual shape: "
        f"status={residual_shape['status']}, "
        f"n_points={residual_shape['n_points']}, "
        f"end_over_start={residual_shape['end_over_start']}, "
        f"nonincreasing_frac={residual_shape['nonincreasing_frac']}, "
        f"min_points_required={residual_shape['min_points_required']}, "
        f"end_over_start_max={residual_shape['end_over_start_max']}, "
        f"nonincreasing_frac_min={residual_shape['nonincreasing_frac_min']}, "
        f"pass={residual_shape['passed']}"
    )
    gates = verify_report["regression_gates"]
    print("regression gates:")
    print(
        "  mach_mae: "
        f"baseline={gates['mach_mae']['baseline']:.6g}, "
        f"current={gates['mach_mae']['current']:.6g}, "
        f"allowed={gates['mach_mae']['allowed']:.6g}, "
        f"pass={gates['mach_mae']['passed']}"
    )
    print(
        "  p_mae_psia: "
        f"baseline={gates['p_mae_psia']['baseline']:.6g}, "
        f"current={gates['p_mae_psia']['current']:.6g}, "
        f"allowed={gates['p_mae_psia']['allowed']:.6g}, "
        f"pass={gates['p_mae_psia']['passed']}"
    )
    print(
        "  wall_p_ratio_mae: "
        f"baseline={gates['wall_p_ratio_mae']['baseline']:.6g}, "
        f"current={gates['wall_p_ratio_mae']['current']:.6g}, "
        f"allowed={gates['wall_p_ratio_mae']['allowed']:.6g}, "
        f"pass={gates['wall_p_ratio_mae']['passed']}"
    )
    print(
        "onedim: "
        f"Mach_inlet={result.metadata['mach_inlet']:.3f}, "
        f"Mach_throat={result.metadata['mach_throat']:.3f}, "
        f"Mach_exit={result.metadata['mach_exit']:.3f}"
    )
