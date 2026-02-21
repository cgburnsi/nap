from dataclasses import dataclass
import argparse
import os

import numpy as np
from ver1.data_structs import History, Solution

class StructuredGrid:
    """Single in-memory structured grid model with geometry hard-coded into build."""

    def __init__(
        self,
        lmax=21,
        mmax=8,
        xi=0.31,
        ri=2.50,
        rt=0.80,
        xe=4.05,
        rci=0.80,
        rct=0.50,
        ani_deg=44.88,
        ane_deg=15.0,
        ndim=1,
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
        self.ndim = ndim

    def build(self):
        """Populate geometry, coordinates, mapping metrics, boundaries, and state arrays."""
        ani = np.deg2rad(self.ani_deg)
        ane = np.deg2rad(self.ane_deg)
        xt, re, lt, yw, nxny = geometry(
            self.lmax,
            self.xi,
            self.ri,
            self.rt,
            self.xe,
            self.rci,
            self.rct,
            ani,
            ane,
        )
        self.xt = xt
        self.re = re
        self.lt = lt
        self.yw = yw
        self.nxny = nxny
        self.xw = np.linspace(self.xi, self.xe, self.lmax)

        self.xcb = self.xw.copy()
        self.ycb = np.zeros(self.lmax)
        self.nxnycb = np.zeros(self.lmax)

        self.x, self.y = _build_structured_coordinates(
            self.xw,
            self.yw,
            self.ycb,
            self.mmax,
            self.ndim,
        )
        self.boundary_lm = _build_boundary_index_sets(self.lmax, self.mmax, self.lt)
        self.jac, self.xi_x, self.xi_y, self.eta_x, self.eta_y = _compute_mapping_metrics(
            self.x, self.y
        )

        shape_3d = (self.lmax, self.mmax, 2)
        self.u = np.zeros(shape_3d, dtype=float, order="F")
        self.v = np.zeros(shape_3d, dtype=float, order="F")
        self.p = np.zeros(shape_3d, dtype=float, order="F")
        self.ro = np.zeros(shape_3d, dtype=float, order="F")

        self.pt = np.zeros(self.mmax)
        self.tt = np.zeros(self.mmax)
        self.theta = np.zeros(self.mmax)
        return self


def geometry(lmax, xi, ri, rt, xe, rci, rct, ani, ane):
    """Generate NAP-style nozzle wall arrays YW(L) and NXNY(L)."""
    yw, nxny = np.zeros(lmax), np.zeros(lmax)
    xt, re, lt = 0.0, 0.0, 0

    xtan = xi + rci * np.sin(ani)
    rtan = ri + rci * (np.cos(ani) - 1.0)
    rt1 = rt - rct * (np.cos(ani) - 1.0)
    xt1 = xtan + (rtan - rt1) / np.tan(ani)
    if xt1 < xtan:
        xt1, rt1 = xtan, rtan
    xt = xt1 + rct * np.sin(ani)
    xt2 = xt + rct * np.sin(ane)

    rt2 = rt + rct * (1.0 - np.cos(ane))
    re = rt2 + (xe - xt2) * np.tan(ane)
    dx = (xe - xi) / (lmax - 1)

    for l in range(lmax):
        x = xi + dx * l
        if xi <= x <= xtan:
            yw[l] = ri + rci * (np.cos(np.arcsin((x - xi) / rci)) - 1.0)
            nxny[l] = (x - xi) / (yw[l] - ri + rci)
        elif xtan < x <= xt1:
            yw[l] = rt1 + (xt1 - x) * np.tan(ani)
            nxny[l] = np.tan(ani)
        elif xt1 < x <= xt:
            yw[l] = rt + rct * (1.0 - np.cos(np.arcsin((xt - x) / rct)))
            nxny[l] = (xt - x) / (rct + rt - yw[l])
        elif xt < x <= xt2:
            yw[l] = rt + rct * (1.0 - np.cos(np.arcsin((x - xt) / rct)))
            nxny[l] = (xt - x) / (rct + rt - yw[l])
        elif xt2 < x <= xe:
            yw[l] = rt2 + (x - xt2) * np.tan(ane)
            nxny[l] = -np.tan(ane)

        if l > 0 and yw[l] < yw[lt]:
            lt = l

    return xt, re, lt, yw, nxny


@dataclass
class NapGrid:
    """Structured grid and NAP-style arrays (Python 0-based indexing)."""

    lmax: int
    mmax: int
    ndim: int
    xi: float
    xe: float
    xt: float
    ri: float
    rt: float
    re: float
    lt: int
    xw: np.ndarray
    yw: np.ndarray
    nxny: np.ndarray
    xcb: np.ndarray
    ycb: np.ndarray
    nxnycb: np.ndarray
    x: np.ndarray
    y: np.ndarray
    u: np.ndarray
    v: np.ndarray
    p: np.ndarray
    ro: np.ndarray
    pt: np.ndarray
    tt: np.ndarray
    theta: np.ndarray
    boundary_lm: dict
    jac: np.ndarray
    xi_x: np.ndarray
    xi_y: np.ndarray
    eta_x: np.ndarray
    eta_y: np.ndarray


def _build_structured_coordinates(xw, yw, ycb, mmax, ndim):
    """Create X(L,M), Y(L,M) from wall/centerbody contours."""
    lmax = xw.shape[0]
    x = np.repeat(xw[:, None], mmax, axis=1)
    y = np.zeros((lmax, mmax), dtype=float, order="F")

    eta = np.linspace(0.0, 1.0, mmax)
    for l in range(lmax):
        if ndim == 0:
            y[l, :] = ycb[l] + eta * (yw[l] - ycb[l])
        else:
            y2 = ycb[l] ** 2 + eta * (yw[l] ** 2 - ycb[l] ** 2)
            y[l, :] = np.sqrt(np.maximum(y2, 0.0))
    return x, y


def _build_boundary_index_sets(lmax, mmax, lt):
    """Build NAP-style boundary/index sets in (L, M) 0-based indexing."""
    boundaries = {
        "inlet": np.column_stack((np.zeros(mmax, dtype=int), np.arange(mmax, dtype=int))),
        "exit": np.column_stack(
            ((lmax - 1) * np.ones(mmax, dtype=int), np.arange(mmax, dtype=int))
        ),
        "centerbody": np.column_stack((np.arange(lmax, dtype=int), np.zeros(lmax, dtype=int))),
        "wall": np.column_stack(
            (np.arange(lmax, dtype=int), (mmax - 1) * np.ones(lmax, dtype=int))
        ),
        "throat": np.column_stack((lt * np.ones(mmax, dtype=int), np.arange(mmax, dtype=int))),
    }

    if lmax > 2 and mmax > 2:
        l_int, m_int = np.meshgrid(
            np.arange(1, lmax - 1, dtype=int),
            np.arange(1, mmax - 1, dtype=int),
            indexing="ij",
        )
        boundaries["interior"] = np.column_stack((l_int.ravel(), m_int.ravel()))
    else:
        boundaries["interior"] = np.zeros((0, 2), dtype=int)

    return boundaries


def _compute_mapping_metrics(x, y):
    """Compute metric terms for the (L,M) -> (x,y) mapping."""
    lmax, mmax = x.shape
    dxi = 1.0 / (lmax - 1)
    deta = 1.0 / (mmax - 1)

    x_xi, x_eta = np.gradient(x, dxi, deta, edge_order=2)
    y_xi, y_eta = np.gradient(y, dxi, deta, edge_order=2)

    jac = x_xi * y_eta - x_eta * y_xi
    xi_x = y_eta / jac
    xi_y = -x_eta / jac
    eta_x = -y_xi / jac
    eta_y = x_xi / jac
    return jac, xi_x, xi_y, eta_x, eta_y


def _solve_mach_from_area_ratio(area_ratio, gamma, mach_guess, max_iter=30, tol=5.0e-4):
    """Newton iteration used by NAP ONEDIM for A/A* -> Mach."""
    acoef = 2.0 / (gamma + 1.0)
    bcoef = (gamma - 1.0) / (gamma + 1.0)
    ccoef = (gamma + 1.0) / (2.0 * (gamma - 1.0))

    mn = max(mach_guess, 1.0e-8)
    for _ in range(max_iter):
        abm = acoef + bcoef * mn * mn
        abmc = abm**ccoef
        fm = abmc / mn - area_ratio
        fpm = abmc * (2.0 * bcoef * ccoef / abm - 1.0 / (mn * mn))
        if abs(fpm) < 1.0e-14:
            break

        old = mn
        mn = old - fm / fpm
        if mn > 1.0 and old < 1.0:
            mn = 0.99
        if mn < 1.0 and old > 1.0:
            mn = 1.01
        if mn < 0.0:
            mn = abs(mn)
        if abs(mn - old) / max(abs(old), 1.0e-8) <= tol:
            break
    return mn


def initialize_1d_surface(
    grid,
    gamma,
    rgas,
    n1d=1,
    level=0,
    unit_system="imperial",
    pt_in_psia=True,
    tt_in_fahrenheit=True,
    gc=32.174,
    verbose=False,
):
    """
    NAP-style ONEDIM initialization for U,V,P,RO on one time level.

    n1d modes follow NAP convention:
    -2 subsonic, -1 supersonic, 1 subsonic-sonic-supersonic,
    2 subsonic-sonic-subsonic, 3 supersonic-sonic-supersonic,
    4 supersonic-sonic-subsonic.
    """
    if level not in (0, 1):
        raise ValueError("level must be 0 or 1")
    if unit_system not in ("imperial", "si"):
        raise ValueError("unit_system must be 'imperial' or 'si'")

    lmax = grid.lmax
    mmax = grid.mmax
    m1 = max(mmax - 1, 1)
    dx = (grid.xe - grid.xi) / max(lmax - 1, 1)
    gam1 = gamma / (gamma - 1.0)
    gam2 = (gamma - 1.0) / 2.0

    # Match NAP preprocessing: in English mode BC inputs are commonly psia and deg-F.
    if unit_system == "imperial":
        pt_total = grid.pt * 144.0 if pt_in_psia else grid.pt
        tt_total = grid.tt + 460.0 if tt_in_fahrenheit else grid.tt
        grgas = 1.0 / rgas
    else:
        pt_total = grid.pt
        tt_total = grid.tt
        grgas = 1.0 / rgas

    if np.any(pt_total <= 0.0) or np.any(tt_total <= 0.0):
        raise ValueError("total pressure/temperature must be positive after unit conversion")

    if verbose:
        if unit_system == "imperial":
            pt_inlet_raw = float(grid.pt[0])
            tt_inlet_raw = float(grid.tt[0])
            pt_inlet_psia = pt_inlet_raw if pt_in_psia else float(pt_total[0]) / 144.0
            tt_inlet_f = tt_inlet_raw if tt_in_fahrenheit else float(tt_total[0]) - 460.0
            print(
                "ONEDIM inputs (imperial): "
                f"PT_in={pt_inlet_psia:.4f} psia -> {float(pt_total[0]):.4f} psf, "
                f"TT_in={tt_inlet_f:.4f} F -> {float(tt_total[0]):.4f} R"
            )
        else:
            print(
                "ONEDIM inputs (si): "
                f"PT_in={float(pt_total[0]):.4f} Pa, TT_in={float(tt_total[0]):.4f} K"
            )

    has_centerbody = np.any(np.abs(grid.ycb) > 0.0)
    if n1d >= 0:
        if has_centerbody:
            rstar = grid.yw[grid.lt] - grid.ycb[grid.lt]
            rstars = grid.yw[grid.lt] ** 2 - grid.ycb[grid.lt] ** 2
        else:
            rstar = grid.rt
            rstars = grid.rt**2

    mn3 = 2.0 if (n1d == -1 or n1d > 2) else 0.01
    nxck = False
    mach_l = np.zeros(lmax)

    for l in range(lmax):
        if l == 0 and (n1d == -1 or n1d > 2):
            continue

        x = grid.xi + dx * l
        if n1d >= 0:
            if has_centerbody:
                if l < grid.lt:
                    pass
                elif l > grid.lt:
                    if not nxck:
                        mn3 = 1.1 if n1d in (1, 3) else 0.9
                        nxck = True
                else:
                    mn3 = 1.0
            else:
                if x < grid.xt:
                    pass
                elif x > grid.xt:
                    if not nxck:
                        mn3 = 1.1 if n1d in (1, 3) else 0.9
                        nxck = True
                else:
                    mn3 = 1.0

        if grid.ndim == 0:
            rad = grid.yw[l] - grid.ycb[l]
            area_ratio = rad / rstar
        else:
            rads = grid.yw[l] ** 2 - grid.ycb[l] ** 2
            area_ratio = rads / rstars

        mn3 = _solve_mach_from_area_ratio(area_ratio, gamma, mn3)
        mach_l[l] = mn3

        dem = 1.0 + gam2 * mn3 * mn3
        demp = dem**gam1
        dnxny = (grid.nxny[l] - grid.nxnycb[l]) / m1
        for m in range(mmax):
            p = pt_total[m] / demp
            temp = tt_total[m] / dem
            rho = p * grgas / temp
            q = mn3 * np.sqrt(gamma * p / rho)

            dn = grid.nxnycb[l] + dnxny * m
            dns = dn * dn
            if dns == 0.0:
                u = q
                v = 0.0
            else:
                sign = -1.0 if dn > 0.0 else 1.0
                u = q / np.sqrt(1.0 + dns)
                v = sign * q / np.sqrt(1.0 + 1.0 / dns)

            grid.p[l, m, level] = p
            grid.ro[l, m, level] = rho
            grid.u[l, m, level] = u
            grid.v[l, m, level] = v

    return mach_l


def plot_grid_with_boundaries(grid, output_path, show=True):
    """Plot structured grid with colored boundaries."""
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 4.5), dpi=150)

    for l in range(grid.lmax):
        ax.plot(grid.x[l, :], grid.y[l, :], color="0.75", lw=0.7, zorder=1)
    for m in range(grid.mmax):
        ax.plot(grid.x[:, m], grid.y[:, m], color="0.75", lw=0.7, zorder=1)

    ax.plot(grid.x[:, grid.mmax - 1], grid.y[:, grid.mmax - 1], color="#d62728", lw=2.2, label="wall", zorder=3)
    ax.plot(grid.x[:, 0], grid.y[:, 0], color="#1f77b4", lw=2.2, label="centerbody/symmetry", zorder=3)
    ax.plot(grid.x[0, :], grid.y[0, :], color="#2ca02c", lw=2.2, label="inlet", zorder=3)
    ax.plot(grid.x[grid.lmax - 1, :], grid.y[grid.lmax - 1, :], color="#ff7f0e", lw=2.2, label="exit", zorder=3)
    ax.plot(grid.x[grid.lt, :], grid.y[grid.lt, :], color="#9467bd", lw=1.8, ls="--", label="throat index", zorder=3)

    ax.set_title("NAP-style Structured Grid and Boundaries")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.15)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    if show and "agg" not in plt.get_backend().lower():
        plt.show()
    plt.close(fig)


def plot_boundary_indices(grid, output_path, show=True):
    """Plot boundary points with (L,M) labels for BC/debug work."""
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 5), dpi=150)
    ax.plot(grid.x, grid.y, color="0.85", lw=0.5, zorder=1)
    ax.plot(grid.x.T, grid.y.T, color="0.85", lw=0.5, zorder=1)

    boundary_specs = [
        ("wall", "#d62728"),
        ("centerbody", "#1f77b4"),
        ("inlet", "#2ca02c"),
        ("exit", "#ff7f0e"),
    ]

    for name, color in boundary_specs:
        lm = grid.boundary_lm[name]
        l_idx = lm[:, 0]
        m_idx = lm[:, 1]
        xpts = grid.x[l_idx, m_idx]
        ypts = grid.y[l_idx, m_idx]
        ax.plot(xpts, ypts, color=color, lw=2.2, label=name, zorder=3)
        ax.scatter(xpts, ypts, color=color, s=14, zorder=4)

        for l, m, xp, yp in zip(l_idx, m_idx, xpts, ypts):
            ax.text(
                xp,
                yp,
                f"({l},{m})",
                fontsize=6,
                color=color,
                ha="left",
                va="bottom",
                zorder=5,
            )

    l_throat = grid.lt
    x_th = grid.x[l_throat, :]
    y_th = grid.y[l_throat, :]
    ax.plot(x_th, y_th, color="#9467bd", ls="--", lw=1.8, label="throat index", zorder=3)

    ax.set_title("Boundary Index Map (0-based L,M)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.15)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    if show and "agg" not in plt.get_backend().lower():
        plt.show()
    plt.close(fig)


def plot_1d_solution_surface(
    grid,
    gamma,
    level,
    output_path,
    show=True,
    unit_system="imperial",
):
    """Plot ONEDIM diagnostics without interpolation-heavy contours."""
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    u = grid.u[:, :, level]
    v = grid.v[:, :, level]
    p = grid.p[:, :, level]
    ro = grid.ro[:, :, level]
    a = np.sqrt(gamma * p / ro)
    mach = np.sqrt(u * u + v * v) / a

    x = grid.x[:, 0]
    j_c = 0
    j_w = grid.mmax - 1
    i_t = grid.lt
    i_e = grid.lmax - 1

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=150)
    ax0, ax1 = axes

    if unit_system == "imperial":
        p_plot = p / 144.0
        p_label = "Pressure (psia)"
    else:
        p_plot = p
        p_label = "Pressure (Pa)"

    # Panel 1: axial traces (centerline and wall)
    ax0.plot(x, mach[:, j_c], label="Mach centerline", color="#1f77b4", lw=2.0)
    ax0.plot(x, mach[:, j_w], label="Mach wall", color="#17becf", lw=2.0)
    ax0.axvline(x[i_t], color="0.35", ls="--", lw=1.2, label="throat x")
    ax0.set_title("Mach vs Axial Position")
    ax0.set_xlabel("x")
    ax0.set_ylabel("Mach")
    ax0.grid(alpha=0.25)
    ax0.legend(loc="best")

    # Panel 2: radial profiles at inlet/throat/exit
    y_in = grid.y[0, :]
    y_th = grid.y[i_t, :]
    y_ex = grid.y[i_e, :]
    ax1.plot(y_in, p_plot[0, :], label="P inlet", color="#d62728", lw=2.0)
    ax1.plot(y_th, p_plot[i_t, :], label="P throat", color="#ff7f0e", lw=2.0)
    ax1.plot(y_ex, p_plot[i_e, :], label="P exit", color="#9467bd", lw=2.0)
    ax1.set_title("Pressure Radial Profiles")
    ax1.set_xlabel("y")
    ax1.set_ylabel(p_label)
    ax1.grid(alpha=0.25)
    ax1.legend(loc="best")

    # Print quick numeric ranges to stdout for sanity checks.
    print(
        "1D plot stats: "
        f"Mach[min,max]=({float(np.min(mach)):.4f},{float(np.max(mach)):.4f}), "
        f"{p_label}[min,max]=({float(np.min(p_plot)):.4f},{float(np.max(p_plot)):.4f})"
    )

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    if show and "agg" not in plt.get_backend().lower():
        plt.show()
    plt.close(fig)


def plot_residual_history(history, output_path, show=True):
    """Plot max relative residual vs iteration."""
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 4.0), dpi=150)
    hist = np.asarray(history, dtype=float)
    hist = np.nan_to_num(hist, nan=np.nan, posinf=np.nan, neginf=np.nan)
    finite = np.isfinite(hist) & (hist > 0.0)
    if np.any(finite):
        it = np.arange(1, len(hist) + 1)[finite]
        ax.semilogy(it, hist[finite], color="#1f77b4", lw=1.8)
    else:
        ax.plot([], [])
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Max relative residual")
    ax.set_title("Transient March to Steady-State")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    if show and "agg" not in plt.get_backend().lower():
        plt.show()
    plt.close(fig)


def _dd_comp(arr, axis, scheme, d):
    """First derivative in computational coordinates with one-sided/centered fallback."""
    out = np.zeros_like(arr)
    if axis == 0:
        if scheme == "forward":
            out[:-1, :] = (arr[1:, :] - arr[:-1, :]) / d
            out[-1, :] = (arr[-1, :] - arr[-2, :]) / d
        else:
            out[1:, :] = (arr[1:, :] - arr[:-1, :]) / d
            out[0, :] = (arr[1, :] - arr[0, :]) / d
    else:
        if scheme == "forward":
            out[:, :-1] = (arr[:, 1:] - arr[:, :-1]) / d
            out[:, -1] = (arr[:, -1] - arr[:, -2]) / d
        else:
            out[:, 1:] = (arr[:, 1:] - arr[:, :-1]) / d
            out[:, 0] = (arr[:, 1] - arr[:, 0]) / d
    return out


def _laplace_comp(arr, dxi, deta):
    """Simple Laplacian in computational space for numerical smoothing."""
    lap = np.zeros_like(arr)
    lap[1:-1, 1:-1] = (
        (arr[2:, 1:-1] - 2.0 * arr[1:-1, 1:-1] + arr[:-2, 1:-1]) / (dxi * dxi)
        + (arr[1:-1, 2:] - 2.0 * arr[1:-1, 1:-1] + arr[1:-1, :-2]) / (deta * deta)
    )
    return lap


def _rhs_nonconservative(grid, rho, u, v, p, gamma, scheme, cav=0.0):
    """
    NAP-like non-conservative RHS in physical space using mapping metrics.
    This is a pragmatic baseline, not a full characteristic-boundary replica.
    """
    lmax, mmax = rho.shape
    dxi = 1.0 / max(lmax - 1, 1)
    deta = 1.0 / max(mmax - 1, 1)

    rho_xi = _dd_comp(rho, axis=0, scheme=scheme, d=dxi)
    rho_eta = _dd_comp(rho, axis=1, scheme=scheme, d=deta)
    u_xi = _dd_comp(u, axis=0, scheme=scheme, d=dxi)
    u_eta = _dd_comp(u, axis=1, scheme=scheme, d=deta)
    v_xi = _dd_comp(v, axis=0, scheme=scheme, d=dxi)
    v_eta = _dd_comp(v, axis=1, scheme=scheme, d=deta)
    p_xi = _dd_comp(p, axis=0, scheme=scheme, d=dxi)
    p_eta = _dd_comp(p, axis=1, scheme=scheme, d=deta)

    rho_x = grid.xi_x * rho_xi + grid.eta_x * rho_eta
    rho_y = grid.xi_y * rho_xi + grid.eta_y * rho_eta
    u_x = grid.xi_x * u_xi + grid.eta_x * u_eta
    u_y = grid.xi_y * u_xi + grid.eta_y * u_eta
    v_x = grid.xi_x * v_xi + grid.eta_x * v_eta
    v_y = grid.xi_y * v_xi + grid.eta_y * v_eta
    p_x = grid.xi_x * p_xi + grid.eta_x * p_eta
    p_y = grid.xi_y * p_xi + grid.eta_y * p_eta

    eps = 1.0 if grid.ndim == 1 else 0.0
    y_safe = np.where(np.abs(grid.y) > 1.0e-10, grid.y, 1.0)
    div = u_x + v_y + eps * v / y_safe

    rho_rhs = -(u * rho_x + v * rho_y + rho * div)
    u_rhs = -(u * u_x + v * u_y + p_x / rho)
    v_rhs = -(u * v_x + v * v_y + p_y / rho)
    p_rhs = -(u * p_x + v * p_y + gamma * p * div)

    if cav > 0.0:
        rho_rhs += cav * _laplace_comp(rho, dxi, deta)
        u_rhs += cav * _laplace_comp(u, dxi, deta)
        v_rhs += cav * _laplace_comp(v, dxi, deta)
        p_rhs += cav * _laplace_comp(p, dxi, deta)

    return rho_rhs, u_rhs, v_rhs, p_rhs


def _solve_inlet_mach_from_characteristic(tt, pt, theta, j_minus, gamma, r_eff):
    """Solve subsonic inlet Mach from total conditions and outgoing characteristic."""
    gam2 = (gamma - 1.0) / 2.0
    gamp = gamma / (gamma - 1.0)

    def f(m):
        dem = 1.0 + gam2 * m * m
        t = tt / dem
        a = np.sqrt(gamma * r_eff * t)
        q = m * a
        return q * np.cos(theta) - 2.0 * a / (gamma - 1.0) - j_minus

    m = 0.15
    for _ in range(30):
        fm = f(m)
        dm = 1.0e-5
        fp = f(m + dm)
        df = (fp - fm) / dm
        if abs(df) < 1.0e-12:
            break
        m_new = m - fm / df
        m = float(np.clip(m_new, 1.0e-4, 0.99))
        if abs(fm) < 1.0e-7:
            break

    dem = 1.0 + gam2 * m * m
    t = tt / dem
    a = np.sqrt(gamma * r_eff * t)
    q = m * a
    p = pt / (dem**gamp)
    rho = p / (r_eff * t)
    u = q * np.cos(theta)
    v = -q * np.sin(theta)
    return rho, u, v, p


def _inlet_update_nap(grid, rho, u, v, p, inlet_state, gamma, dt, p_floor, rho_floor):
    """Approximate line-by-line port of NAP INLET subsonic characteristic update."""
    lmax, mmax = rho.shape
    if lmax < 2:
        return

    x3 = grid.xi
    dx = (grid.xe - grid.xi) / max(lmax - 1, 1)
    dxr = 1.0 / max(dx, 1.0e-12)
    dy = 1.0 / max(mmax - 1, 1)
    dyr = 1.0 / max(dy, 1.0e-12)
    p_cap = 1.0e7
    rho_cap = 1.0e3
    u_cap = 1.0e4
    span = max(grid.xe - grid.xi, 1.0)

    for ichar in (0, 1):
        for m in range(mmax):
            al, be, _ = _map_terms_nap(grid, 0, m, dy)
            al1, be1, _ = _map_terms_nap(grid, 1, m, dy)

            u2 = u[1, m]
            v2 = v[1, m]
            p2 = max(p[1, m], p_floor)
            ro2 = max(rho[1, m], rho_floor)
            a2 = np.sqrt(max(gamma * p2 / ro2, 1.0e-12))

            if ichar == 0:
                u[0, m] = u2
                v[0, m] = v2
                a3 = a2
            else:
                a3 = np.sqrt(max(gamma * p[0, m] / max(rho[0, m], rho_floor), 1.0e-12))

            bu = (u[1, m] - u[0, m]) * dxr
            bv = (v[1, m] - v[0, m]) * dxr
            bp = (p[1, m] - p[0, m]) * dxr
            bro = (rho[1, m] - rho[0, m]) * dxr
            bal = (al1 - al) * dxr
            bbe = (be1 - be) * dxr
            cu = u[0, m] - bu * x3
            cv = v[0, m] - bv * x3
            cp = p[0, m] - bp * x3
            cro = rho[0, m] - bro * x3
            cal = al - bal * x3
            cbe = be - bbe * x3

            if m > 0:
                du = (u[1, m] - u[1, m - 1]) * dyr
                dv = (v[1, m] - v[1, m - 1]) * dyr
                dp = (p[1, m] - p[1, m - 1]) * dyr
                dro = (rho[1, m] - rho[1, m - 1]) * dyr
                du1 = (u[0, m] - u[0, m - 1]) * dyr
                dv1 = (v[0, m] - v[0, m - 1]) * dyr
                dp1 = (p[0, m] - p[0, m - 1]) * dyr
                dro1 = (rho[0, m] - rho[0, m - 1]) * dyr
            else:
                du = 0.0
                dv = v[1, min(1, mmax - 1)] * dyr if mmax > 1 else 0.0
                dp = 0.0
                dro = 0.0
                du1 = 0.0
                dv1 = v[0, min(1, mmax - 1)] * dyr if mmax > 1 else 0.0
                dp1 = 0.0
                dro1 = 0.0

            bdu = (du - du1) * dxr
            bdv = (dv - dv1) * dxr
            bdp = (dp - dp1) * dxr
            bdro = (dro - dro1) * dxr
            cdu = du1 - bdu * x3
            cdv = dv1 - bdv * x3
            cdp = dp1 - bdp * x3
            cdro = dro1 - bdro * x3

            # Characteristic foot interpolation for x2.
            for _ in range(2):
                x2 = x3 - (u[0, m] - a3 + u2 - a2) * 0.5 * dt
                x2 = float(np.clip(x2, grid.xi - 2.0 * span, grid.xe + 2.0 * span))
                u2 = bu * x2 + cu
                p2 = float(np.clip(bp * x2 + cp, p_floor, p_cap))
                ro2 = float(np.clip(bro * x2 + cro, rho_floor, rho_cap))
                a2 = np.sqrt(max(gamma * p2 / ro2, 1.0e-12))

            v2 = bv * x2 + cv
            al2 = bal * x2 + cal
            be2 = bbe * x2 + cbe
            uv2 = u2 * al2 + v2 * be2
            du2 = bdu * x2 + cdu
            dv2 = bdv * x2 + cdv
            dp2 = bdp * x2 + cdp
            dro2 = bdro * x2 + cdro

            aterm2 = 0.0
            if grid.ndim == 1:
                ycb2 = grid.ycb[0]
                den2 = max((dy * m) / max(be2, 1.0e-12) + ycb2, 1.0e-12)
                aterm2 = ro2 * v2 / den2

            psi12 = -uv2 * dro2 - ro2 * al2 * du2 - ro2 * be2 * dv2 - aterm2
            psi22 = -uv2 * du2 - al2 * dp2 / ro2
            psi42 = -uv2 * dp2 + a2 * a2 * uv2 * dro2

            if ichar == 1:
                if m == 0:
                    du3 = 0.0
                    dv3 = v[0, min(1, mmax - 1)] * dyr if mmax > 1 else 0.0
                    dp3 = 0.0
                    dro3 = 0.0
                elif m == mmax - 1:
                    du3 = (u[0, m] - u[0, m - 1]) * dyr
                    dv3 = (v[0, m] - v[0, m - 1]) * dyr
                    dp3 = (p[0, m] - p[0, m - 1]) * dyr
                    dro3 = (rho[0, m] - rho[0, m - 1]) * dyr
                else:
                    du3 = (u[0, m + 1] - u[0, m]) * dyr
                    dv3 = (v[0, m + 1] - v[0, m]) * dyr
                    dp3 = (p[0, m + 1] - p[0, m]) * dyr
                    dro3 = (rho[0, m + 1] - rho[0, m]) * dyr

                aterm3 = 0.0
                if grid.ndim == 1:
                    den3 = max((dy * m) / max(be, 1.0e-12) + grid.ycb[0], 1.0e-12)
                    aterm3 = rho[0, m] * v[0, m] / den3
                uv3 = u[0, m] * al + v[0, m] * be
                psi13 = -uv3 * dro3 - rho[0, m] * al * du3 - rho[0, m] * be * dv3 - aterm3
                psi23 = -uv3 * du3 - al * dp3 / max(rho[0, m], rho_floor)
                psi43 = -uv3 * dp3 + a3 * a3 * uv3 * dro3
            else:
                psi13, psi23, psi43 = psi12, psi22, psi42

            psi1b = 0.5 * (psi12 + psi13)
            psi2b = 0.5 * (psi22 + psi23)
            psi4b = 0.5 * (psi42 + psi43)
            gpsi1b = gamma * psi1b
            theta = inlet_state["theta"][m]
            ttheta = np.tan(theta)
            ucorr = 0.5 + 0.5 / np.sqrt(1.0 + ttheta * ttheta)

            mn3 = np.sqrt(max((u[0, m] * u[0, m] + v[0, m] * v[0, m]), 0.0)) / max(a3, 1.0e-8)
            t2 = p2 / (ro2 * inlet_state["r_eff"])
            for _ in range(20):
                dem = 1.0 + 0.5 * (gamma - 1.0) * mn3 * mn3
                p_new = inlet_state["pt"][m] / (dem ** (gamma / (gamma - 1.0)))
                t3 = inlet_state["tt"][m] / dem
                pb = 0.5 * (p2 + p_new)
                rtb = inlet_state["r_eff"] * 0.5 * (t2 + t3)
                u_new = u2 + dt * psi2b + (p_new - p2 - (psi4b + rtb * gpsi1b) * dt) * np.sqrt(
                    max(rtb / gamma, 1.0e-12)
                ) / max(pb, p_floor)
                u_new = u_new * ucorr
                v_new = -u_new * ttheta
                u_new = float(np.clip(u_new, -u_cap, u_cap))
                v_new = float(np.clip(v_new, -u_cap, u_cap))
                mn_old = mn3
                mn3 = np.sqrt(max(u_new * u_new + v_new * v_new, 0.0) / max(t3 * gamma * inlet_state["r_eff"], 1.0e-12))
                if mn_old == 0.0:
                    if abs(mn3 - mn_old) <= 1.0e-4:
                        break
                elif abs((mn3 - mn_old) / mn_old) <= 1.0e-3:
                    break

            p[0, m] = float(np.clip(p_new, p_floor, p_cap))
            rho[0, m] = float(np.clip(p[0, m] / max(inlet_state["r_eff"] * t3, 1.0e-12), rho_floor, rho_cap))
            u[0, m] = u_new
            v[0, m] = v_new


def _exit_update_nap(grid, rho, u, v, p, p_exit, gamma, dt, p_floor, rho_floor):
    """Approximate line-by-line port of NAP EXITT subsonic characteristic update."""
    lmax, mmax = rho.shape
    if lmax < 3:
        return

    x3 = grid.xe
    dx = (grid.xe - grid.xi) / max(lmax - 1, 1)
    dxr = 1.0 / max(dx, 1.0e-12)
    dy = 1.0 / max(mmax - 1, 1)
    dyr = 1.0 / max(dy, 1.0e-12)
    p_cap = 1.0e7
    rho_cap = 1.0e3
    u_cap = 1.0e4
    span = max(grid.xe - grid.xi, 1.0)
    l = lmax - 1
    l1 = lmax - 2

    for ichar in (0, 1):
        for m in range(mmax):
            al, be, de = _map_terms_nap(grid, l, m, dy)
            al1, be1, de1 = _map_terms_nap(grid, l1, m, dy)

            u1 = u[l, m]
            u2 = u1
            p1 = max(p[l, m], p_floor)
            ro1 = max(rho[l, m], rho_floor)
            a1 = np.sqrt(max(gamma * p1 / ro1, 1.0e-12))
            a2 = a1

            if ichar == 0:
                u[l, m] = u1
                a3 = a1
            else:
                a3 = np.sqrt(max(gamma * max(p[l, m], p_floor) / max(rho[l, m], rho_floor), 1.0e-12))

            # Property interpolation coefficients along x.
            bu = (u[l, m] - u[l1, m]) * dxr
            bv = (v[l, m] - v[l1, m]) * dxr
            bp = (p[l, m] - p[l1, m]) * dxr
            bro = (rho[l, m] - rho[l1, m]) * dxr
            bal = (al - al1) * dxr
            bbe = (be - be1) * dxr
            bde = (de - de1) * dxr
            cu = u[l, m] - bu * x3
            cv = v[l, m] - bv * x3
            cp = p[l, m] - bp * x3
            cro = rho[l, m] - bro * x3
            cal = al - bal * x3
            cbe = be - bbe * x3
            cde = de - bde * x3

            # Cross-derivative interpolation coefficients.
            if m > 0:
                du = (u[l, m] - u[l, m - 1]) * dyr
                dv = (v[l, m] - v[l, m - 1]) * dyr
                dp = (p[l, m] - p[l, m - 1]) * dyr
                dro = (rho[l, m] - rho[l, m - 1]) * dyr
                du1 = (u[l1, m] - u[l1, m - 1]) * dyr
                dv1 = (v[l1, m] - v[l1, m - 1]) * dyr
                dp1 = (p[l1, m] - p[l1, m - 1]) * dyr
                dro1 = (rho[l1, m] - rho[l1, m - 1]) * dyr
            else:
                du = 0.0
                dv = v[l, min(1, mmax - 1)] * dyr if mmax > 1 else 0.0
                dp = 0.0
                dro = 0.0
                du1 = 0.0
                dv1 = v[l1, min(1, mmax - 1)] * dyr if mmax > 1 else 0.0
                dp1 = 0.0
                dro1 = 0.0

            bdu = (du - du1) * dxr
            bdv = (dv - dv1) * dxr
            bdp = (dp - dp1) * dxr
            bdro = (dro - dro1) * dxr
            cdu = du - bdu * x3
            cdv = dv - bdv * x3
            cdp = dp - bdp * x3
            cdro = dro - bdro * x3

            # Compute x1/x2 characteristic feet.
            for _ in range(2):
                x1 = x3 - (u[l, m] - u1) * 0.5 * dt
                x2 = x3 - (u[l, m] - a3 + u2 - a2) * 0.5 * dt
                x1 = float(np.clip(x1, grid.xi - 2.0 * span, grid.xe + 2.0 * span))
                x2 = float(np.clip(x2, grid.xi - 2.0 * span, grid.xe + 2.0 * span))
                u1 = bu * x1 + cu
                u2 = bu * x2 + cu
                p2 = float(np.clip(bp * x2 + cp, p_floor, p_cap))
                ro2 = float(np.clip(bro * x2 + cro, rho_floor, rho_cap))
                a2 = np.sqrt(max(gamma * p2 / ro2, 1.0e-12))

            v1 = bv * x1 + cv
            p1 = float(np.clip(bp * x1 + cp, p_floor, p_cap))
            ro1 = float(np.clip(bro * x1 + cro, rho_floor, rho_cap))
            ycb1 = grid.ycb[l]
            al_1 = bal * x1 + cal
            be_1 = bbe * x1 + cbe
            de_1 = bde * x1 + cde
            uv1 = u1 * al_1 + v1 * be_1 + de_1
            a1 = np.sqrt(max(gamma * p1 / ro1, 1.0e-12))

            v2 = bv * x2 + cv
            ycb2 = grid.ycb[l]
            al_2 = bal * x2 + cal
            be_2 = bbe * x2 + cbe
            de_2 = bde * x2 + cde
            uv2 = u2 * al_2 + v2 * be_2 + de_2

            du1 = bdu * x1 + cdu
            dv1 = bdv * x1 + cdv
            dp1 = bdp * x1 + cdp
            dro1 = bdro * x1 + cdro
            du2 = bdu * x2 + cdu
            dv2 = bdv * x2 + cdv
            dp2 = bdp * x2 + cdp
            dro2 = bdro * x2 + cdro

            # PSI terms at interpolation points.
            if grid.ndim == 1:
                den1 = max((dy * m) / max(be_1, 1.0e-12) + ycb1, 1.0e-12)
                den2 = max((dy * m) / max(be_2, 1.0e-12) + ycb2, 1.0e-12)
                aterm1 = ro1 * v1 / den1
                aterm2 = ro2 * v2 / den2
            else:
                aterm1 = 0.0
                aterm2 = 0.0

            psi31 = -uv1 * dv1 - be_1 * dp1 / ro1
            psi41 = -uv1 * dp1 + a1 * a1 * uv1 * dro1
            psi12 = -uv2 * dro2 - ro2 * al_2 * du2 - ro2 * be_2 * dv2 - aterm2
            psi22 = -uv2 * du2 - al_2 * dp2 / ro2
            psi42 = -uv2 * dp2 + a2 * a2 * uv2 * dro2

            if ichar == 1:
                if m == 0:
                    du3 = 0.0
                    dv3 = v[l, min(1, mmax - 1)] * dyr if mmax > 1 else 0.0
                    dp3 = 0.0
                    dro3 = 0.0
                elif m == mmax - 1:
                    du3 = (u[l, m] - u[l, m - 1]) * dyr
                    dv3 = (v[l, m] - v[l, m - 1]) * dyr
                    dp3 = (p[l, m] - p[l, m - 1]) * dyr
                    dro3 = (rho[l, m] - rho[l, m - 1]) * dyr
                else:
                    du3 = (u[l, m + 1] - u[l, m]) * dyr
                    dv3 = (v[l, m + 1] - v[l, m]) * dyr
                    dp3 = (p[l, m + 1] - p[l, m]) * dyr
                    dro3 = (rho[l, m + 1] - rho[l, m]) * dyr

                if grid.ndim == 1:
                    den3 = max((dy * m) / max(be, 1.0e-12) + grid.ycb[l], 1.0e-12)
                    aterm3 = rho[l, m] * v[l, m] / den3
                else:
                    aterm3 = 0.0
                uv3 = u[l, m] * al + v[l, m] * be + de
                psi13 = -uv3 * dro3 - rho[l, m] * (al * du3 + be * dv3) - aterm3
                psi23 = -uv3 * du3 - al * dp3 / max(rho[l, m], rho_floor)
                psi33 = -uv3 * dv3 - be * dp3 / max(rho[l, m], rho_floor)
                psi43 = -uv3 * dp3 + a3 * a3 * uv3 * dro3
            else:
                psi13, psi23, psi33, psi43 = psi12, psi22, psi31, psi42

            if ichar == 1:
                psi31b = 0.5 * (psi31 + psi33)
                psi41b = 0.5 * (psi41 + psi43)
                psi12b = 0.5 * (psi12 + psi13)
                psi22b = 0.5 * (psi22 + psi23)
                psi42b = 0.5 * (psi42 + psi43)
            else:
                psi31b, psi41b = psi31, psi41
                psi12b, psi22b, psi42b = psi12, psi22, psi42

            # Solve compatibility equations with prescribed exit pressure.
            p[l, m] = float(np.clip(max(p_exit, p_floor), p_floor, p_cap))
            rho[l, m] = float(
                np.clip(ro1 + (p[l, m] - p1 - dt * psi41b) / max(a3 * a3 + a1 * a1, 1.0e-12), rho_floor, rho_cap)
            )
            num = (
                psi42b
                + (ro2 - rho[l, m]) * (a2 + a3) * psi22b / 4.0
                - (a2 * a2 + a3 * a3) * psi12b / 2.0
            ) * dt - (p[l, m] - p2)
            den = max((ro2 + rho[l, m]) * (a2 + a3), 1.0e-12)
            u[l, m] = float(np.clip(u2 + 4.0 * num / den, -u_cap, u_cap))
            v[l, m] = float(np.clip(v1 + dt * psi31b, -u_cap, u_cap))


def _wall_update_nap(grid, rho, u, v, p, gamma, dt, p_floor, rho_floor):
    """
    NAP WALL-style compatibility update for the nozzle wall (upper boundary only).
    This excludes exhaust-jet moving-boundary logic for now.
    """
    lmax, mmax = rho.shape
    if lmax < 3 or mmax < 2:
        return

    dx = (grid.xe - grid.xi) / max(lmax - 1, 1)
    dxr = 1.0 / max(dx, 1.0e-12)
    dy = 1.0 / max(mmax - 1, 1)
    dyr = 1.0 / max(dy, 1.0e-12)
    m_w = mmax - 1
    m_i = mmax - 2
    sign = -1.0
    xwid = 0.0
    p_cap = 1.0e7
    rho_cap = 1.0e3
    u_cap = 1.0e4

    for l in range(1, lmax - 1):
        al, be, de = _map_terms_nap(grid, l, m_w, dy)
        als = np.sqrt(max(al * al + be * be, 1.0e-12))
        abr = grid.nxny[l]

        u1 = u[l, m_w]
        v1 = v[l, m_w]
        p1 = max(p[l, m_w], p_floor)
        ro1 = max(rho[l, m_w], rho_floor)
        a1 = np.sqrt(max(gamma * p1 / ro1, 1.0e-12))

        u2 = u[l, m_i]
        v2 = v[l, m_i]
        p2 = max(p[l, m_i], p_floor)
        ro2 = max(rho[l, m_i], rho_floor)
        a2 = np.sqrt(max(gamma * p2 / ro2, 1.0e-12))

        u3 = u1
        v3 = v1
        p3 = p1
        ro3 = ro1
        a3 = a1

        # X-derivatives at wall and adjacent interior row.
        du1 = (u1 - u[l - 1, m_w]) * dxr
        dv1 = (v1 - v[l - 1, m_w]) * dxr
        dp1 = (p1 - p[l - 1, m_w]) * dxr
        dro1 = (ro1 - rho[l - 1, m_w]) * dxr
        du2 = (u2 - u[l - 1, m_i]) * dxr
        dv2 = (v2 - v[l - 1, m_i]) * dxr
        dp2 = (p2 - p[l - 1, m_i]) * dxr
        dro2 = (ro2 - rho[l - 1, m_i]) * dxr

        aterm2 = 0.0
        if grid.ndim == 1:
            den2 = max(grid.ycb[l] + 1.0 / max(be, 1.0e-12), 1.0e-12)
            aterm2 = ro2 * v2 / den2

        psi21 = -u1 * du1 - dp1 / ro1
        psi31 = -u1 * dv1
        psi41 = -u1 * dp1 + a1 * a1 * u1 * dro1

        psi12 = -u2 * dro2 - ro2 * du2 - aterm2
        psi22 = -u2 * du2 - dp2 / ro2
        psi32 = -u2 * dv2
        psi42 = -u2 * dp2 + a2 * a2 * u2 * dro2

        # Approximate solution-point x-derivatives for second pass.
        du3 = (u[l + 1, m_w] - u3) * dxr
        dv3 = (v[l + 1, m_w] - v3) * dxr
        dp3 = (p[l + 1, m_w] - p3) * dxr
        dro3 = (rho[l + 1, m_w] - ro3) * dxr

        aterm3 = 0.0
        if grid.ndim == 1:
            den3 = max(grid.ycb[l] + 1.0 / max(be, 1.0e-12), 1.0e-12)
            aterm3 = ro3 * v3 / den3
        psi13 = -u3 * dro3 - ro3 * du3 - aterm3
        psi23 = -u3 * du3 - dp3 / ro3
        psi33 = -u3 * dv3
        psi43 = -u3 * dp3 + a3 * a3 * u3 * dro3

        psi21b = 0.5 * (psi21 + psi23)
        psi31b = 0.5 * (psi31 + psi33)
        psi41b = 0.5 * (psi41 + psi43)
        psi12b = 0.5 * (psi12 + psi13)
        psi22b = 0.5 * (psi22 + psi23)
        psi32b = 0.5 * (psi32 + psi33)
        psi42b = 0.5 * (psi42 + psi43)

        alb = al / als
        beb = be / als
        a1b = 0.5 * (a1 + a3)
        a2b = 0.5 * (a2 + a3)
        ro2b = 0.5 * (ro2 + ro3)

        u_new = (u1 - abr * (v1 - xwid) + (psi21b - abr * psi31b) * dt) / (1.0 + abr * abr)
        v_new = -u_new * abr + xwid
        u_new = float(np.clip(u_new, -u_cap, u_cap))
        v_new = float(np.clip(v_new, -u_cap, u_cap))

        p_new = p2 - sign * ro2b * a2b * (alb * (u_new - u2) + beb * (v_new - v2))
        p_new += (psi42b + a2b * a2b * psi12b + sign * ro2b * a2b * (alb * psi22b + beb * psi32b)) * dt
        p_new = float(np.clip(p_new, p_floor, p_cap))

        ro_new = ro1 + (p_new - p1 - psi41b * dt) / max(a1b * a1b, 1.0e-12)
        ro_new = float(np.clip(ro_new, rho_floor, rho_cap))

        u[l, m_w] = u_new
        v[l, m_w] = v_new
        p[l, m_w] = p_new
        rho[l, m_w] = ro_new

    # End points: tangency/extrapolation fallback.
    u[0, m_w] = u[0, m_i]
    v[0, m_w] = -u[0, m_w] * grid.nxny[0]
    p[0, m_w] = p[0, m_i]
    rho[0, m_w] = rho[0, m_i]
    u[-1, m_w] = u[-1, m_i]
    v[-1, m_w] = -u[-1, m_w] * grid.nxny[-1]
    p[-1, m_w] = p[-1, m_i]
    rho[-1, m_w] = rho[-1, m_i]


def _apply_boundaries(
    grid,
    rho,
    u,
    v,
    p,
    inlet_state,
    p_exit,
    p_floor,
    rho_floor,
    gamma,
    use_characteristics=False,
    dt=None,
):
    """Apply practical NAP-like boundaries on primitive fields."""
    lmax, mmax = rho.shape
    rho[:, :] = np.nan_to_num(rho, nan=rho_floor, posinf=rho_floor, neginf=rho_floor)
    p[:, :] = np.nan_to_num(p, nan=p_floor, posinf=p_floor, neginf=p_floor)
    u[:, :] = np.nan_to_num(u, nan=0.0, posinf=0.0, neginf=0.0)
    v[:, :] = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
    rho[:, :] = np.maximum(rho, rho_floor)
    p[:, :] = np.maximum(p, p_floor)

    if use_characteristics:
        if dt is not None and dt > 0.0:
            _inlet_update_nap(grid, rho, u, v, p, inlet_state, gamma, dt, p_floor, rho_floor)
        else:
            # Fallback inlet characteristic when no timestep is available.
            for m in range(mmax):
                j_minus = u[1, m] - 2.0 * np.sqrt(max(gamma * p[1, m] / rho[1, m], 1.0e-12)) / (gamma - 1.0)
                rho_i, u_i, v_i, p_i = _solve_inlet_mach_from_characteristic(
                    inlet_state["tt"][m],
                    inlet_state["pt"][m],
                    inlet_state["theta"][m],
                    j_minus,
                    gamma,
                    inlet_state["r_eff"],
                )
                rho[0, m] = rho_i
                u[0, m] = u_i
                v[0, m] = v_i
                p[0, m] = p_i

        if dt is not None and dt > 0.0:
            _exit_update_nap(grid, rho, u, v, p, p_exit, gamma, dt, p_floor, rho_floor)
        else:
            # Exit fallback if dt isn't available.
            p[-1, :] = p_exit
            for m in range(mmax):
                j_plus = u[-2, m] + 2.0 * np.sqrt(max(gamma * p[-2, m] / rho[-2, m], 1.0e-12)) / (gamma - 1.0)
                rho_m = max(rho[-2, m], rho_floor)
                p_m = max(p[-2, m], p_floor)
                s = p_m / max(rho_m**gamma, 1.0e-12)
                rho_e = max((p_exit / max(s, 1.0e-12)) ** (1.0 / gamma), rho_floor)
                a_e = np.sqrt(max(gamma * p_exit / rho_e, 1.0e-12))
                u_e = j_plus - 2.0 * a_e / (gamma - 1.0)
                rho[-1, m] = rho_e
                u[-1, m] = u_e
                v[-1, m] = v[-2, m]
    else:
        # Inlet: hold to initial surface profile.
        rho[0, :] = inlet_state["rho"]
        u[0, :] = inlet_state["u"]
        v[0, :] = inlet_state["v"]
        p[0, :] = inlet_state["p"]

        # Exit: fixed static pressure, extrapolate others.
        p[-1, :] = p_exit
        rho[-1, :] = rho[-2, :]
        u[-1, :] = u[-2, :]
        v[-1, :] = v[-2, :]

    # Centerline/symmetry.
    v[:, 0] = 0.0
    rho[:, 0] = rho[:, 1]
    u[:, 0] = u[:, 1]
    p[:, 0] = p[:, 1]

    # Wall update.
    if use_characteristics and dt is not None and dt > 0.0:
        _wall_update_nap(grid, rho, u, v, p, gamma, dt, p_floor, rho_floor)
    else:
        wall_slope = grid.nxny
        u[:, -1] = u[:, -2]
        v[:, -1] = -u[:, -1] * wall_slope
        rho[:, -1] = rho[:, -2]
        p[:, -1] = p[:, -2]

    # Positivity floors.
    rho[:, :] = np.maximum(rho, rho_floor)
    p[:, :] = np.maximum(p, p_floor)

    # Guard corners after all BCs.
    rho[0, 0] = max(rho[0, 0], rho_floor)
    p[0, 0] = max(p[0, 0], p_floor)
    rho[lmax - 1, mmax - 1] = max(rho[lmax - 1, mmax - 1], rho_floor)
    p[lmax - 1, mmax - 1] = max(p[lmax - 1, mmax - 1], p_floor)


def _compute_cfl_timestep(grid, rho, u, v, p, gamma, fdt):
    """Compute explicit CFL time step from local wave speeds."""
    a = np.sqrt(np.maximum(gamma * p / rho, 1.0e-12))
    lmax, mmax = rho.shape

    dx = np.zeros_like(rho)
    dy = np.zeros_like(rho)
    dx[:-1, :] = np.sqrt((grid.x[1:, :] - grid.x[:-1, :]) ** 2 + (grid.y[1:, :] - grid.y[:-1, :]) ** 2)
    dx[-1, :] = dx[-2, :]
    dy[:, :-1] = np.sqrt((grid.x[:, 1:] - grid.x[:, :-1]) ** 2 + (grid.y[:, 1:] - grid.y[:, :-1]) ** 2)
    dy[:, -1] = dy[:, -2]

    dt_x = dx / (np.abs(u) + a + 1.0e-12)
    dt_y = dy / (np.abs(v) + a + 1.0e-12)
    dt_raw = np.minimum(dt_x, dt_y)
    dt_raw = np.nan_to_num(dt_raw, nan=1.0e-10, posinf=1.0e-10, neginf=1.0e-10)
    dt = fdt * np.min(dt_raw)
    return max(dt, 1.0e-10)


def _map_terms_nap(grid, l, m, dy):
    """NAP MAP-equivalent terms at (L,M) in 0-based indexing."""
    be = 1.0 / max(grid.yw[l] - grid.ycb[l], 1.0e-12)
    y = m * dy
    al = be * (grid.nxnycb[l] + y * (grid.nxny[l] - grid.nxnycb[l]))
    de = -be * y * 0.0  # XWI(L) currently zero in this implementation.
    return al, be, de


def _interior_step_nap(grid, rho, u, v, p, gamma, dt, p_floor, rho_floor):
    """
    Port of NAP INTER predictor/corrector structure for interior points.
    Returns updated primitive fields.
    """
    lmax, mmax = rho.shape
    dx = (grid.xe - grid.xi) / max(lmax - 1, 1)
    dxr = 1.0 / max(dx, 1.0e-12)
    dy = 1.0 / max(mmax - 1, 1)
    dyr = 1.0 / max(dy, 1.0e-12)

    u_t = u.copy()
    v_t = v.copy()
    p_t = p.copy()
    rho_t = rho.copy()

    # Predictor: backward differencing in computational directions.
    for l in range(1, lmax):
        for m in range(0, mmax - 1):
            al, be, de = _map_terms_nap(grid, l, m, dy)

            ub = u[l, m]
            vb = v[l, m]
            pb = p[l, m]
            rob = max(rho[l, m], rho_floor)
            asb = gamma * pb / rob

            if m == 0:
                dudx = (ub - u[l - 1, m]) * dxr
                dpdx = (pb - p[l - 1, m]) * dxr
                drodx = (rob - rho[l - 1, m]) * dxr
                dvd_y = (4.0 * v[l, 1] - v[l, 2]) * 0.5 * dyr if mmax > 2 else 0.0

                urhs = -ub * dudx - dpdx / rob
                rorhs = -ub * drodx - rob * dudx - (1 + grid.ndim) * rob * be * dvd_y
                prhs = -ub * dpdx + asb * (rorhs + ub * drodx)
                vrhs = 0.0
            else:
                yden = max((m * dy) / max(be, 1.0e-12) + grid.ycb[l], 1.0e-12)
                aterm = rob * vb / yden if grid.ndim == 1 else 0.0
                uvb = ub * al + vb * be + de
                dudx = (ub - u[l - 1, m]) * dxr
                dvdx = (vb - v[l - 1, m]) * dxr
                dpdx = (pb - p[l - 1, m]) * dxr
                drodx = (rob - rho[l - 1, m]) * dxr
                dudy = (ub - u[l, m - 1]) * dyr
                dvdy = (vb - v[l, m - 1]) * dyr
                dpdy = (pb - p[l, m - 1]) * dyr
                drody = (rob - rho[l, m - 1]) * dyr

                urhs = -ub * dudx - uvb * dudy - (dpdx + al * dpdy) / rob
                vrhs = -ub * dvdx - uvb * dvdy - be * dpdy / rob
                rorhs = -ub * drodx - uvb * drody - rob * (dudx + al * dudy + be * dvdy) - aterm
                prhs = -ub * dpdx - uvb * dpdy + asb * (rorhs + ub * drodx + uvb * drody)

            u_t[l, m] = u[l, m] + urhs * dt
            v_t[l, m] = 0.0 if m == 0 else v[l, m] + vrhs * dt
            p_t[l, m] = max(p[l, m] + prhs * dt, p_floor)
            rho_t[l, m] = max(rho[l, m] + rorhs * dt, rho_floor)

    # Corrector: forward differencing and average with predictor.
    u_n = u.copy()
    v_n = v.copy()
    p_n = p.copy()
    rho_n = rho.copy()

    for l in range(1, lmax - 1):
        for m in range(0, mmax - 1):
            al, be, de = _map_terms_nap(grid, l, m, dy)

            ub = u_t[l, m]
            vb = v_t[l, m]
            pb = p_t[l, m]
            rob = max(rho_t[l, m], rho_floor)
            asb = gamma * pb / rob

            if m == 0:
                dudx = (u_t[l + 1, m] - ub) * dxr
                dpdx = (p_t[l + 1, m] - pb) * dxr
                drodx = (rho_t[l + 1, m] - rob) * dxr
                dvd_y = (4.0 * v_t[l, 1] - v_t[l, 2]) * 0.5 * dyr if mmax > 2 else 0.0

                urhs = -ub * dudx - dpdx / rob
                rorhs = -ub * drodx - rob * dudx - (1 + grid.ndim) * rob * be * dvd_y
                prhs = -ub * dpdx + asb * (rorhs + ub * drodx)
                vrhs = 0.0
            else:
                yden = max((m * dy) / max(be, 1.0e-12) + grid.ycb[l], 1.0e-12)
                aterm = rob * vb / yden if grid.ndim == 1 else 0.0
                uvb = ub * al + vb * be + de
                dudx = (u_t[l + 1, m] - ub) * dxr
                dvdx = (v_t[l + 1, m] - vb) * dxr
                dpdx = (p_t[l + 1, m] - pb) * dxr
                drodx = (rho_t[l + 1, m] - rob) * dxr
                dudy = (u_t[l, m + 1] - ub) * dyr if m + 1 < mmax else 0.0
                dvdy = (v_t[l, m + 1] - vb) * dyr if m + 1 < mmax else 0.0
                dpdy = (p_t[l, m + 1] - pb) * dyr if m + 1 < mmax else 0.0
                drody = (rho_t[l, m + 1] - rob) * dyr if m + 1 < mmax else 0.0

                urhs = -ub * dudx - uvb * dudy - (dpdx + al * dpdy) / rob
                vrhs = -ub * dvdx - uvb * dvdy - be * dpdy / rob
                rorhs = -ub * drodx - uvb * drody - rob * (dudx + al * dudy + be * dvdy) - aterm
                prhs = -ub * dpdx - uvb * dpdy + asb * (rorhs + ub * drodx + uvb * drody)

            u_n[l, m] = 0.5 * (u[l, m] + u_t[l, m] + urhs * dt)
            v_n[l, m] = 0.0 if m == 0 else 0.5 * (v[l, m] + v_t[l, m] + vrhs * dt)
            p_n[l, m] = max(0.5 * (p[l, m] + p_t[l, m] + prhs * dt), p_floor)
            rho_n[l, m] = max(0.5 * (rho[l, m] + rho_t[l, m] + rorhs * dt), rho_floor)

    return rho_n, u_n, v_n, p_n


def _solution_from_grid(grid, level=0):
    """Create a Solution view (copied arrays) from a NapGrid time level."""
    return Solution(
        d=grid.ro[:, :, level].copy(),
        u=grid.u[:, :, level].copy(),
        v=grid.v[:, :, level].copy(),
        p=grid.p[:, :, level].copy(),
    )


def _write_solution_to_grid(grid, solution, level):
    """Write a Solution back into a NapGrid time level."""
    grid.ro[:, :, level] = solution.d
    grid.u[:, :, level] = solution.u
    grid.v[:, :, level] = solution.v
    grid.p[:, :, level] = solution.p


def _allocate_history(nmax, lmax, mmax):
    """Allocate storage for all intermediate solution surfaces and diagnostics."""
    shape = (nmax + 1, lmax, mmax)
    return History(
        d=np.full(shape, np.nan, dtype=float),
        u=np.full(shape, np.nan, dtype=float),
        v=np.full(shape, np.nan, dtype=float),
        p=np.full(shape, np.nan, dtype=float),
        res=np.full(nmax, np.nan, dtype=float),
        dt=np.full(nmax, np.nan, dtype=float),
    )


def march_to_steady(
    grid,
    gamma,
    nmax=800,
    tconv=5.0e-3,
    fdt=0.35,
    p_exit=None,
    cav=0.0,
    p_floor=1.0e-6,
    rho_floor=1.0e-8,
    nprint=25,
    use_characteristics=False,
    interior_scheme="nap_inter",
    inlet_profile=None,
):
    """
    Transient MacCormack march on primitive variables to steady state.
    """
    n0, n1 = 0, 1
    solution = _solution_from_grid(grid, level=n0)
    history = _allocate_history(nmax, grid.lmax, grid.mmax)
    history.d[0, :, :] = solution.d
    history.u[0, :, :] = solution.u
    history.v[0, :, :] = solution.v
    history.p[0, :, :] = solution.p

    d = solution.d
    u = solution.u
    v = solution.v
    p = solution.p

    if p_exit is None:
        p_exit = float(np.mean(p[-1, :]))

    inlet_state = {"rho": d[0, :].copy(), "u": u[0, :].copy(), "v": v[0, :].copy(), "p": p[0, :].copy()}
    if inlet_profile is not None:
        inlet_state["rho"] = inlet_profile["rho"].copy()
        inlet_state["u"] = inlet_profile["u"].copy()
        inlet_state["v"] = inlet_profile["v"].copy()
        inlet_state["p"] = inlet_profile["p"].copy()
        d[0, :] = inlet_state["rho"]
        u[0, :] = inlet_state["u"]
        v[0, :] = inlet_state["v"]
        p[0, :] = inlet_state["p"]
    if use_characteristics:
        pt_psf = grid.pt * 144.0
        tt_r = grid.tt + 460.0
        theta_rad = np.deg2rad(grid.theta)
        inlet_state["pt"] = pt_psf
        inlet_state["tt"] = tt_r
        inlet_state["theta"] = theta_rad
        inlet_state["r_eff"] = 53.35
    iterations = 0
    converged = False
    dt_last = float("nan")

    for it in range(1, nmax + 1):
        dt = _compute_cfl_timestep(grid, d, u, v, p, gamma, fdt)
        dt_last = dt
        history.dt[it - 1] = dt
        _apply_boundaries(
            grid,
            d,
            u,
            v,
            p,
            inlet_state,
            p_exit,
            p_floor,
            rho_floor,
            gamma,
            use_characteristics=use_characteristics,
            dt=dt,
        )

        if interior_scheme == "nap_inter":
            d_new, u_new, v_new, p_new = _interior_step_nap(
                grid, d, u, v, p, gamma, dt, p_floor, rho_floor
            )
        else:
            rr, ur, vr, pr = _rhs_nonconservative(grid, d, u, v, p, gamma, scheme="forward", cav=cav)
            d_p = d + dt * rr
            u_p = u + dt * ur
            v_p = v + dt * vr
            p_p = p + dt * pr
            _apply_boundaries(
                grid,
                d_p,
                u_p,
                v_p,
                p_p,
                inlet_state,
                p_exit,
                p_floor,
                rho_floor,
                gamma,
                use_characteristics=use_characteristics,
                dt=dt,
            )

            rr2, ur2, vr2, pr2 = _rhs_nonconservative(grid, d_p, u_p, v_p, p_p, gamma, scheme="backward", cav=cav)
            d_new = 0.5 * (d + d_p + dt * rr2)
            u_new = 0.5 * (u + u_p + dt * ur2)
            v_new = 0.5 * (v + v_p + dt * vr2)
            p_new = 0.5 * (p + p_p + dt * pr2)
        _apply_boundaries(
            grid,
            d_new,
            u_new,
            v_new,
            p_new,
            inlet_state,
            p_exit,
            p_floor,
            rho_floor,
            gamma,
            use_characteristics=use_characteristics,
            dt=dt,
        )

        denom = np.maximum(np.abs(p), p_floor)
        res = float(np.nanmax(np.abs((p_new - p) / denom)))
        history.res[it - 1] = res

        d, u, v, p = d_new, u_new, v_new, p_new
        history.d[it, :, :] = d
        history.u[it, :, :] = u
        history.v[it, :, :] = v
        history.p[it, :, :] = p
        iterations = it
        if not np.isfinite(res):
            print(f"Terminating march at iter {it}: non-finite residual.")
            break

        if it % nprint == 0 or it == 1:
            print(f"ITER {it:5d}  dt={dt:.3e}  res={res:.3e}")

        if res <= tconv:
            converged = True
            print(f"Converged at iter {it} with residual {res:.3e}")
            break

    solution = Solution(d=d, u=u, v=v, p=p)
    _write_solution_to_grid(grid, solution, level=n1)

    return {
        "solution": solution,
        "history": history,
        "residuals": history.res[:iterations].copy(),
        "iterations": iterations,
        "converged": converged,
        "dt_last": dt_last,
    }


def verify_against_table22(grid, csv_path):
    """Compare steady solution plane (level 1) against Table 2.2 CSV data."""
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    if data.size == 0:
        raise ValueError(f"No rows found in verification table: {csv_path}")

    u_ref = np.full((grid.lmax, grid.mmax), np.nan)
    v_ref = np.full((grid.lmax, grid.mmax), np.nan)
    p_ref = np.full((grid.lmax, grid.mmax), np.nan)
    mach_ref = np.full((grid.lmax, grid.mmax), np.nan)

    for row in np.atleast_1d(data):
        l = int(round(float(row["L"]))) - 1
        m = int(round(float(row["M"]))) - 1
        if 0 <= l < grid.lmax and 0 <= m < grid.mmax:
            u_ref[l, m] = float(row["U"])
            v_ref[l, m] = float(row["V"])
            p_ref[l, m] = float(row["P"])  # psia
            mach_ref[l, m] = float(row["Mach"])

    u = grid.u[:, :, 1]
    v = grid.v[:, :, 1]
    p = grid.p[:, :, 1] / 144.0  # psia
    a = np.sqrt(np.maximum(1.4 * grid.p[:, :, 1] / grid.ro[:, :, 1], 1.0e-12))
    mach = np.sqrt(u * u + v * v) / a

    mask = ~np.isnan(u_ref)

    def mae(a_, b_):
        vals = np.abs(a_[mask] - b_[mask])
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return float("nan")
        return float(np.mean(vals))

    def mape(a_, b_, floor=1.0e-8):
        vals_a = a_[mask]
        vals_b = b_[mask]
        finite = np.isfinite(vals_a) & np.isfinite(vals_b)
        if np.sum(finite) == 0:
            return float("nan")
        denom = np.maximum(np.abs(vals_b[finite]), floor)
        return float(np.mean(np.abs((vals_a[finite] - vals_b[finite]) / denom)) * 100.0)

    return {
        "count": int(np.sum(mask)),
        "u_mae_fps": mae(u, u_ref),
        "u_mape_pct": mape(u, u_ref),
        "v_mae_fps": mae(v, v_ref),
        "v_mape_pct": mape(v, v_ref),
        "p_mae_psia": mae(p, p_ref),
        "p_mape_pct": mape(p, p_ref),
        "mach_mae": mae(mach, mach_ref),
        "mach_mape_pct": mape(mach, mach_ref),
    }


def load_table22_inlet_profile(csv_path, l_target=1):
    """Load inlet profile (given L index in 1-based table indexing) from Table 2.2 CSV."""
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    rows = []
    for row in np.atleast_1d(data):
        if int(round(float(row["L"]))) == l_target:
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows found for L={l_target} in {csv_path}")

    rows = sorted(rows, key=lambda r: int(round(float(r["M"]))))
    u = np.array([float(r["U"]) for r in rows], dtype=float)
    v = np.array([float(r["V"]) for r in rows], dtype=float)
    p = np.array([float(r["P"]) for r in rows], dtype=float) * 144.0  # psia -> psf
    rho = np.array([float(r["RHO"]) for r in rows], dtype=float)
    return {"u": u, "v": v, "p": p, "rho": rho}


def build_nap_grid(
    lmax,
    mmax,
    xi,
    ri,
    rt,
    xe,
    rci,
    rct,
    ani_deg=44.88,
    ane_deg=15.0,
    ndim=1,
):
    """Build and return the working structured grid object."""
    grid = StructuredGrid(
        lmax=lmax,
        mmax=mmax,
        xi=xi,
        ri=ri,
        rt=rt,
        xe=xe,
        rci=rci,
        rct=rct,
        ani_deg=ani_deg,
        ane_deg=ane_deg,
        ndim=ndim,
    )
    return grid.build()

def main():
    
    print('NAP (Nozzle Analysis Program) Solver')
    
    grid = StructuredGrid()
    
    sol = Solver(grid)
    
    sol.solve()

    plots(sol)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NAP geometry and solver driver")
    parser.add_argument(
        "--run-clean-port",
        action="store_true",
        help="Also run the experimental clean Fortran-port march (may be unstable).",
    )
    args = parser.parse_args()

    grid1 = StructuredGrid()
    
    grid = build_nap_grid(
        lmax=21,
        mmax=8,
        xi=0.31,
        ri=2.50,
        rt=0.80,
        xe=4.05,
        rci=0.80,
        rct=0.50,
        ani_deg=44.88,
        ane_deg=15.0,
        ndim=1,
    )

    print(f"LMAX={grid.lmax}, MMAX={grid.mmax}, LT={grid.lt}, XT={grid.xt:.4f}")
    print(f"X shape={grid.x.shape}, Y shape={grid.y.shape}")
    print(f"U/V/P/RO shape={grid.u.shape}")
    print(f"Metric shape (J)={grid.jac.shape}")
    print("Boundaries:", ", ".join(grid.boundary_lm.keys()))

    # Example ONEDIM setup using NAP-style English inputs.
    # PT in psia and TT in deg-F are converted internally to psf and deg-R.
    grid.pt[:] = 70.00
    grid.tt[:] = 80.0
    mach_l = initialize_1d_surface(
        grid,
        gamma=1.4,
        rgas=53.35,
        n1d=1,
        level=0,
        unit_system="imperial",
        pt_in_psia=True,
        tt_in_fahrenheit=True,
        gc=32.174,
        verbose=True,
    )
    print(f"ONEDIM centerline Mach: inlet={mach_l[0]:.3f}, throat={mach_l[grid.lt]:.3f}, exit={mach_l[-1]:.3f}")

    plot_path = "artifacts/plots/nap_grid_boundaries.png"
    plot_grid_with_boundaries(grid, plot_path, show=True)
    print(f"Saved plot: {plot_path}")

    solution_plot_path = "artifacts/plots/nap_1d_solution_surface.png"
    plot_1d_solution_surface(
        grid,
        gamma=1.4,
        level=0,
        output_path=solution_plot_path,
        show=True,
        unit_system="imperial",
    )
    print(f"Saved plot: {solution_plot_path}")

    inlet_profile = load_table22_inlet_profile("artifacts/fig8_qc_latest/parsed_table.csv", l_target=1)

    transient = march_to_steady(
        grid,
        gamma=1.4,
        nmax=400,
        tconv=5.0e-3,
        fdt=0.25,
        p_exit=14.7 * 144.0,
        cav=0.0015,
        nprint=20,
        use_characteristics=False,
        inlet_profile=inlet_profile,
    )
    res_hist = transient["history"].res[: transient["iterations"]]
    final_res = float(res_hist[-1]) if transient["iterations"] > 0 else float("nan")
    print(
        "Transient summary: "
        f"iters={transient['iterations']}, converged={transient['converged']}, "
        f"final_res={final_res:.3e}"
    )

    steady_plot_path = "artifacts/plots/nap_steady_solution_surface.png"
    plot_1d_solution_surface(
        grid,
        gamma=1.4,
        level=1,
        output_path=steady_plot_path,
        show=True,
        unit_system="imperial",
    )
    print(f"Saved plot: {steady_plot_path}")

    residual_plot_path = "artifacts/plots/nap_residual_history.png"
    plot_residual_history(res_hist, residual_plot_path, show=True)
    print(f"Saved plot: {residual_plot_path}")

    verify = verify_against_table22(grid, "artifacts/fig8_qc_latest/parsed_table.csv")
    print(
        "Table 2.2 verification: "
        f"N={verify['count']}, "
        f"U MAE={verify['u_mae_fps']:.2f} fps ({verify['u_mape_pct']:.2f}%), "
        f"V MAE={verify['v_mae_fps']:.2f} fps ({verify['v_mape_pct']:.2f}%), "
        f"P MAE={verify['p_mae_psia']:.3f} psia ({verify['p_mape_pct']:.2f}%), "
        f"Mach MAE={verify['mach_mae']:.4f} ({verify['mach_mape_pct']:.2f}%)"
    )

    if args.run_clean_port:
        # Also run the clean Fortran-style port march from nap_fortran_port.py.
        from nap_fortran_port import build_state_from_grid, march_case1, verify_table22_state

        state = build_state_from_grid(grid, gamma=1.4)
        pt_port = np.zeros(state.mmax + 1, dtype=float)
        tt_port = np.zeros(state.mmax + 1, dtype=float)
        th_port = np.zeros(state.mmax + 1, dtype=float)
        pt_port[1:] = grid.pt * 144.0
        tt_port[1:] = grid.tt + 460.0
        th_port[1:] = grid.theta

        res_port = march_case1(
            state,
            pt=pt_port,
            tt=tt_port,
            theta_deg=th_port,
            pe_psf=14.7 * 144.0,
            nmax=200,
            dt=1.0e-5,
            tconv=5.0e-3,
        )
        verify_port = verify_table22_state(state, "artifacts/fig8_qc_latest/parsed_table.csv")
        print(
            "Clean-port verification: "
            f"iters={len(res_port)}, final_res={res_port[-1]:.3e}, "
            f"N={verify_port['count']}, "
            f"U MAE={verify_port['u_mae_fps']:.2f} fps ({verify_port['u_mape_pct']:.2f}%), "
            f"V MAE={verify_port['v_mae_fps']:.2f} fps ({verify_port['v_mape_pct']:.2f}%), "
            f"P MAE={verify_port['p_mae_psia']:.3f} psia ({verify_port['p_mape_pct']:.2f}%), "
            f"Mach MAE={verify_port['mach_mae']:.4f} ({verify_port['mach_mape_pct']:.2f}%)"
        )
