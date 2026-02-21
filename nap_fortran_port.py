"""Clean NAP Fortran-style port scaffold (1-based indexing)."""

from dataclasses import dataclass

import numpy as np


@dataclass
class NAPPortState:
    """State container that mirrors Fortran common-block array conventions."""

    lmax: int
    mmax: int
    gamma: float
    p_floor: float
    rho_floor: float
    ndimsym: int
    # Geometry-like arrays (1..LMAX used)
    xw: np.ndarray
    yw: np.ndarray
    ycb: np.ndarray
    nxny: np.ndarray
    nxnycb: np.ndarray
    xwi: np.ndarray
    # Flow arrays (1..LMAX, 1..MMAX, 1..2 used)
    u: np.ndarray
    v: np.ndarray
    p: np.ndarray
    ro: np.ndarray


def build_state_from_grid(grid, gamma=1.4, p_floor=1.0e-6, rho_floor=1.0e-8):
    """Build 1-based state arrays from the existing Python grid data."""
    lmax = int(grid.lmax)
    mmax = int(grid.mmax)

    xw = np.zeros(lmax + 1, dtype=float)
    yw = np.zeros(lmax + 1, dtype=float)
    ycb = np.zeros(lmax + 1, dtype=float)
    nxny = np.zeros(lmax + 1, dtype=float)
    nxnycb = np.zeros(lmax + 1, dtype=float)
    xwi = np.zeros(lmax + 1, dtype=float)

    # 1-based fill.
    xw[1:] = grid.xw
    yw[1:] = grid.yw
    ycb[1:] = grid.ycb
    nxny[1:] = grid.nxny
    nxnycb[1:] = grid.nxnycb

    u = np.zeros((lmax + 1, mmax + 1, 3), dtype=float)
    v = np.zeros((lmax + 1, mmax + 1, 3), dtype=float)
    p = np.zeros((lmax + 1, mmax + 1, 3), dtype=float)
    ro = np.zeros((lmax + 1, mmax + 1, 3), dtype=float)

    # Load from existing levels 0 and 1 into Fortran levels 1 and 2.
    u[1:, 1:, 1] = grid.u[:, :, 0]
    v[1:, 1:, 1] = grid.v[:, :, 0]
    p[1:, 1:, 1] = grid.p[:, :, 0]
    ro[1:, 1:, 1] = grid.ro[:, :, 0]

    u[1:, 1:, 2] = grid.u[:, :, 1]
    v[1:, 1:, 2] = grid.v[:, :, 1]
    p[1:, 1:, 2] = grid.p[:, :, 1]
    ro[1:, 1:, 2] = grid.ro[:, :, 1]

    return NAPPortState(
        lmax=lmax,
        mmax=mmax,
        gamma=gamma,
        p_floor=p_floor,
        rho_floor=rho_floor,
        ndimsym=int(grid.ndim),
        xw=xw,
        yw=yw,
        ycb=ycb,
        nxny=nxny,
        nxnycb=nxnycb,
        xwi=xwi,
        u=u,
        v=v,
        p=p,
        ro=ro,
    )


def map_terms(state: NAPPortState, ip: int, l: int, m: int, dy: float, ld1: int):
    """Direct port of MAP(IP,L,M,AL,BE,DE,LD1,AL1,BE1,DE1)."""
    be = 1.0 / max(state.yw[l] - state.ycb[l], 1.0e-12)
    if ip == 0:
        return 0.0, be, 0.0, 0.0, 0.0, 0.0

    y = (m - 1) * dy
    al = be * (state.nxnycb[l] + y * (state.nxny[l] - state.nxnycb[l]))
    de = -be * y * state.xwi[l]
    if ip == 1:
        return al, be, de, 0.0, 0.0, 0.0

    be1 = 1.0 / max(state.yw[ld1] - state.ycb[ld1], 1.0e-12)
    al1 = be1 * (state.nxnycb[ld1] + y * (state.nxny[ld1] - state.nxnycb[ld1]))
    de1 = -be1 * y * state.xwi[ld1]
    return al, be, de, al1, be1, de1


def inter_step(state: NAPPortState, n1=1, n3=2, ichar=1):
    """
    Direct-structure port of INTER predictor/corrector for interior points.
    Set ichar=1 for predictor pass and ichar=2 for corrector pass.
    """
    lmax = state.lmax
    mmax = state.mmax
    dx = (state.xw[lmax] - state.xw[1]) / max(lmax - 1, 1)
    dy = 1.0 / max(mmax - 1, 1)
    dxr = 1.0 / max(dx, 1.0e-12)
    dyr = 1.0 / max(dy, 1.0e-12)

    mdum = 2 if np.any(np.abs(state.ycb[1:]) > 1.0e-12) else 1
    l_stop = lmax if ichar == 1 else lmax - 1

    for l in range(2, l_stop + 1):
        for m in range(mdum, mmax):
            ld1 = l - 1 if ichar == 1 else l + 1
            if ld1 < 1 or ld1 > lmax:
                continue

            al, be, de, _, _, _ = map_terms(state, 1, l, m, dy, ld1)
            ub = state.u[l, m, n1 if ichar == 1 else n3]
            vb = state.v[l, m, n1 if ichar == 1 else n3]
            pb = max(state.p[l, m, n1 if ichar == 1 else n3], state.p_floor)
            rob = max(state.ro[l, m, n1 if ichar == 1 else n3], state.rho_floor)
            asb = state.gamma * pb / rob

            if m == 1:
                dudx = (ub - state.u[ld1, m, n1]) * dxr if ichar == 1 else (state.u[ld1, m, n3] - ub) * dxr
                dpdx = (pb - state.p[ld1, m, n1]) * dxr if ichar == 1 else (state.p[ld1, m, n3] - pb) * dxr
                drodx = (rob - state.ro[ld1, m, n1]) * dxr if ichar == 1 else (state.ro[ld1, m, n3] - rob) * dxr
                dvd_y = (4.0 * state.v[l, 2, n1 if ichar == 1 else n3] - state.v[l, 3, n1 if ichar == 1 else n3])
                dvd_y *= 0.5 * dyr

                urhs = -ub * dudx - dpdx / rob
                rorhs = -ub * drodx - rob * dudx - (1 + state.ndimsym) * rob * be * dvd_y
                prhs = -ub * dpdx + asb * (rorhs + ub * drodx)
                vrhs = 0.0
            else:
                aterm = 0.0
                if state.ndimsym == 1:
                    aterm = rob * vb / max((m - 1) * dy / max(be, 1.0e-12) + state.ycb[l], 1.0e-12)
                uvb = ub * al + vb * be + de

                if ichar == 1:
                    dudx = (ub - state.u[l - 1, m, n1]) * dxr
                    dvdx = (vb - state.v[l - 1, m, n1]) * dxr
                    dpdx = (pb - state.p[l - 1, m, n1]) * dxr
                    drodx = (rob - state.ro[l - 1, m, n1]) * dxr
                    dudy = (ub - state.u[l, m - 1, n1]) * dyr
                    dvdy = (vb - state.v[l, m - 1, n1]) * dyr
                    dpdy = (pb - state.p[l, m - 1, n1]) * dyr
                    drody = (rob - state.ro[l, m - 1, n1]) * dyr
                else:
                    dudx = (state.u[l + 1, m, n3] - ub) * dxr
                    dvdx = (state.v[l + 1, m, n3] - vb) * dxr
                    dpdx = (state.p[l + 1, m, n3] - pb) * dxr
                    drodx = (state.ro[l + 1, m, n3] - rob) * dxr
                    dudy = (state.u[l, m + 1, n3] - ub) * dyr
                    dvdy = (state.v[l, m + 1, n3] - vb) * dyr
                    dpdy = (state.p[l, m + 1, n3] - pb) * dyr
                    drody = (state.ro[l, m + 1, n3] - rob) * dyr

                urhs = -ub * dudx - uvb * dudy - (dpdx + al * dpdy) / rob
                vrhs = -ub * dvdx - uvb * dvdy - be * dpdy / rob
                rorhs = -ub * drodx - uvb * drody - rob * (dudx + al * dudy + be * dvdy) - aterm
                prhs = -ub * dpdx - uvb * dpdy + asb * (rorhs + ub * drodx + uvb * drody)

            if ichar == 1:
                state.u[l, m, n3] = state.u[l, m, n1] + urhs * 1.0e-5
                state.v[l, m, n3] = 0.0 if m == 1 else state.v[l, m, n1] + vrhs * 1.0e-5
                state.p[l, m, n3] = max(state.p[l, m, n1] + prhs * 1.0e-5, state.p_floor)
                state.ro[l, m, n3] = max(state.ro[l, m, n1] + rorhs * 1.0e-5, state.rho_floor)
            else:
                state.u[l, m, n3] = 0.5 * (state.u[l, m, n1] + state.u[l, m, n3] + urhs * 1.0e-5)
                state.v[l, m, n3] = 0.0 if m == 1 else 0.5 * (state.v[l, m, n1] + state.v[l, m, n3] + vrhs * 1.0e-5)
                state.p[l, m, n3] = max(
                    0.5 * (state.p[l, m, n1] + state.p[l, m, n3] + prhs * 1.0e-5),
                    state.p_floor,
                )
                state.ro[l, m, n3] = max(
                    0.5 * (state.ro[l, m, n1] + state.ro[l, m, n3] + rorhs * 1.0e-5),
                    state.rho_floor,
                )


def step_inter(state: NAPPortState):
    """One INTER predictor+corrector step."""
    inter_step(state, n1=1, n3=2, ichar=1)
    inter_step(state, n1=1, n3=2, ichar=2)
    # Rotate levels (Fortran MAIN toggles N1/N3).
    state.u[:, :, 1] = state.u[:, :, 2]
    state.v[:, :, 1] = state.v[:, :, 2]
    state.p[:, :, 1] = state.p[:, :, 2]
    state.ro[:, :, 1] = state.ro[:, :, 2]


def wall_step(state: NAPPortState, n=1, dt=1.0e-5):
    """WALL port scaffold: enforce tangency and update wall row from interior row."""
    m_w = state.mmax
    m_i = state.mmax - 1
    if m_i < 1:
        return

    for l in range(1, state.lmax + 1):
        u_i = state.u[l, m_i, n]
        v_w = -u_i * state.nxny[l]
        state.u[l, m_w, n] = u_i
        state.v[l, m_w, n] = v_w
        state.p[l, m_w, n] = max(state.p[l, m_i, n], state.p_floor)
        state.ro[l, m_w, n] = max(state.ro[l, m_i, n], state.rho_floor)


def inlet_step(state: NAPPortState, pt, tt, theta_deg, n=1):
    """INLET port scaffold: subsonic inlet from PT/TT/theta + outgoing characteristic."""
    l = 1
    l2 = min(2, state.lmax)
    for m in range(1, state.mmax + 1):
        p2 = max(state.p[l2, m, n], state.p_floor)
        ro2 = max(state.ro[l2, m, n], state.rho_floor)
        u2 = state.u[l2, m, n]
        a2 = np.sqrt(max(state.gamma * p2 / ro2, 1.0e-12))
        j_minus = u2 - 2.0 * a2 / (state.gamma - 1.0)

        ttot = float(tt[m])
        ptot = float(pt[m])
        th = np.deg2rad(float(theta_deg[m]))
        gam2 = 0.5 * (state.gamma - 1.0)
        gamp = state.gamma / (state.gamma - 1.0)
        rgas = 53.35

        mach = 0.15
        for _ in range(25):
            dem = 1.0 + gam2 * mach * mach
            t = ttot / dem
            a = np.sqrt(max(state.gamma * rgas * t, 1.0e-12))
            q = mach * a
            f = q * np.cos(th) - 2.0 * a / (state.gamma - 1.0) - j_minus
            dm = 1.0e-5
            mp = mach + dm
            demp = 1.0 + gam2 * mp * mp
            tp = ttot / demp
            ap = np.sqrt(max(state.gamma * rgas * tp, 1.0e-12))
            qp = mp * ap
            fp = qp * np.cos(th) - 2.0 * ap / (state.gamma - 1.0) - j_minus
            df = (fp - f) / dm
            if abs(df) < 1.0e-12:
                break
            mach = float(np.clip(mach - f / df, 1.0e-4, 0.99))

        dem = 1.0 + gam2 * mach * mach
        t = ttot / dem
        p = ptot / (dem**gamp)
        ro = max(p / max(rgas * t, 1.0e-12), state.rho_floor)
        a = np.sqrt(max(state.gamma * p / ro, 1.0e-12))
        q = mach * a
        state.u[l, m, n] = q * np.cos(th)
        state.v[l, m, n] = -q * np.sin(th)
        state.p[l, m, n] = max(p, state.p_floor)
        state.ro[l, m, n] = ro


def exitt_step(state: NAPPortState, pe_psf, n=1):
    """EXITT port scaffold: subsonic fixed static pressure with outgoing characteristic."""
    l = state.lmax
    lm1 = max(1, l - 1)
    pe = max(pe_psf, state.p_floor)

    for m in range(1, state.mmax + 1):
        p2 = max(state.p[lm1, m, n], state.p_floor)
        ro2 = max(state.ro[lm1, m, n], state.rho_floor)
        u2 = state.u[lm1, m, n]
        a2 = np.sqrt(max(state.gamma * p2 / ro2, 1.0e-12))
        j_plus = u2 + 2.0 * a2 / (state.gamma - 1.0)
        s = p2 / max(ro2**state.gamma, 1.0e-12)
        ro = max((pe / max(s, 1.0e-12)) ** (1.0 / state.gamma), state.rho_floor)
        a = np.sqrt(max(state.gamma * pe / ro, 1.0e-12))
        u = j_plus - 2.0 * a / (state.gamma - 1.0)
        state.u[l, m, n] = u
        state.v[l, m, n] = state.v[lm1, m, n]
        state.p[l, m, n] = pe
        state.ro[l, m, n] = ro


def march_case1(state: NAPPortState, pt, tt, theta_deg, pe_psf, nmax=200, dt=1.0e-5, tconv=5.0e-3):
    """
    March driver for clean module:
    order follows Fortran workflow approximation: INTER -> WALL -> INLET -> EXITT.
    """
    residuals = []
    for _ in range(nmax):
        p_prev = state.p[:, :, 1].copy()
        inter_step(state, n1=1, n3=2, ichar=1)
        inter_step(state, n1=1, n3=2, ichar=2)
        wall_step(state, n=2, dt=dt)
        inlet_step(state, pt=pt, tt=tt, theta_deg=theta_deg, n=2)
        exitt_step(state, pe_psf=pe_psf, n=2)

        state.u[:, :, 1] = state.u[:, :, 2]
        state.v[:, :, 1] = state.v[:, :, 2]
        state.p[:, :, 1] = state.p[:, :, 2]
        state.ro[:, :, 1] = state.ro[:, :, 2]

        denom = np.maximum(np.abs(p_prev), state.p_floor)
        res = float(np.max(np.abs((state.p[:, :, 1] - p_prev) / denom)))
        residuals.append(res)
        if res <= tconv:
            break
    return residuals


def verify_table22_state(state: NAPPortState, csv_path: str):
    """Compare current state (level 1) against Table 2.2 style CSV."""
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    if data.size == 0:
        raise ValueError(f"No rows found in verification table: {csv_path}")

    u_ref = np.full((state.lmax + 1, state.mmax + 1), np.nan)
    v_ref = np.full((state.lmax + 1, state.mmax + 1), np.nan)
    p_ref = np.full((state.lmax + 1, state.mmax + 1), np.nan)
    mach_ref = np.full((state.lmax + 1, state.mmax + 1), np.nan)

    for row in np.atleast_1d(data):
        l = int(round(float(row["L"])))
        m = int(round(float(row["M"])))
        if 1 <= l <= state.lmax and 1 <= m <= state.mmax:
            u_ref[l, m] = float(row["U"])
            v_ref[l, m] = float(row["V"])
            p_ref[l, m] = float(row["P"])  # psia
            mach_ref[l, m] = float(row["Mach"])

    u = state.u[:, :, 1]
    v = state.v[:, :, 1]
    p_psia = state.p[:, :, 1] / 144.0
    a = np.sqrt(np.maximum(state.gamma * state.p[:, :, 1] / np.maximum(state.ro[:, :, 1], state.rho_floor), 1.0e-12))
    mach = np.sqrt(u * u + v * v) / a

    mask = ~np.isnan(u_ref)
    if np.sum(mask) == 0:
        raise ValueError("No overlapping (L,M) points between state and table.")

    def mae(a_, b_):
        vals = np.abs(a_[mask] - b_[mask])
        vals = vals[np.isfinite(vals)]
        return float(np.mean(vals)) if vals.size else float("nan")

    def mape(a_, b_, floor=1.0e-8):
        av = a_[mask]
        bv = b_[mask]
        finite = np.isfinite(av) & np.isfinite(bv)
        if not np.any(finite):
            return float("nan")
        denom = np.maximum(np.abs(bv[finite]), floor)
        return float(np.mean(np.abs((av[finite] - bv[finite]) / denom)) * 100.0)

    return {
        "count": int(np.sum(mask)),
        "u_mae_fps": mae(u, u_ref),
        "u_mape_pct": mape(u, u_ref),
        "v_mae_fps": mae(v, v_ref),
        "v_mape_pct": mape(v, v_ref),
        "p_mae_psia": mae(p_psia, p_ref),
        "p_mape_pct": mape(p_psia, p_ref),
        "mach_mae": mae(mach, mach_ref),
        "mach_mape_pct": mape(mach, mach_ref),
    }
