import numpy as np


def default_case1_config():
    # NAP-style baseline case scaffold: 21 x 8 converging-diverging nozzle.
    return {
        "name": "case1_cd_nozzle",
        "lmax": 21,
        "mmax": 8,
        "x_inlet": 0.31,
        "x_exit": 4.05,
        "r_inlet": 2.50,
        "r_throat": 1.45,
        "r_exit": 3.20,
        "throat_x": 1.25,
        "gamma": 1.4,
        "rgas": 53.35,
        "pt_inlet": 70.0,
        "tt_inlet": 80.0,
        "pe_exit": 69.2,
        "cfl": 0.6,
        "nmax": 200,
        "fdt": 5,
        "tconv": 1.0e-5,
        "artificial_viscosity": 0.02,
        "onedim_supersonic_after_throat": False,
        "baseline_csv": "docs/verification/baselines/cd_nozzle_converged_solution.csv",
    }


def build_geometry(cfg):
    x = np.linspace(cfg["x_inlet"], cfg["x_exit"], cfg["lmax"])

    # Simple piecewise-linear wall profile placeholder for case-1 style CD nozzle.
    r_wall = np.empty_like(x)
    left = x <= cfg["throat_x"]
    r_wall[left] = np.interp(
        x[left], [cfg["x_inlet"], cfg["throat_x"]], [cfg["r_inlet"], cfg["r_throat"]]
    )
    r_wall[~left] = np.interp(
        x[~left], [cfg["throat_x"], cfg["x_exit"]], [cfg["r_throat"], cfg["r_exit"]]
    )

    eta = np.linspace(0.0, 1.0, cfg["mmax"])
    X = np.repeat(x[:, None], cfg["mmax"], axis=1)
    Y = np.repeat((eta[None, :]), cfg["lmax"], axis=0) * r_wall[:, None]
    return {"x_line": x, "r_wall": r_wall, "X": X, "Y": Y}


def area_mach_ratio(gamma, mach):
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    term = (2.0 / gp1) * (1.0 + 0.5 * gm1 * mach * mach)
    exp = gp1 / (2.0 * gm1)
    return (1.0 / mach) * (term**exp)


def solve_mach_from_area_ratio(gamma, area_ratio, branch):
    area_ratio = max(area_ratio, 1.0)

    if abs(area_ratio - 1.0) < 1.0e-10:
        return 1.0

    if branch == "subsonic":
        lo, hi = 1.0e-6, 1.0 - 1.0e-8
        guess = 0.2
    else:
        lo, hi = 1.0 + 1.0e-8, 8.0
        guess = 2.0

    target = area_ratio
    m = guess

    for _ in range(40):
        f = area_mach_ratio(gamma, m) - target
        if abs(f) < 1.0e-10:
            return m

        h = max(1.0e-6, 1.0e-4 * m)
        fp = area_mach_ratio(gamma, m + h) - target
        fm = area_mach_ratio(gamma, m - h) - target
        dfdm = (fp - fm) / (2.0 * h)
        if abs(dfdm) < 1.0e-12:
            break

        m_new = m - f / dfdm
        if m_new <= lo or m_new >= hi:
            break
        m = m_new

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fmid = area_mach_ratio(gamma, mid) - target
        if abs(fmid) < 1.0e-10 or (hi - lo) < 1.0e-10:
            return mid

        flo = area_mach_ratio(gamma, lo) - target
        if np.sign(fmid) == np.sign(flo):
            lo = mid
        else:
            hi = mid

    return 0.5 * (lo + hi)


def onedim_isentropic_profile(cfg, geom):
    gamma = cfg["gamma"]
    rgas = cfg["rgas"]
    pt = cfg["pt_inlet"]
    tt = cfg["tt_inlet"]

    r = geom["r_wall"]
    throat_idx = int(np.argmin(r))
    area = np.pi * r * r
    area_star = np.min(area)
    area_ratio = area / area_star

    mach = np.empty_like(r)
    for i in range(r.size):
        if i < throat_idx:
            branch = "subsonic"
        elif i > throat_idx and cfg["onedim_supersonic_after_throat"]:
            branch = "supersonic"
        else:
            branch = "subsonic"
        mach[i] = solve_mach_from_area_ratio(gamma, area_ratio[i], branch)

    t = tt / (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
    p = pt / (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (gamma / (gamma - 1.0))
    rho = p / (rgas * t)
    a = np.sqrt(gamma * p / np.maximum(rho, 1.0e-12))
    u = mach * a
    return u, p, rho


def initialize_state(cfg, geom):
    shape = geom["X"].shape
    u_1d, p_1d, rho_1d = onedim_isentropic_profile(cfg, geom)

    u = np.repeat(u_1d[:, None], shape[1], axis=1)
    v = np.zeros(shape)
    p = np.repeat(p_1d[:, None], shape[1], axis=1)
    rho = np.repeat(rho_1d[:, None], shape[1], axis=1)

    return {"u": u, "v": v, "p": p, "rho": rho, "time": 0.0}


def update_timestep(cfg, geom, state):
    x = geom["x_line"]
    y = geom["Y"]
    dx = np.min(np.diff(x))
    dy = np.min(np.diff(y, axis=1))

    a = np.sqrt(cfg["gamma"] * state["p"] / np.maximum(state["rho"], 1.0e-12))
    vel = np.sqrt(state["u"] ** 2 + state["v"] ** 2)
    spectral = np.max((vel + a) * np.sqrt((1.0 / dx**2) + (1.0 / dy**2)))
    return cfg["cfl"] / max(spectral, 1.0e-12)


def apply_boundary_conditions(cfg, geom, state):
    # Inlet: keep simple static values for now; replace with characteristic BCs.
    state["p"][0, :] = cfg["pt_inlet"] / 1.01
    state["rho"][0, :] = 0.35
    state["u"][0, :] = 145.0
    state["v"][0, :] = 0.0

    # Exit: pressure condition for subsonic outlet.
    state["p"][-1, :] = cfg["pe_exit"]

    # Axis symmetry and wall slip placeholders.
    state["v"][:, 0] = 0.0
    state["v"][:, -1] = 0.0


def maccormack_predictor(cfg, geom, state, dt):
    pred = {k: np.array(v, copy=True) for k, v in state.items() if isinstance(v, np.ndarray)}
    pred["time"] = state["time"] + dt

    # Placeholder predictor with low-order smoothing to keep the scaffold numerically calm.
    for key in ("u", "v", "p", "rho"):
        arr = pred[key]
        arr[1:-1, 1:-1] = 0.5 * arr[1:-1, 1:-1] + 0.125 * (
            arr[:-2, 1:-1] + arr[2:, 1:-1] + arr[1:-1, :-2] + arr[1:-1, 2:]
        )
    return pred


def maccormack_corrector(cfg, geom, state, pred, dt):
    corr = {k: np.array(v, copy=True) for k, v in pred.items() if isinstance(v, np.ndarray)}
    corr["time"] = pred["time"]

    # Explicit artificial viscosity placeholder (replace with NAP von Neumann-Richtmyer form).
    nu = cfg["artificial_viscosity"]
    for key in ("u", "v", "p", "rho"):
        arr = corr[key]
        lap = (
            arr[:-2, 1:-1] + arr[2:, 1:-1] + arr[1:-1, :-2] + arr[1:-1, 2:] - 4.0 * arr[1:-1, 1:-1]
        )
        arr[1:-1, 1:-1] = arr[1:-1, 1:-1] + nu * lap
    return corr


def compute_residual(prev, curr):
    p0 = np.maximum(np.abs(prev["p"]), 1.0e-8)
    return np.max(np.abs(curr["p"] - prev["p"]) / p0)


def march_solver(cfg, geom, state0):
    state = state0
    dt = update_timestep(cfg, geom, state)
    history = []

    for n in range(1, cfg["nmax"] + 1):
        if n % cfg["fdt"] == 0:
            dt = update_timestep(cfg, geom, state)

        pred = maccormack_predictor(cfg, geom, state, dt)
        apply_boundary_conditions(cfg, geom, pred)

        curr = maccormack_corrector(cfg, geom, state, pred, dt)
        apply_boundary_conditions(cfg, geom, curr)

        residual = compute_residual(state, curr)
        history.append((n, curr["time"], dt, residual))
        state = curr

        if residual < cfg["tconv"]:
            break

    return state, history


def flatten_solution(cfg, geom, state):
    gamma = cfg["gamma"]
    rgas = cfg["rgas"]

    u = state["u"]
    v = state["v"]
    p = state["p"]
    rho = state["rho"]

    q = np.sqrt(u**2 + v**2)
    a = np.sqrt(gamma * p / np.maximum(rho, 1.0e-12))
    mach = q / np.maximum(a, 1.0e-12)
    temp = p / np.maximum(rho * rgas, 1.0e-12)

    data = []
    for l in range(cfg["lmax"]):
        for m in range(cfg["mmax"]):
            data.append(
                [
                    float(l + 1),
                    float(m + 1),
                    geom["X"][l, m],
                    geom["Y"][l, m],
                    u[l, m],
                    v[l, m],
                    p[l, m],
                    rho[l, m],
                    q[l, m],
                    mach[l, m],
                    temp[l, m],
                ]
            )
    return np.array(data)


def load_baseline_table(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def compare_to_baseline(solution_table, baseline):
    cols = {"U": 4, "V": 5, "P": 6, "RHO": 7, "Mach": 9, "T": 10}
    metrics = {}
    for key, idx in cols.items():
        sim = solution_table[:, idx]
        ref = np.asarray(baseline[key], dtype=float)
        rmse = np.sqrt(np.mean((sim - ref) ** 2))
        rel = rmse / max(np.sqrt(np.mean(ref**2)), 1.0e-12)
        metrics[key] = {"rmse": rmse, "relative_rmse": rel}
    return metrics


def print_summary(history, metrics):
    if history:
        n, t, dt, r = history[-1]
        print(f"Final step: n={n}, time={t:.6f}, dt={dt:.6e}, residual={r:.3e}")
    else:
        print("Solver did not run any steps.")

    print("\nBaseline comparison (RMSE / relative RMSE):")
    for key in ("U", "V", "P", "RHO", "Mach", "T"):
        m = metrics[key]
        print(f"  {key:>4}: {m['rmse']:.6e} / {m['relative_rmse']:.6e}")

    print("\nNext implementation targets:")
    print("  1) Replace predictor/corrector placeholders with NAP Eqns (interior update).")
    print("  2) Implement characteristic BCs (inlet, exit, wall, jet boundary).")
    print("  3) Replace viscosity placeholder with NAP shock model.")
    print("  4) Track mass flow/thrust and tune to case-1 baseline.")


def main():
    cfg = default_case1_config()
    geom = build_geometry(cfg)
    state0 = initialize_state(cfg, geom)
    apply_boundary_conditions(cfg, geom, state0)
    state, history = march_solver(cfg, geom, state0)

    solution_table = flatten_solution(cfg, geom, state)
    baseline = load_baseline_table(cfg["baseline_csv"])
    metrics = compare_to_baseline(solution_table, baseline)
    print_summary(history, metrics)


if __name__ == "__main__":
    main()
