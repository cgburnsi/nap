#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import main


def load_baseline(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})

    by_x = defaultdict(list)
    for row in rows:
        by_x[row["X"]].append(row)

    x_stations = np.array(sorted(by_x.keys()))
    mach_cl = []
    mach_w = []
    p_cl = []
    p_w = []
    for x in x_stations:
        station = by_x[x]
        center = min(station, key=lambda r: r["Y"])
        wall = max(station, key=lambda r: r["Y"])
        mach_cl.append(center["Mach"])
        mach_w.append(wall["Mach"])
        p_cl.append(center["P"])
        p_w.append(wall["P"])

    return {
        "x": x_stations,
        "mach_cl": np.array(mach_cl),
        "mach_w": np.array(mach_w),
        "p_cl": np.array(p_cl),
        "p_w": np.array(p_w),
    }


def rmse(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def run_case(args):
    u, v, p, d, params, geometry = main.initalize_simulation()
    params["nmax"] = args.nmax
    params["fdt"] = args.fdt
    params["iav"] = args.iav
    params["smp"] = args.smp
    params["exit_bc_mode"] = args.exit_mode
    params["pe"] = args.pe

    u, v, p, d, params, geometry = main.run_simulation(u, v, p, d, params, geometry)

    n_idx = (params.get("iteration_count", 0) + 1) % 2
    gamma = params["gamma"]
    g = params["g"]
    a = np.sqrt(gamma * p[:, :, n_idx] * 144.0 * g / d[:, :, n_idx])
    mach = np.sqrt(u[:, :, n_idx] ** 2 + v[:, :, n_idx] ** 2) / a

    x = geometry["xw"]
    run = {
        "x": x,
        "mach_cl": mach[:, 0],
        "mach_w": mach[:, -1],
        "p_cl": p[:, 0, n_idx],
        "p_w": p[:, -1, n_idx],
    }
    return run, params


def compare(run, base):
    x_ref = base["x"]
    run_mach_cl = np.interp(x_ref, run["x"], run["mach_cl"])
    run_mach_w = np.interp(x_ref, run["x"], run["mach_w"])
    run_p_cl = np.interp(x_ref, run["x"], run["p_cl"])
    run_p_w = np.interp(x_ref, run["x"], run["p_w"])

    metrics = {
        "rmse_mach_centerline": rmse(run_mach_cl, base["mach_cl"]),
        "rmse_mach_wall": rmse(run_mach_w, base["mach_w"]),
        "rmse_p_centerline_psia": rmse(run_p_cl, base["p_cl"]),
        "rmse_p_wall_psia": rmse(run_p_w, base["p_w"]),
        "exit_mach_centerline": float(run_mach_cl[-1]),
        "exit_mach_wall": float(run_mach_w[-1]),
        "exit_p_centerline_psia": float(run_p_cl[-1]),
        "exit_p_wall_psia": float(run_p_w[-1]),
        "baseline_exit_mach_centerline": float(base["mach_cl"][-1]),
        "baseline_exit_mach_wall": float(base["mach_w"][-1]),
        "baseline_exit_p_centerline_psia": float(base["p_cl"][-1]),
        "baseline_exit_p_wall_psia": float(base["p_w"][-1]),
    }
    return metrics


def main_cli():
    parser = argparse.ArgumentParser(description="Run NAP solver and compare against baseline.")
    parser.add_argument("--baseline", default="docs/verification/baselines/cd_nozzle_converged_solution.csv")
    parser.add_argument("--nmax", type=int, default=1200)
    parser.add_argument("--fdt", type=float, default=0.15)
    parser.add_argument("--smp", type=float, default=0.03)
    parser.add_argument("--iav", type=int, default=0, choices=[0, 1])
    parser.add_argument("--exit-mode", default="pressure", choices=["pressure", "extrapolate"])
    parser.add_argument("--pe", type=float, default=5.3)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    baseline = load_baseline(args.baseline)
    run, params = run_case(args)
    metrics = compare(run, baseline)
    metrics["iterations"] = int(params.get("iteration_count", -1))
    metrics["fdt"] = args.fdt
    metrics["smp"] = args.smp
    metrics["iav"] = args.iav
    metrics["exit_mode"] = args.exit_mode
    metrics["pe"] = args.pe

    if args.json:
        print(json.dumps(metrics, indent=2, sort_keys=True))
        return

    print("Regression Metrics")
    print(f"iterations: {metrics['iterations']}")
    print(f"rmse Mach centerline: {metrics['rmse_mach_centerline']:.4f}")
    print(f"rmse Mach wall      : {metrics['rmse_mach_wall']:.4f}")
    print(f"rmse P centerline   : {metrics['rmse_p_centerline_psia']:.4f} psia")
    print(f"rmse P wall         : {metrics['rmse_p_wall_psia']:.4f} psia")
    print(
        "exit Mach (run/base) centerline: "
        f"{metrics['exit_mach_centerline']:.4f} / {metrics['baseline_exit_mach_centerline']:.4f}"
    )
    print(
        "exit Mach (run/base) wall      : "
        f"{metrics['exit_mach_wall']:.4f} / {metrics['baseline_exit_mach_wall']:.4f}"
    )
    print(
        "exit P (run/base) centerline   : "
        f"{metrics['exit_p_centerline_psia']:.4f} / {metrics['baseline_exit_p_centerline_psia']:.4f} psia"
    )
    print(
        "exit P (run/base) wall         : "
        f"{metrics['exit_p_wall_psia']:.4f} / {metrics['baseline_exit_p_wall_psia']:.4f} psia"
    )


if __name__ == "__main__":
    main_cli()
