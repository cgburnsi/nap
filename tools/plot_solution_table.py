#!/usr/bin/env python3
"""Parse OCR-normalized LaTeX longtable solution data and generate QA plots.

Usage example:
  python3 tools/plot_solution_table.py \
    --input docs/conversion/fig8/fig8_surface301_longtable.tex \
    --expected-l 21 --expected-m 8 \
    --field Mach --outdir artifacts/fig8_qc
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

COLUMNS = ["L", "M", "X", "Y", "U", "V", "P", "RHO", "Q", "Mach", "T"]


@dataclass
class ParseIssue:
    line_no: int
    reason: str
    raw: str


def _clean_token(token: str) -> str:
    t = token.strip()
    t = t.replace("−", "-")
    t = t.replace("\\", "")
    t = t.replace("\u00a0", "")
    t = t.strip()
    # OCR occasionally inserts internal spaces in numeric literals.
    t = re.sub(r"\s+", "", t)
    return t


def _merge_split_decimal(tokens: List[str]) -> List[str]:
    """Merge obvious OCR splits like '0.' '&' '187128' -> '0.187128'."""
    out: List[str] = []
    i = 0
    while i < len(tokens):
        cur = tokens[i]
        nxt = tokens[i + 1] if i + 1 < len(tokens) else None
        if nxt is not None and re.fullmatch(r"[-+]?\d+\.", cur) and re.fullmatch(r"\d+", nxt):
            out.append(cur + nxt)
            i += 2
            continue
        out.append(cur)
        i += 1
    return out


def _to_float(token: str) -> float:
    if token in {"", ".", "-", "+"}:
        raise ValueError(f"empty numeric token: {token!r}")
    return float(token)


def parse_longtable(path: Path) -> tuple[np.ndarray, List[ParseIssue]]:
    rows: List[List[float]] = []
    issues: List[ParseIssue] = []

    for line_no, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = raw.strip()
        if "&" not in stripped:
            continue
        if any(k in stripped for k in ("\\hline", "\\caption", "\\end", "\\begin", "\\textbf", "\\endhead", "\\endfirsthead")):
            continue
        if not re.match(r"^\s*\d+\s*&", stripped):
            continue

        parts = [_clean_token(p) for p in stripped.split("&")]
        parts = _merge_split_decimal(parts)

        # Drop trailing empty fragments caused by OCR/backslash debris.
        while parts and parts[-1] == "":
            parts.pop()

        if len(parts) != len(COLUMNS):
            # Fallback: try extracting numeric literals from the raw line.
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", stripped)
            if len(nums) == len(COLUMNS):
                parts = nums
            else:
                issues.append(ParseIssue(line_no, f"column_count={len(parts)}", raw))
                continue

        try:
            values = [_to_float(tok) for tok in parts]
        except Exception as exc:  # noqa: BLE001
            issues.append(ParseIssue(line_no, f"float_parse_error={exc}", raw))
            continue

        rows.append(values)

    if not rows:
        raise RuntimeError(f"No data rows parsed from {path}")

    return np.array(rows, dtype=float), issues


def summarize_grid(data: np.ndarray, expected_l: Optional[int], expected_m: Optional[int]) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    l_vals = data[:, 0].astype(int)
    m_vals = data[:, 1].astype(int)
    seen = {(int(l), int(m)) for l, m in zip(l_vals, m_vals)}

    if expected_l is None:
        expected_l = int(l_vals.max())
    if expected_m is None:
        expected_m = int(m_vals.max())

    expected = {(l, m) for l in range(1, expected_l + 1) for m in range(1, expected_m + 1)}
    missing = sorted(expected - seen)

    dup_count: dict[tuple[int, int], int] = {}
    for pair in zip(l_vals, m_vals):
        p = (int(pair[0]), int(pair[1]))
        dup_count[p] = dup_count.get(p, 0) + 1
    duplicates = sorted([k for k, v in dup_count.items() if v > 1])

    return missing, duplicates


def write_csv(data: np.ndarray, out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(COLUMNS)
        for row in data:
            writer.writerow(row.tolist())


def plot_field(data: np.ndarray, field: str, outdir: Path) -> None:
    idx = COLUMNS.index(field)
    x = data[:, 2]
    y = data[:, 3]
    z = data[:, idx]

    outdir.mkdir(parents=True, exist_ok=True)

    # Scatter map in physical coordinates.
    plt.figure(figsize=(8, 5))
    sc = plt.scatter(x, y, c=z, cmap="viridis", s=26)
    plt.colorbar(sc, label=field)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"{field} scatter on (X,Y)")
    plt.tight_layout()
    plt.savefig(outdir / f"{field.lower()}_scatter.png", dpi=170)
    plt.close()

    # Triangulated contour for smoother view.
    plt.figure(figsize=(8, 5))
    try:
        ctf = plt.tricontourf(x, y, z, levels=24, cmap="viridis")
        plt.tricontour(x, y, z, levels=12, colors="k", linewidths=0.35, alpha=0.45)
        plt.colorbar(ctf, label=field)
    except Exception:  # noqa: BLE001
        plt.scatter(x, y, c=z, cmap="viridis", s=26)
        plt.colorbar(label=field)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"{field} contour on (X,Y)")
    plt.tight_layout()
    plt.savefig(outdir / f"{field.lower()}_contour.png", dpi=170)
    plt.close()


def plot_lm_heatmap(data: np.ndarray, field: str, expected_l: Optional[int], expected_m: Optional[int], outdir: Path) -> None:
    idx = COLUMNS.index(field)
    l = data[:, 0].astype(int)
    m = data[:, 1].astype(int)
    z = data[:, idx]

    l_max = expected_l or int(l.max())
    m_max = expected_m or int(m.max())

    grid = np.full((m_max, l_max), np.nan, dtype=float)
    for li, mi, zi in zip(l, m, z):
        if 1 <= li <= l_max and 1 <= mi <= m_max:
            grid[mi - 1, li - 1] = zi

    outdir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 4))
    im = plt.imshow(grid, aspect="auto", origin="lower", cmap="viridis")
    plt.colorbar(im, label=field)
    plt.xlabel("L index")
    plt.ylabel("M index")
    plt.title(f"{field} on (L,M) index grid")
    plt.tight_layout()
    plt.savefig(outdir / f"{field.lower()}_lm_heatmap.png", dpi=170)
    plt.close()


def detect_outliers(data: np.ndarray, field: str, z_thresh: float = 4.0) -> list[int]:
    idx = COLUMNS.index(field)
    vals = data[:, idx]
    mean = float(np.nanmean(vals))
    std = float(np.nanstd(vals))
    if std == 0.0 or math.isnan(std):
        return []
    z = np.abs((vals - mean) / std)
    return [int(i) for i in np.where(z > z_thresh)[0]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot NAP longtable solution data for OCR QA")
    parser.add_argument("--input", required=True, type=Path, help="Path to LaTeX longtable file")
    parser.add_argument("--field", default="Mach", choices=COLUMNS[2:], help="Single field to plot")
    parser.add_argument("--all-fields", action="store_true", help="Plot all fields from X through T")
    parser.add_argument("--expected-l", type=int, default=None, help="Expected max L index")
    parser.add_argument("--expected-m", type=int, default=None, help="Expected max M index")
    parser.add_argument("--outdir", type=Path, default=Path("artifacts/table_qc"), help="Output directory for plots/reports")
    args = parser.parse_args()

    data, issues = parse_longtable(args.input)
    missing, duplicates = summarize_grid(data, args.expected_l, args.expected_m)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(data, args.outdir / "parsed_table.csv")
    fields = COLUMNS[2:] if args.all_fields else [args.field]
    for field in fields:
        plot_field(data, field, args.outdir)
        plot_lm_heatmap(data, field, args.expected_l, args.expected_m, args.outdir)

    outliers = {field: detect_outliers(data, field) for field in fields}

    report = args.outdir / "qa_report.txt"
    with report.open("w", encoding="utf-8") as f:
        f.write(f"Input: {args.input}\n")
        f.write(f"Rows parsed: {len(data)}\n")
        f.write(f"Parse issues: {len(issues)}\n")
        f.write(f"Missing (L,M) pairs: {len(missing)}\n")
        f.write(f"Duplicate (L,M) pairs: {len(duplicates)}\n")
        for field in fields:
            f.write(f"Outliers in {field} (|z|>4): {len(outliers[field])}\n")
        if missing:
            f.write(f"Missing pairs: {missing}\n")
        if duplicates:
            f.write(f"Duplicate pairs: {duplicates}\n")
        for field in fields:
            if outliers[field]:
                f.write(
                    f"Outlier rows for {field} (0-based in parsed_table.csv): "
                    + ", ".join(map(str, outliers[field]))
                    + "\n"
                )
        if issues:
            f.write("\nParse issues (line_no | reason | raw):\n")
            for issue in issues[:200]:
                f.write(f"{issue.line_no} | {issue.reason} | {issue.raw}\n")

    print(f"Wrote: {args.outdir / 'parsed_table.csv'}")
    for field in fields:
        print(f"Wrote: {args.outdir / (field.lower() + '_scatter.png')}")
        print(f"Wrote: {args.outdir / (field.lower() + '_contour.png')}")
        print(f"Wrote: {args.outdir / (field.lower() + '_lm_heatmap.png')}")
    print(f"Wrote: {report}")


if __name__ == "__main__":
    main()
