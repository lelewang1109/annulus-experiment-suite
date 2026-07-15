#!/usr/bin/env python3
"""Render the canonical center-polar input as fixed numbered cells.

The original experiment derived these records from an external PM2.5/GeoJSON
pipeline.  For a portable repository, ``data/input/beijing_grid_cells.csv`` is
the versioned interface between data preparation and the annulus experiments.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_INPUT_CSV = PROJECT_ROOT / "data" / "input" / "beijing_grid_cells.csv"
DEFAULT_OUTPUT_PNG = PROJECT_ROOT / "results" / "workflow" / "step01_center_polar.png"


def read_cells(path: Path) -> list[dict[str, float | int]]:
    records: list[dict[str, float | int]] = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            records.append(
                {
                    "cell_id": int(row["cell_id"]),
                    "x": float(row["X_center"]),
                    "y": float(row["Y_center"]),
                    "side": float(row["fixed_cell_side"]),
                }
            )
    if not records:
        raise ValueError(f"No cell records found in {path}")
    return records


def square(record: dict[str, float | int]) -> np.ndarray:
    x, y, half = float(record["x"]), float(record["y"]), float(record["side"]) / 2.0
    return np.asarray(
        [[x - half, y + half], [x + half, y + half], [x + half, y - half], [x - half, y - half]],
        dtype=float,
    )


def render(input_csv: Path, output_png: Path, *, inner_radius: float, outer_radius: float) -> None:
    records = read_cells(input_csv)
    polygons = [square(record) for record in records]
    centers = np.asarray([[float(record["x"]), float(record["y"])] for record in records], dtype=float)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.2, 10.2), dpi=280)
    theta = np.linspace(0.0, 2.0 * math.pi, 721)
    for radius in (inner_radius, outer_radius):
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.add_collection(PolyCollection(polygons, facecolors="none", edgecolors="#696969", linewidths=0.38, zorder=2))
    ax.scatter(centers[:, 0], centers[:, 1], s=4.5, color="#111111", linewidths=0, zorder=3)
    for record, (x, y) in zip(records, centers):
        ax.text(x, y, str(record["cell_id"]), ha="center", va="center", fontsize=3.65, color="#111111", zorder=4)
    for angle, label, ha, va in (
        (math.pi / 2.0, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3.0 * math.pi / 2.0, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ):
        radius = outer_radius + 0.13
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=10, color="#444444")
    ax.set_title("Beijing Grid Center Polar Mapping to Annulus", loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    limit = max(float(np.abs(np.vstack(polygons)).max()), outer_radius) + 0.35
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.22)
    fig.savefig(output_png, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"cells: {len(records)}")
    print(f"wrote: {output_png}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--output-png", type=Path, default=DEFAULT_OUTPUT_PNG)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=1.85)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.outer_radius <= args.inner_radius:
        raise ValueError("--outer-radius must be greater than --inner-radius")
    render(args.input_csv, args.output_png, inner_radius=args.inner_radius, outer_radius=args.outer_radius)
