#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
DEFAULT_INPUT_CSV = PROJECT_ROOT / "data" / "input" / "beijing_grid_cells.csv"
DEFAULT_OUTPUT_PNG = PROJECT_ROOT / "results" / "workflow" / "step02_center_topology.png"


def read_center_polar_cells(path: Path) -> list[dict[str, float | int]]:
    records: list[dict[str, float | int]] = []
    with path.open(encoding="utf-8") as f:
        for item in csv.DictReader(f):
            records.append(
                {
                    "cell_id": int(item["cell_id"]),
                    "row": int(item["row"]),
                    "col": int(item["col"]),
                    "x": float(item["X_center"]),
                    "y": float(item["Y_center"]),
                    "side": float(item["fixed_cell_side"]),
                }
            )
    return records


def square_from_center(x: float, y: float, side: float) -> np.ndarray:
    half = side / 2.0
    return np.array(
        [
            [x - half, y + half],
            [x + half, y + half],
            [x + half, y - half],
            [x - half, y - half],
        ],
        dtype=float,
    )


def cell_center(record: dict[str, float | int]) -> np.ndarray:
    return np.array([float(record["x"]), float(record["y"])], dtype=float)


def topology_edges(records: list[dict[str, float | int]]) -> list[tuple[int, int]]:
    index_by_grid = {(int(item["row"]), int(item["col"])): i for i, item in enumerate(records)}
    edges: list[tuple[int, int]] = []
    for i, item in enumerate(records):
        row = int(item["row"])
        col = int(item["col"])
        for key in ((row, col + 1), (row + 1, col)):
            j = index_by_grid.get(key)
            if j is not None:
                edges.append((i, j))
    return edges


def save_topology_figure(
    output_path: Path,
    records: list[dict[str, float | int]],
    *,
    inner_radius: float,
    outer_radius: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    polygons = [square_from_center(float(item["x"]), float(item["y"]), float(item["side"])) for item in records]
    centers = np.asarray([cell_center(item) for item in records], dtype=float)
    edges = topology_edges(records)
    segments = [(centers[i], centers[j]) for i, j in edges]

    fig, ax = plt.subplots(figsize=(10.2, 10.2), dpi=280)
    fig.patch.set_facecolor("white")
    theta = np.linspace(0.0, 2.0 * math.pi, 721)
    ax.plot(inner_radius * np.cos(theta), inner_radius * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.plot(outer_radius * np.cos(theta), outer_radius * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)

    cell_collection = PolyCollection(polygons, facecolors="none", edgecolors="#b8b8b8", linewidths=0.32, zorder=2)
    ax.add_collection(cell_collection)
    edge_collection = LineCollection(segments, colors="#b64235", linewidths=0.62, alpha=0.88, zorder=3)
    ax.add_collection(edge_collection)
    ax.scatter(centers[:, 0], centers[:, 1], s=7.5, color="#1f5d9a", edgecolors="white", linewidths=0.22, zorder=4)

    ax.set_title("Center Polar Annulus - Original Grid Topology Connected by Centers", loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    coords = np.vstack(polygons + [centers])
    limit = max(float(np.abs(coords).max()), outer_radius) + 0.18
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.22)
    fig.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    print(f"cells: {len(records)}")
    print(f"topology edges: {len(edges)}")
    print(f"wrote: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Connect center-polar annulus cell centers by original grid topology.")
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--output-png", type=Path, default=DEFAULT_OUTPUT_PNG)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=1.85)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    save_topology_figure(
        args.output_png,
        read_center_polar_cells(args.input_csv),
        inner_radius=args.inner_radius,
        outer_radius=args.outer_radius,
    )
