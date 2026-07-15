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
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "workflow" / "step04_initial"


def normalize_angle(theta: np.ndarray | float) -> np.ndarray | float:
    return np.mod(theta, 2.0 * math.pi)


def angle_delta(target: float, source: float) -> float:
    return float(np.angle(np.exp(1j * (target - source))))


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


def square_from_center(record: dict[str, float | int]) -> np.ndarray:
    x = float(record["x"])
    y = float(record["y"])
    side = float(record["side"])
    half = side / 2.0
    return np.array([[x - half, y + half], [x + half, y + half], [x + half, y - half], [x - half, y - half]], dtype=float)


def cell_centers(records: list[dict[str, float | int]]) -> np.ndarray:
    return np.asarray([[float(record["x"]), float(record["y"])] for record in records], dtype=float)


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


def assign_layers(
    radius: np.ndarray, layer_count: int, *, inner_radius: float = 1.0, outer_radius: float = 1.85
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    quantiles = np.quantile(radius, np.linspace(0.0, 1.0, layer_count + 1))
    quantiles[0] -= 1e-9
    quantiles[-1] += 1e-9
    layer_ids = np.searchsorted(quantiles, radius, side="right") - 1
    layer_ids = np.clip(layer_ids, 0, layer_count - 1).astype(int)
    layer_radii = np.asarray([float(np.median(radius[layer_ids == layer])) for layer in range(layer_count)], dtype=float)
    layer_edges = np.empty(layer_count + 1, dtype=float)
    layer_edges[0] = inner_radius
    layer_edges[-1] = outer_radius
    for layer in range(1, layer_count):
        layer_edges[layer] = float(np.clip(0.5 * (layer_radii[layer - 1] + layer_radii[layer]), inner_radius, outer_radius))
    return layer_ids, layer_radii, layer_edges


def arc_points(radius: float, theta_a: float, theta_b: float, samples: int = 24) -> np.ndarray:
    delta = angle_delta(theta_b, theta_a)
    theta = theta_a + np.linspace(0.0, delta, samples)
    return np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])


def radial_curve_points(r_a: float, theta_a: float, r_b: float, theta_b: float, samples: int = 24) -> np.ndarray:
    delta = angle_delta(theta_b, theta_a)
    t = np.linspace(0.0, 1.0, samples)
    radius = (1.0 - t) * r_a + t * r_b
    theta = theta_a + delta * t
    return np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])


def point_to_point_polar_curve(
    point_a: np.ndarray,
    point_b: np.ndarray,
    *,
    samples: int = 24,
) -> np.ndarray:
    r_a = float(np.hypot(point_a[0], point_a[1]))
    r_b = float(np.hypot(point_b[0], point_b[1]))
    theta_a = float(np.arctan2(point_a[1], point_a[0]))
    theta_b = float(np.arctan2(point_b[1], point_b[0]))
    curve = radial_curve_points(r_a, theta_a, r_b, theta_b, samples=samples)
    curve[0] = point_a
    curve[-1] = point_b
    return curve


def classify_topology_edges(
    edges: list[tuple[int, int]],
    centers: np.ndarray,
    radius: np.ndarray,
    theta: np.ndarray,
    layer_ids: np.ndarray,
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    tangential: list[np.ndarray] = []
    radial: list[np.ndarray] = []
    ignored: list[np.ndarray] = []
    for i, j in edges:
        layer_delta = abs(int(layer_ids[i]) - int(layer_ids[j]))
        if layer_delta == 0:
            tangential.append(point_to_point_polar_curve(centers[i], centers[j], samples=28))
        elif layer_delta == 1:
            radial.append(point_to_point_polar_curve(centers[i], centers[j], samples=28))
        else:
            ignored.append(np.vstack([centers[i], centers[j]]))
    return tangential, radial, ignored


def annular_cell_polygon(r_inner: float, r_outer: float, theta_left: float, theta_right: float, samples: int = 10) -> np.ndarray:
    if theta_right <= theta_left:
        theta_right += 2.0 * math.pi
    outer_theta = np.linspace(theta_left, theta_right, samples)
    inner_theta = np.linspace(theta_right, theta_left, samples)
    outer = np.column_stack([r_outer * np.cos(outer_theta), r_outer * np.sin(outer_theta)])
    inner = np.column_stack([r_inner * np.cos(inner_theta), r_inner * np.sin(inner_theta)])
    return np.vstack([outer, inner])


def circular_order_with_largest_gap(theta_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(theta_values)
    sorted_theta = np.asarray(theta_values[order], dtype=float)
    if len(sorted_theta) <= 1:
        return order, sorted_theta
    gaps = np.diff(np.r_[sorted_theta, sorted_theta[0] + 2.0 * math.pi])
    start = int((np.argmax(gaps) + 1) % len(sorted_theta))
    rotated_order = np.r_[order[start:], order[:start]]
    if start == 0:
        return rotated_order, sorted_theta.copy()
    rotated_theta = np.r_[sorted_theta[start:], sorted_theta[:start] + 2.0 * math.pi]
    return rotated_order, rotated_theta


def spread_close_angles(unwrapped_theta: np.ndarray, *, min_gap_deg: float = 1.0) -> np.ndarray:
    if len(unwrapped_theta) <= 1:
        return unwrapped_theta.copy()
    min_gap = math.radians(min_gap_deg)
    adjusted = unwrapped_theta.copy()
    raw_gaps = np.diff(unwrapped_theta)
    start = 0
    while start < len(unwrapped_theta):
        end = start
        while end < len(raw_gaps) and raw_gaps[end] < min_gap:
            end += 1
        if end > start:
            count = end - start + 1
            cluster_mean = float(np.mean(unwrapped_theta[start : end + 1]))
            proposed = cluster_mean + (np.arange(count) - (count - 1) / 2.0) * min_gap
            if start > 0:
                left_limit = adjusted[start - 1] + min_gap
                if proposed[0] < left_limit:
                    proposed += left_limit - proposed[0]
            if end + 1 < len(unwrapped_theta):
                right_limit = unwrapped_theta[end + 1] - min_gap
                if proposed[-1] > right_limit:
                    proposed += right_limit - proposed[-1]
            adjusted[start : end + 1] = proposed
        start = end + 1
    return adjusted


def partition_polygons(
    radius: np.ndarray,
    theta: np.ndarray,
    layer_ids: np.ndarray,
    layer_edges: np.ndarray,
    cell_ids: np.ndarray,
) -> tuple[list[np.ndarray], list[int], list[dict[str, float | int]]]:
    polygons_by_index: list[np.ndarray | None] = [None] * len(radius)
    ids: list[int] = []
    rows: list[dict[str, float | int]] = []
    for layer in range(len(layer_edges) - 1):
        indices = np.flatnonzero(layer_ids == layer)
        if len(indices) == 0:
            continue
        local_order, unwrapped = circular_order_with_largest_gap(theta[indices])
        order = indices[local_order]
        layout_theta = spread_close_angles(unwrapped)
        if len(order) == 1:
            left = float(layout_theta[0] - math.pi / 18.0)
            right = float(layout_theta[0] + math.pi / 18.0)
            bounds = [(left, right)]
        else:
            bounds = []
            for pos, idx in enumerate(order):
                current = float(layout_theta[pos])
                prev_theta = float(layout_theta[pos - 1] if pos > 0 else layout_theta[-1] - 2.0 * math.pi)
                next_theta = float(layout_theta[pos + 1] if pos < len(order) - 1 else layout_theta[0] + 2.0 * math.pi)
                left = 0.5 * (prev_theta + current)
                right = 0.5 * (current + next_theta)
                bounds.append((left, right))
        for slot, (idx, (left, right)) in enumerate(zip(order, bounds), start=1):
            poly = annular_cell_polygon(float(layer_edges[layer]), float(layer_edges[layer + 1]), left, right)
            polygons_by_index[int(idx)] = poly
            rows.append(
                {
                    "cell_id": int(cell_ids[idx]),
                    "layer": int(layer + 1),
                    "slot_in_layer": int(slot),
                    "source_radius": float(radius[idx]),
                    "source_theta_deg": math.degrees(float(theta[idx])),
                    "theta_left_deg": math.degrees(left % (2.0 * math.pi)),
                    "theta_right_deg": math.degrees(right % (2.0 * math.pi)),
                    "r_inner": float(layer_edges[layer]),
                    "r_outer": float(layer_edges[layer + 1]),
                }
            )
    polygons = []
    for index, poly in enumerate(polygons_by_index):
        if poly is None:
            continue
        polygons.append(poly)
        ids.append(int(cell_ids[index]))
    return polygons, ids, rows


def add_direction_labels(ax, radius: float) -> None:
    for angle, label, ha, va in [
        (math.pi / 2, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3 * math.pi / 2, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ]:
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=11, color="#444444")


def setup_axis(ax, title: str, limit: float) -> None:
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.22)


def save_scaffold_figure(
    output_path: Path,
    records: list[dict[str, float | int]],
    centers: np.ndarray,
    radius: np.ndarray,
    theta: np.ndarray,
    layer_edges: np.ndarray,
    tangential: list[np.ndarray],
    radial: list[np.ndarray],
    ignored: list[np.ndarray],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=280)
    fig.patch.set_facecolor("white")
    circle = np.linspace(0.0, 2.0 * math.pi, 721)
    for r in layer_edges:
        ax.plot(r * np.cos(circle), r * np.sin(circle), color="#d0d0d0", linewidth=0.45, zorder=1)
    if ignored:
        ax.add_collection(LineCollection(ignored, colors="#cfcfcf", linewidths=0.35, alpha=0.45, zorder=2))
    if tangential:
        ax.add_collection(LineCollection(tangential, colors="#2870b8", linewidths=0.72, alpha=0.88, zorder=3))
    if radial:
        ax.add_collection(LineCollection(radial, colors="#c9483a", linewidths=0.72, alpha=0.88, zorder=4))
    ax.scatter(centers[:, 0], centers[:, 1], s=6.5, color="#111111", linewidths=0, zorder=5)
    add_direction_labels(ax, float(np.max(layer_edges)) + 0.12)
    patches = [
        mpl.patches.Patch(facecolor="#2870b8", edgecolor="none", label="same-layer arc"),
        mpl.patches.Patch(facecolor="#c9483a", edgecolor="none", label="adjacent-layer radial curve"),
        mpl.patches.Patch(facecolor="#cfcfcf", edgecolor="none", label="ignored long diagonal"),
    ]
    ax.legend(handles=patches, loc="upper right", frameon=True, fontsize=8)
    setup_axis(ax, "Line.md Step 1 - Center-polar Annulus-adapted Scaffold", float(np.max(layer_edges)) + 0.22)
    fig.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_partition_figure(
    output_path: Path,
    polygons: list[np.ndarray],
    ids: list[int],
    layer_edges: np.ndarray,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=280)
    fig.patch.set_facecolor("white")
    circle = np.linspace(0.0, 2.0 * math.pi, 721)
    for r in layer_edges:
        ax.plot(
            r * np.cos(circle),
            r * np.sin(circle),
            color="#222222" if r in (layer_edges[0], layer_edges[-1]) else "#d0d0d0",
            linewidth=0.7,
            zorder=1,
        )
    collection = PolyCollection(polygons, facecolors="none", edgecolors="#565656", linewidths=0.38, zorder=2)
    ax.add_collection(collection)
    for poly, cell_id in zip(polygons, ids):
        xy = np.mean(poly, axis=0)
        ax.text(float(xy[0]), float(xy[1]), str(cell_id), ha="center", va="center", fontsize=3.3, color="#111111", zorder=3)
    add_direction_labels(ax, float(np.max(layer_edges)) + 0.12)
    setup_axis(ax, "Step 4 - Center-polar Layer x Angular Midline Partition", float(np.max(layer_edges)) + 0.22)
    fig.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def write_assignment_csv(output_path: Path, rows: list[dict[str, float | int]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cell_id",
        "layer",
        "slot_in_layer",
        "source_radius",
        "source_theta_deg",
        "theta_left_deg",
        "theta_right_deg",
        "r_inner",
        "r_outer",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{value:.10f}" if isinstance(value, float) else value for key, value in row.items()})


def run(args: argparse.Namespace) -> None:
    records = read_center_polar_cells(args.input_csv)
    centers = cell_centers(records)
    radius = np.hypot(centers[:, 0], centers[:, 1])
    theta = normalize_angle(np.arctan2(centers[:, 1], centers[:, 0]))
    layer_ids, _, layer_edges = assign_layers(radius, args.layers)
    cell_ids = np.asarray([int(record["cell_id"]) for record in records], dtype=int)
    polygons, ids, rows = partition_polygons(radius, theta, layer_ids, layer_edges, cell_ids)
    save_partition_figure(args.output_dir / "step04_line_md_layer_column_partition.png", polygons, ids, layer_edges)
    write_assignment_csv(args.output_dir / "step04_line_md_layer_column_assignment.csv", rows)
    print(f"cells: {len(records)}")
    print(f"layers: {args.layers}")
    print(f"partition cells: {len(polygons)}")
    print(f"wrote: {args.output_dir / 'step04_line_md_layer_column_partition.png'}")
    print(f"wrote: {args.output_dir / 'step04_line_md_layer_column_assignment.csv'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step 4: convert center-polar layers and angular slots into annular partition cells.")
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--layers", type=int, default=6)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
