#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

TAU = 2.0 * math.pi


@dataclass(frozen=True)
class CellRecord:
    index: int
    cell_id: int
    row: int
    col: int
    x: float
    y: float
    radius: float
    theta: float
    side: float


def normalize_angle(theta: np.ndarray | float) -> np.ndarray | float:
    return np.mod(theta, TAU)


def signed_angle_delta(target: float, source: float) -> float:
    return float(np.angle(np.exp(1j * (target - source))))


def circular_angle_distance(a: float, b: float) -> float:
    return abs(signed_angle_delta(a, b))


def unwrap_angle(theta: float, seam_angle: float) -> float:
    """Return theta in the common [seam, seam + 2pi) coordinate system."""
    theta = float(normalize_angle(theta))
    seam = float(normalize_angle(seam_angle))
    if theta < seam:
        theta += TAU
    return theta


def unwrap_angles(theta: np.ndarray, seam_angle: float) -> np.ndarray:
    seam = float(normalize_angle(seam_angle))
    out = np.asarray(normalize_angle(theta), dtype=float).copy()
    out[out < seam] += TAU
    return out


def choose_global_seam_angle(theta: np.ndarray) -> float:
    values = np.sort(np.asarray(normalize_angle(theta), dtype=float))
    if len(values) == 0:
        raise ValueError("Cannot choose a seam angle from an empty theta array.")
    if len(values) == 1:
        return float(values[0])
    gaps = np.diff(np.r_[values, values[0] + TAU])
    start = int((np.argmax(gaps) + 1) % len(values))
    return float(values[start])


def read_center_polar_cells(path: Path) -> list[CellRecord]:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV does not exist: {path}")
    records: list[CellRecord] = []
    with path.open(encoding="utf-8") as f:
        for index, item in enumerate(csv.DictReader(f)):
            cell_id = int(item["cell_id"])
            row = int(item["row"])
            col = int(item["col"])
            x = float(item.get("X_center", item.get("x", 0.0)))
            y = float(item.get("Y_center", item.get("y", 0.0)))
            radius = float(item["R"]) if item.get("R") not in (None, "") else float(math.hypot(x, y))
            if item.get("theta_rad") not in (None, ""):
                theta = float(item["theta_rad"])
            elif item.get("theta_deg") not in (None, ""):
                theta = math.radians(float(item["theta_deg"]))
            else:
                theta = math.atan2(y, x)
            side = float(item.get("fixed_cell_side", item.get("side", 0.0)) or 0.0)
            records.append(
                CellRecord(
                    index=index,
                    cell_id=cell_id,
                    row=row,
                    col=col,
                    x=x,
                    y=y,
                    radius=radius,
                    theta=float(normalize_angle(theta)),
                    side=side,
                )
            )
    validate_records(records)
    return records


def validate_records(records: list[CellRecord]) -> None:
    if not records:
        raise ValueError("Input contains no cells.")
    cell_ids = [item.cell_id for item in records]
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError("cell_id values must be unique.")
    grid_keys = [(item.row, item.col) for item in records]
    if len(grid_keys) != len(set(grid_keys)):
        raise ValueError("row/col grid coordinates must be unique.")


def arrays_from_records(records: list[CellRecord]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cell_ids = np.asarray([item.cell_id for item in records], dtype=int)
    rows = np.asarray([item.row for item in records], dtype=int)
    cols = np.asarray([item.col for item in records], dtype=int)
    radius = np.asarray([item.radius for item in records], dtype=float)
    theta = np.asarray([item.theta for item in records], dtype=float)
    centers = np.asarray([[item.x, item.y] for item in records], dtype=float)
    return cell_ids, rows, cols, radius, theta, centers


def build_original_edges(records: list[CellRecord]) -> set[tuple[int, int]]:
    index_by_grid = {(item.row, item.col): item.index for item in records}
    edges: set[tuple[int, int]] = set()
    for item in records:
        for key in ((item.row, item.col + 1), (item.row + 1, item.col)):
            j = index_by_grid.get(key)
            if j is not None:
                a, b = sorted((item.index, int(j)))
                edges.add((a, b))
    return edges


def build_neighbors(node_count: int, edges: Iterable[tuple[int, int]]) -> list[list[int]]:
    neighbors = [[] for _ in range(node_count)]
    for i, j in sorted(edges):
        neighbors[i].append(j)
        neighbors[j].append(i)
    return neighbors


def read_topology_constraints(path: Path, records: list[CellRecord]) -> set[tuple[int, int]]:
    if not path.exists():
        return build_original_edges(records)
    index_by_cell_id = {item.cell_id: item.index for item in records}
    edges: set[tuple[int, int]] = set()
    with path.open(encoding="utf-8") as f:
        for item in csv.DictReader(f):
            if item.get("index_a") not in (None, "") and item.get("index_b") not in (None, ""):
                i = int(item["index_a"])
                j = int(item["index_b"])
            else:
                i = index_by_cell_id[int(item["cell_id_a"])]
                j = index_by_cell_id[int(item["cell_id_b"])]
            if i == j:
                continue
            a, b = sorted((i, j))
            edges.add((a, b))
    return edges


def assign_layers_by_radius(
    radius: np.ndarray,
    layer_count: int,
    *,
    inner_radius: float,
    outer_radius: float,
) -> tuple[np.ndarray, np.ndarray]:
    if layer_count <= 0:
        raise ValueError("--layers must be positive.")
    if layer_count > len(radius):
        raise ValueError(f"--layers ({layer_count}) cannot exceed the number of cells ({len(radius)}).")
    quantiles = np.quantile(radius, np.linspace(0.0, 1.0, layer_count + 1))
    quantiles[0] -= 1e-9
    quantiles[-1] += 1e-9
    layer_ids = np.searchsorted(quantiles, radius, side="right") - 1
    layer_ids = np.clip(layer_ids, 0, layer_count - 1).astype(int)
    if len(set(layer_ids.tolist())) != layer_count:
        raise ValueError("Initial quantile layer assignment produced an empty layer.")
    layer_edges = np.linspace(inner_radius, outer_radius, layer_count + 1, dtype=float)
    return layer_ids, layer_edges


def edge_layer_type(layer_delta: int) -> str:
    if layer_delta == 0:
        return "same_layer"
    if layer_delta == 1:
        return "adjacent_layer"
    return "long_layer_conflict"


def count_long_layer_conflicts(edges: set[tuple[int, int]], layer_ids: np.ndarray) -> int:
    return sum(1 for i, j in edges if abs(int(layer_ids[i]) - int(layer_ids[j])) > 1)


def annular_sector_polygon(r_inner: float, r_outer: float, theta_left: float, theta_right: float, samples: int = 12) -> np.ndarray:
    if theta_right <= theta_left:
        raise ValueError(f"Invalid annular sector with non-positive width: left={theta_left}, right={theta_right}")
    outer_theta = np.linspace(theta_left, theta_right, samples)
    inner_theta = np.linspace(theta_right, theta_left, samples)
    outer = np.column_stack([r_outer * np.cos(outer_theta), r_outer * np.sin(outer_theta)])
    inner = np.column_stack([r_inner * np.cos(inner_theta), r_inner * np.sin(inner_theta)])
    return np.vstack([outer, inner])


def split_circular_interval(left: float, right: float) -> list[tuple[float, float]]:
    left_mod = float(normalize_angle(left))
    width = right - left
    if width <= 0.0:
        return []
    if width >= TAU - 1e-12:
        return [(0.0, TAU)]
    right_mod = left_mod + width
    if right_mod <= TAU:
        return [(left_mod, right_mod)]
    return [(left_mod, TAU), (0.0, right_mod - TAU)]


def interval_overlap(left_a: float, right_a: float, left_b: float, right_b: float) -> float:
    overlap = 0.0
    for a0, a1 in split_circular_interval(left_a, right_a):
        for b0, b1 in split_circular_interval(left_b, right_b):
            overlap += max(0.0, min(a1, b1) - max(a0, b0))
    return float(overlap)


def linear_interval_overlap(left_a: float, right_a: float, left_b: float, right_b: float) -> float:
    return float(max(0.0, min(right_a, right_b) - max(left_a, left_b)))


def build_layer_orders(layer_ids: np.ndarray, unwrapped_theta: np.ndarray, layer_count: int) -> dict[int, list[int]]:
    orders: dict[int, list[int]] = {}
    for layer in range(layer_count):
        indices = np.flatnonzero(layer_ids == layer)
        if len(indices) == 0:
            raise ValueError(f"Layer {layer} is empty.")
        ordered = sorted(indices.tolist(), key=lambda idx: (float(unwrapped_theta[idx]), int(idx)))
        orders[layer] = ordered
    return orders


def initialize_boundaries_from_orders(
    orders: dict[int, list[int]],
    unwrapped_theta: np.ndarray,
    seam_angle: float,
    *,
    min_cell_angle: float,
) -> dict[int, np.ndarray]:
    seam = float(normalize_angle(seam_angle))
    end = seam + TAU
    boundaries: dict[int, np.ndarray] = {}
    for layer, order in sorted(orders.items()):
        n = len(order)
        if n * min_cell_angle > TAU + 1e-12:
            raise ValueError(f"Layer {layer} has too many cells for min_cell_angle={min_cell_angle}.")
        if n == 1:
            boundaries[layer] = np.asarray([seam, end], dtype=float)
            continue
        centers = np.asarray([float(unwrapped_theta[idx]) for idx in order], dtype=float)
        internal = 0.5 * (centers[:-1] + centers[1:])
        candidate = np.r_[seam, internal, end].astype(float)
        boundaries[layer] = enforce_min_width_boundaries(candidate, min_cell_angle, seam, end)
    return boundaries


def enforce_min_width_boundaries(boundaries: np.ndarray, min_width: float, start: float, end: float) -> np.ndarray:
    widths = np.diff(boundaries)
    if np.all(widths >= min_width - 1e-12):
        return boundaries.copy()
    n = len(boundaries) - 1
    available = (end - start) - n * min_width
    if available < -1e-12:
        raise ValueError("Not enough angular span to satisfy min_cell_angle.")
    centers = 0.5 * (boundaries[:-1] + boundaries[1:])
    order = np.argsort(centers)
    raw = np.diff(np.r_[start, centers[order], end])
    positive = np.maximum(raw[:-1] + raw[1:], 1e-9)
    extra = available * positive / float(np.sum(positive))
    widths_new = min_width + extra
    adjusted = np.r_[start, start + np.cumsum(widths_new)]
    adjusted[-1] = end
    return adjusted.astype(float)


def intervals_by_index(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
) -> dict[int, tuple[float, float]]:
    intervals: dict[int, tuple[float, float]] = {}
    for layer, order in sorted(orders.items()):
        b = boundaries[layer]
        if len(b) != len(order) + 1:
            raise ValueError(f"Layer {layer} has {len(order)} cells but {len(b)} boundaries.")
        for pos, idx in enumerate(order):
            if int(layer_ids[idx]) != layer:
                raise ValueError(f"Cell index {idx} is in order for layer {layer}, but layer_ids says {layer_ids[idx]}.")
            intervals[int(idx)] = (float(b[pos]), float(b[pos + 1]))
    return intervals


def build_partition_adjacency(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    min_shared_angle: float,
) -> set[tuple[int, int]]:
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    final_edges: set[tuple[int, int]] = set()
    for _layer, order in sorted(orders.items()):
        if len(order) <= 1:
            continue
        for pos, i in enumerate(order):
            j = order[(pos + 1) % len(order)]
            final_edges.add(tuple(sorted((int(i), int(j)))))
    layers = sorted(orders)
    for layer in layers:
        next_layer = layer + 1
        if next_layer not in orders:
            continue
        for i in orders[layer]:
            left_i, right_i = intervals[i]
            for j in orders[next_layer]:
                left_j, right_j = intervals[j]
                if linear_interval_overlap(left_i, right_i, left_j, right_j) >= min_shared_angle:
                    final_edges.add(tuple(sorted((int(i), int(j)))))
    return final_edges


def compute_adjacency_metrics(
    original_edges: set[tuple[int, int]],
    final_edges: set[tuple[int, int]],
) -> dict[str, float | int]:
    preserved = original_edges & final_edges
    lost = original_edges - final_edges
    new = final_edges - original_edges
    precision = len(preserved) / len(final_edges) if final_edges else 0.0
    recall = len(preserved) / len(original_edges) if original_edges else 0.0
    f1 = 0.0 if precision + recall == 0.0 else 2.0 * precision * recall / (precision + recall)
    return {
        "original_edge_count": len(original_edges),
        "final_edge_count": len(final_edges),
        "preserved_edge_count": len(preserved),
        "lost_edge_count": len(lost),
        "new_edge_count": len(new),
        "adjacency_precision": precision,
        "adjacency_recall": recall,
        "adjacency_f1": f1,
    }
