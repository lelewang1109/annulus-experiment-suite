#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
WORKFLOW_DIR = THIS_DIR
DEFAULT_INPUT_CSV = PROJECT_ROOT / "data" / "input" / "beijing_grid_cells.csv"
DEFAULT_TOPOLOGY_CSV = PROJECT_ROOT / "results" / "workflow" / "step03" / "step03_topology_constraints.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "workflow" / "step04_topology_aware_v2"

for candidate in (THIS_DIR, WORKFLOW_DIR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from annulus_topology_utils import (  # noqa: E402
    TAU,
    arrays_from_records,
    assign_layers_by_radius,
    build_layer_orders,
    build_neighbors,
    choose_global_seam_angle,
    circular_angle_distance,
    compute_adjacency_metrics,
    count_long_layer_conflicts,
    interval_overlap,
    intervals_by_index,
    linear_interval_overlap,
    normalize_angle,
    read_center_polar_cells,
    read_topology_constraints,
    unwrap_angles,
)
from step04_build_topology_aware_partition import (  # noqa: E402
    PartitionState,
    angular_errors,
    old_step4_intervals,
    radial_inversion_count,
    rank_maps,
    setup_axis,
)


@dataclass(frozen=True)
class V2Weights:
    w_lost: float
    w_new: float
    w_layer_jump: float
    w_capacity: float
    w_angle: float
    w_angle_cap: float
    w_rank: float
    w_width_balance: float
    w_aspect: float
    w_degree_inflation: float
    w_nonedge_overlap: float
    w_radial_deviation: float
    w_radial_inversion: float
    w_empty_layer: float


@dataclass
class V2Metrics:
    original_edge_count: int
    raw_final_edge_count: int
    effective_final_edge_count: int
    preserved_edge_count: int
    lost_edge_count: int
    new_edge_count: int
    adjacency_precision: float
    adjacency_recall: float
    adjacency_f1: float
    raw_adjacency_precision: float
    raw_adjacency_recall: float
    raw_adjacency_f1: float
    long_layer_conflict_count: int
    long_layer_conflict_rate: float
    min_angular_width_deg: float
    max_angular_width_deg: float
    p05_angular_width_deg: float
    p95_angular_width_deg: float
    lower_bound_cell_count: int
    exceed_max_target_count: int
    mean_aspect_ratio: float
    min_aspect_ratio: float
    max_aspect_ratio: float
    aspect_ratio_violation_count: int
    mean_degree_inflation: float
    max_degree_inflation: int
    cells_degree_inflation_ge_3: int
    cells_degree_inflation_ge_5: int
    mean_angular_error_deg: float
    median_angular_error_deg: float
    p90_angular_error_deg: float
    p95_angular_error_deg: float
    max_angular_error_deg: float
    angle_error_above_10_count: int
    angle_error_above_15_count: int
    capacity_error: float
    max_capacity_deviation: float
    nonedge_overlap_error: float
    objective: float


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def portable_parameter(value: Any) -> Any:
    if not isinstance(value, Path):
        return value
    try:
        return str(value.resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        return str(value)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def target_layer_counts(total_cells: int, layer_edges: np.ndarray) -> np.ndarray:
    weights = np.diff(layer_edges**2)
    raw = total_cells * weights / float(np.sum(weights))
    counts = np.floor(raw).astype(int)
    remainder = int(total_cells - np.sum(counts))
    order = sorted(range(len(counts)), key=lambda k: (-(raw[k] - counts[k]), k))
    for k in order[:remainder]:
        counts[k] += 1
    return counts.astype(int)


def layer_counts(layer_ids: np.ndarray, layer_count: int) -> np.ndarray:
    return np.asarray([int(np.sum(layer_ids == k)) for k in range(layer_count)], dtype=int)


def capacity_error(layer_ids: np.ndarray, targets: np.ndarray) -> float:
    counts = layer_counts(layer_ids, len(targets))
    return float(np.sum(((counts - targets) / np.maximum(1, targets)) ** 2))


def capacity_bounds(targets: np.ndarray, tolerance: float) -> tuple[np.ndarray, np.ndarray]:
    lower = np.floor((1.0 - tolerance) * targets).astype(int)
    upper = np.ceil((1.0 + tolerance) * targets).astype(int)
    lower = np.maximum(1, lower)
    upper = np.maximum(lower, upper)
    return lower, upper


def adaptive_circle_points(radius: float, start: float, end: float, max_step_deg: float) -> np.ndarray:
    width_deg = abs(math.degrees(end - start))
    samples = max(4, int(math.ceil(width_deg / max_step_deg)) + 1)
    theta = np.linspace(start, end, samples)
    return np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])


def adaptive_annular_sector_polygon(
    r_inner: float,
    r_outer: float,
    theta_left: float,
    theta_right: float,
    max_step_deg: float,
) -> np.ndarray:
    if theta_right <= theta_left:
        raise ValueError("Annular polygon requires theta_right > theta_left.")
    outer = adaptive_circle_points(r_outer, theta_left, theta_right, max_step_deg)
    inner = adaptive_circle_points(r_inner, theta_right, theta_left, max_step_deg)
    return np.vstack([outer, inner])


def polygon_area(poly: np.ndarray) -> float:
    x = poly[:, 0]
    y = poly[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def polygon_centroid(poly: np.ndarray) -> np.ndarray:
    area = polygon_area(poly)
    if abs(area) < 1e-14:
        return np.mean(poly, axis=0)
    x = poly[:, 0]
    y = poly[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    cx = float(np.sum((x + np.roll(x, -1)) * cross) / (6.0 * area))
    cy = float(np.sum((y + np.roll(y, -1)) * cross) / (6.0 * area))
    return np.asarray([cx, cy], dtype=float)


def classify_edges_by_layer(
    edges: set[tuple[int, int]], layer_ids: np.ndarray
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], list[tuple[int, int]]]:
    same: list[tuple[int, int]] = []
    adjacent: list[tuple[int, int]] = []
    long: list[tuple[int, int]] = []
    for i, j in sorted(edges):
        delta = abs(int(layer_ids[i]) - int(layer_ids[j]))
        if delta == 0:
            same.append((i, j))
        elif delta == 1:
            adjacent.append((i, j))
        else:
            long.append((i, j))
    return same, adjacent, long


def build_adjacencies_v2(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    effective_overlap_ratio: float,
    *,
    contact_tolerance: float = 1e-10,
) -> tuple[set[tuple[int, int]], set[tuple[int, int]], dict[tuple[int, int], str]]:
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    raw_edges: set[tuple[int, int]] = set()
    effective_edges: set[tuple[int, int]] = set()
    edge_kind: dict[tuple[int, int], str] = {}
    for layer, order in sorted(orders.items()):
        if len(order) <= 1:
            continue
        for pos, i in enumerate(order):
            j = order[(pos + 1) % len(order)]
            edge = tuple(sorted((int(i), int(j))))
            raw_edges.add(edge)
            effective_edges.add(edge)
            edge_kind[edge] = "same_layer"
    for layer in sorted(orders):
        if layer + 1 not in orders:
            continue
        for i in orders[layer]:
            li, ri = intervals[i]
            wi = ri - li
            for j in orders[layer + 1]:
                lj, rj = intervals[j]
                wj = rj - lj
                ov = linear_interval_overlap(li, ri, lj, rj)
                if ov > contact_tolerance:
                    edge = tuple(sorted((int(i), int(j))))
                    raw_edges.add(edge)
                    edge_kind[edge] = "cross_layer"
                    if ov / max(1e-12, min(wi, wj)) >= effective_overlap_ratio:
                        effective_edges.add(edge)
    return raw_edges, effective_edges, edge_kind


def initialize_balanced_boundaries(
    orders: dict[int, list[int]],
    unwrapped_theta: np.ndarray,
    seam_angle: float,
    *,
    min_cell_angle: float,
    blend_to_equal: float = 0.25,
) -> dict[int, np.ndarray]:
    seam = float(normalize_angle(seam_angle))
    end = seam + TAU
    boundaries: dict[int, np.ndarray] = {}
    for layer, order in sorted(orders.items()):
        n = len(order)
        if n * min_cell_angle > TAU + 1e-12:
            raise ValueError(f"Layer {layer} has {n} cells, too many for min angle.")
        equal = np.linspace(seam, end, n + 1)
        if n == 1:
            boundaries[layer] = equal
            continue
        centers = np.asarray([float(unwrapped_theta[idx]) for idx in order], dtype=float)
        midpoint = np.r_[seam, 0.5 * (centers[:-1] + centers[1:]), end]
        candidate = (1.0 - blend_to_equal) * midpoint + blend_to_equal * equal
        candidate[0] = seam
        candidate[-1] = end
        boundaries[layer] = project_boundaries(candidate, min_cell_angle, seam, end)
    return boundaries


def project_boundaries(boundaries: np.ndarray, min_width: float, start: float, end: float) -> np.ndarray:
    b = boundaries.copy()
    b[0] = start
    b[-1] = end
    n = len(b) - 1
    if n * min_width > end - start + 1e-12:
        raise ValueError("min_cell_angle is too large for this layer.")
    for i in range(1, len(b)):
        b[i] = max(b[i], b[i - 1] + min_width)
    overflow = b[-1] - end
    if overflow > 0:
        b[1:-1] -= overflow * np.linspace(1.0 / n, (n - 1.0) / n, n - 1)
        b[-1] = end
        for i in range(len(b) - 2, -1, -1):
            b[i] = min(b[i], b[i + 1] - min_width)
        b[0] = start
    return b


def width_stats(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    *,
    min_cell_angle: float,
    max_angle_factor: float,
) -> dict[str, Any]:
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    widths = np.asarray([right - left for left, right in intervals.values()], dtype=float)
    max_targets = []
    for idx in range(len(layer_ids)):
        layer = int(layer_ids[idx])
        max_targets.append(max_angle_factor * TAU / len(orders[layer]))
    max_targets_arr = np.asarray(max_targets, dtype=float)
    return {
        "widths": widths,
        "max_targets": max_targets_arr,
        "min_width": float(np.min(widths)),
        "max_width": float(np.max(widths)),
        "p05_width": float(np.percentile(widths, 5)),
        "p95_width": float(np.percentile(widths, 95)),
        "lower_bound_count": int(np.sum(widths <= min_cell_angle + math.radians(0.02))),
        "exceed_max_target_count": int(np.sum(widths > max_targets_arr + 1e-12)),
        "width_error": float(
            np.mean(np.maximum(0.0, min_cell_angle - widths) ** 2 + np.maximum(0.0, widths - max_targets_arr) ** 2) / (math.pi**2)
        ),
    }


def aspect_stats(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    layer_edges: np.ndarray,
    *,
    min_aspect: float,
    max_aspect: float,
) -> dict[str, Any]:
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    aspects = []
    penalties = []
    for idx, (left, right) in intervals.items():
        layer = int(layer_ids[idx])
        mean_radius = 0.5 * (layer_edges[layer] + layer_edges[layer + 1])
        thickness = layer_edges[layer + 1] - layer_edges[layer]
        aspect = mean_radius * (right - left) / max(thickness, 1e-12)
        aspects.append(aspect)
        penalties.append(max(0.0, min_aspect - aspect) ** 2 + max(0.0, aspect - max_aspect) ** 2)
    arr = np.asarray(aspects, dtype=float)
    pen = np.asarray(penalties, dtype=float)
    return {
        "aspects": arr,
        "mean_aspect": float(np.mean(arr)),
        "min_aspect": float(np.min(arr)),
        "max_aspect": float(np.max(arr)),
        "violation_count": int(np.sum((arr < min_aspect) | (arr > max_aspect))),
        "aspect_error": float(np.mean(pen) / max(1.0, max_aspect**2)),
    }


def degree_inflation_stats(node_count: int, original_edges: set[tuple[int, int]], final_edges: set[tuple[int, int]]) -> dict[str, Any]:
    original_degree = np.asarray([len(x) for x in build_neighbors(node_count, original_edges)], dtype=int)
    final_degree = np.asarray([len(x) for x in build_neighbors(node_count, final_edges)], dtype=int)
    inflation = np.maximum(0, final_degree - original_degree)
    return {
        "original_degree": original_degree,
        "final_degree": final_degree,
        "inflation": inflation,
        "mean": float(np.mean(inflation)),
        "max": int(np.max(inflation)),
        "ge3": int(np.sum(inflation >= 3)),
        "ge5": int(np.sum(inflation >= 5)),
        "error": float(np.mean(inflation.astype(float) ** 2) / 25.0),
    }


def nonedge_overlap_error(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    original_edges: set[tuple[int, int]],
) -> float:
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    values = []
    for layer in sorted(orders):
        if layer + 1 not in orders:
            continue
        for i in orders[layer]:
            li, ri = intervals[i]
            wi = ri - li
            for j in orders[layer + 1]:
                edge = tuple(sorted((int(i), int(j))))
                if edge in original_edges:
                    continue
                lj, rj = intervals[j]
                wj = rj - lj
                ov = linear_interval_overlap(li, ri, lj, rj)
                values.append((ov / max(1e-12, min(wi, wj))) ** 2)
    return float(np.mean(values)) if values else 0.0


def comprehensive_objective(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    layer_edges: np.ndarray,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    targets: np.ndarray,
    weights: V2Weights,
    *,
    min_cell_angle: float,
    max_angle_factor: float,
    min_aspect: float,
    max_aspect: float,
    effective_overlap_ratio: float,
    max_angle_displacement: float,
) -> tuple[float, dict[str, Any]]:
    raw_edges, effective_edges, _edge_kind = build_adjacencies_v2(layer_ids, orders, boundaries, effective_overlap_ratio)
    eff_metrics = compute_adjacency_metrics(original_edges, effective_edges)
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    angle_err = angular_errors(theta, intervals)
    current_ranks = rank_maps(orders)
    rank_disp = [
        abs(current_ranks[idx] - original_ranks.get(idx, current_ranks[idx])) / max(1, len(orders[int(layer_ids[idx])]) - 1)
        for idx in range(len(layer_ids))
    ]
    raw_metrics = compute_adjacency_metrics(original_edges, raw_edges)
    width = width_stats(layer_ids, orders, boundaries, min_cell_angle=min_cell_angle, max_angle_factor=max_angle_factor)
    aspect = aspect_stats(layer_ids, orders, boundaries, layer_edges, min_aspect=min_aspect, max_aspect=max_aspect)
    degree = degree_inflation_stats(len(layer_ids), original_edges, effective_edges)
    long_rate = count_long_layer_conflicts(original_edges, layer_ids) / max(1, len(original_edges))
    cap_error = capacity_error(layer_ids, targets)
    angle_cap_error = float(np.mean(np.maximum(0.0, angle_err - max_angle_displacement) ** 2) / (math.pi**2))
    nonedge_error = nonedge_overlap_error(layer_ids, orders, boundaries, original_edges)
    objective = (
        weights.w_lost * (eff_metrics["lost_edge_count"] / max(1, len(original_edges)))
        + weights.w_new * (eff_metrics["new_edge_count"] / max(1, len(effective_edges)))
        + weights.w_layer_jump * long_rate
        + weights.w_capacity * cap_error
        + weights.w_angle * float(np.mean(angle_err / math.pi))
        + weights.w_angle_cap * angle_cap_error
        + weights.w_rank * float(np.mean(rank_disp))
        + weights.w_width_balance * width["width_error"]
        + weights.w_aspect * aspect["aspect_error"]
        + weights.w_degree_inflation * degree["error"]
        + weights.w_nonedge_overlap * nonedge_error
    )
    details = {
        "raw_edges": raw_edges,
        "effective_edges": effective_edges,
        "raw_metrics": raw_metrics,
        "effective_metrics": eff_metrics,
        "angle_errors": angle_err,
        "rank_displacement": float(np.mean(rank_disp)),
        "width": width,
        "aspect": aspect,
        "degree": degree,
        "long_layer_conflict_count": count_long_layer_conflicts(original_edges, layer_ids),
        "long_layer_conflict_rate": long_rate,
        "capacity_error": cap_error,
        "capacity_deviation": ((layer_counts(layer_ids, len(targets)) - targets) / np.maximum(1, targets)).astype(float),
        "angle_cap_error": angle_cap_error,
        "nonedge_overlap_error": nonedge_error,
        "objective": float(objective),
    }
    return float(objective), details


def layer_objective_v2(
    layer_ids: np.ndarray,
    initial_layer_ids: np.ndarray,
    radius: np.ndarray,
    original_edges: set[tuple[int, int]],
    targets: np.ndarray,
    weights: V2Weights,
) -> float:
    jump = sum(max(0, abs(int(layer_ids[i]) - int(layer_ids[j])) - 1) ** 2 for i, j in original_edges) / max(1, len(original_edges))
    deviation = float(np.mean((layer_ids - initial_layer_ids) ** 2))
    inversion = radial_inversion_count(radius, layer_ids) / max(1, len(layer_ids) * (len(layer_ids) - 1) / 2)
    empty = len(targets) - len(set(layer_ids.tolist()))
    return (
        weights.w_layer_jump * jump
        + weights.w_radial_deviation * deviation
        + weights.w_radial_inversion * inversion
        + weights.w_capacity * capacity_error(layer_ids, targets)
        + weights.w_empty_layer * max(0, empty)
    )


def optimize_layer_assignment_v2(
    initial_layer_ids: np.ndarray,
    radius: np.ndarray,
    original_edges: set[tuple[int, int]],
    targets: np.ndarray,
    *,
    max_layer_shift: int,
    max_iters: int,
    capacity_tolerance: float,
    weights: V2Weights,
) -> tuple[np.ndarray, dict[str, Any]]:
    layer_count = len(targets)
    lower, upper = capacity_bounds(targets, capacity_tolerance)
    current = initial_layer_ids.copy()
    allowed = [set(range(max(0, int(k) - max_layer_shift), min(layer_count - 1, int(k) + max_layer_shift) + 1)) for k in initial_layer_ids]
    current_obj = layer_objective_v2(current, initial_layer_ids, radius, original_edges, targets, weights)
    best = current.copy()
    best_obj = current_obj
    history = [{"iteration": 0, "objective": current_obj, "counts": layer_counts(current, layer_count).tolist()}]

    for iteration in range(1, max_iters + 1):
        changed = False
        for idx in range(len(current)):
            old = int(current[idx])
            candidates = sorted(({old, old - 1, old + 1} & allowed[idx]))
            local_best_layer = old
            local_best_obj = current_obj
            for candidate in candidates:
                if candidate == old:
                    obj = current_obj
                else:
                    trial = current.copy()
                    trial[idx] = candidate
                    if len(set(trial.tolist())) != layer_count:
                        continue
                    obj = layer_objective_v2(trial, initial_layer_ids, radius, original_edges, targets, weights)
                if obj < local_best_obj - 1e-12:
                    local_best_obj = obj
                    local_best_layer = candidate
                elif abs(obj - local_best_obj) <= 1e-12:
                    if (abs(candidate - int(initial_layer_ids[idx])), candidate) < (
                        abs(local_best_layer - int(initial_layer_ids[idx])),
                        local_best_layer,
                    ):
                        local_best_layer = candidate
            if local_best_layer != old:
                current[idx] = local_best_layer
                current_obj = local_best_obj
                changed = True
                if current_obj < best_obj - 1e-12:
                    best = current.copy()
                    best_obj = current_obj

        # Deterministic repair for layers outside the soft capacity interval.
        repaired = True
        while repaired:
            repaired = False
            counts = layer_counts(current, layer_count)
            over_layers = [k for k in range(layer_count) if counts[k] > upper[k]]
            under_layers = [k for k in range(layer_count) if counts[k] < lower[k]]
            if not over_layers or not under_layers:
                break
            best_move: tuple[float, int, int] | None = None
            for src in over_layers:
                for idx in np.flatnonzero(current == src):
                    for dst in under_layers:
                        if dst not in allowed[int(idx)] or abs(dst - src) > 1:
                            continue
                        trial = current.copy()
                        trial[int(idx)] = dst
                        obj = layer_objective_v2(trial, initial_layer_ids, radius, original_edges, targets, weights)
                        key = (obj, int(idx), dst)
                        if best_move is None or key < best_move:
                            best_move = key
            if best_move is not None:
                _obj, idx, dst = best_move
                current[idx] = dst
                current_obj = layer_objective_v2(current, initial_layer_ids, radius, original_edges, targets, weights)
                repaired = True
                changed = True
                if current_obj < best_obj - 1e-12:
                    best = current.copy()
                    best_obj = current_obj
        history.append({"iteration": iteration, "objective": current_obj, "counts": layer_counts(current, layer_count).tolist()})
        if not changed:
            break

    # If historical best violates tolerance but current repair satisfies it, prefer the feasible state.
    counts_best = layer_counts(best, layer_count)
    counts_current = layer_counts(current, layer_count)
    if np.any((counts_best < lower) | (counts_best > upper)) and np.all((counts_current >= lower) & (counts_current <= upper)):
        best = current.copy()
        best_obj = current_obj
    return best, {
        "objective_before": history[0]["objective"],
        "objective_after": best_obj,
        "history": history,
        "target_layer_counts": targets.tolist(),
        "capacity_lower_bounds": lower.tolist(),
        "capacity_upper_bounds": upper.tolist(),
        "initial_layer_counts": layer_counts(initial_layer_ids, layer_count).tolist(),
        "optimized_layer_counts": layer_counts(best, layer_count).tolist(),
        "layer_capacity_deviation": ((layer_counts(best, layer_count) - targets) / np.maximum(1, targets)).tolist(),
        "layer_changed_cell_count": int(np.sum(best != initial_layer_ids)),
    }


def make_state_v2(
    name: str,
    layer_ids: np.ndarray,
    layer_edges: np.ndarray,
    unwrapped_theta: np.ndarray,
    seam_angle: float,
    *,
    min_cell_angle: float,
) -> PartitionState:
    orders = build_layer_orders(layer_ids, unwrapped_theta, len(layer_edges) - 1)
    boundaries = initialize_balanced_boundaries(orders, unwrapped_theta, seam_angle, min_cell_angle=min_cell_angle)
    return PartitionState(name, layer_ids.copy(), orders, boundaries, layer_edges.copy(), math.inf)


def evaluate_state(
    state: PartitionState,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    targets: np.ndarray,
    weights: V2Weights,
    args: argparse.Namespace,
) -> tuple[float, dict[str, Any]]:
    obj, details = comprehensive_objective(
        state.layer_ids,
        state.orders,
        state.boundaries,
        state.layer_edges,
        theta,
        original_edges,
        original_ranks,
        targets,
        weights,
        min_cell_angle=math.radians(args.min_cell_angle_deg),
        max_angle_factor=args.max_angle_factor,
        min_aspect=args.min_aspect,
        max_aspect=args.max_aspect,
        effective_overlap_ratio=args.target_overlap_ratio,
        max_angle_displacement=math.radians(args.max_angle_displacement_deg),
    )
    state.objective = obj
    return obj, details


def optimize_orders_v2(
    state: PartitionState,
    unwrapped_theta: np.ndarray,
    seam_angle: float,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    targets: np.ndarray,
    weights: V2Weights,
    args: argparse.Namespace,
) -> tuple[PartitionState, dict[str, Any]]:
    current = state.copy(name="order_initial")
    evaluate_state(current, theta, original_edges, original_ranks, targets, weights, args)
    best = current.copy(name="order_best")
    history = [{"iteration": 0, "objective": current.objective}]
    min_cell_angle = math.radians(args.min_cell_angle_deg)
    base_ranks = rank_maps(state.orders)

    def eval_orders(orders: dict[int, list[int]]) -> tuple[float, dict[int, np.ndarray]]:
        trial_ranks = rank_maps(orders)
        if any(abs(trial_ranks[idx] - base_ranks[idx]) > args.max_rank_shift for idx in trial_ranks):
            return math.inf, current.boundaries
        boundaries = initialize_balanced_boundaries(orders, unwrapped_theta, seam_angle, min_cell_angle=min_cell_angle)
        trial = PartitionState("trial", current.layer_ids.copy(), orders, boundaries, current.layer_edges.copy(), math.inf)
        obj, _details = evaluate_state(trial, theta, original_edges, original_ranks, targets, weights, args)
        return obj, boundaries

    for iteration in range(1, args.order_opt_iters + 1):
        improved = False
        for layer in sorted(current.orders):
            n = len(current.orders[layer])
            if n <= 1:
                continue
            for pos in range(n - 1):
                trial_orders = {k: v.copy() for k, v in current.orders.items()}
                trial_orders[layer][pos], trial_orders[layer][pos + 1] = trial_orders[layer][pos + 1], trial_orders[layer][pos]
                obj, boundaries = eval_orders(trial_orders)
                if obj < current.objective - 1e-12:
                    current.orders = trial_orders
                    current.boundaries = boundaries
                    current.objective = obj
                    improved = True
                    if obj < best.objective - 1e-12:
                        best = current.copy(name="order_best")
            for pos in range(n):
                accepted = False
                for delta in range(-args.max_rank_shift, args.max_rank_shift + 1):
                    if delta == 0:
                        continue
                    new_pos = pos + delta
                    if new_pos < 0 or new_pos >= n:
                        continue
                    trial_orders = {k: v.copy() for k, v in current.orders.items()}
                    item = trial_orders[layer].pop(pos)
                    trial_orders[layer].insert(new_pos, item)
                    obj, boundaries = eval_orders(trial_orders)
                    if obj < current.objective - 1e-12:
                        current.orders = trial_orders
                        current.boundaries = boundaries
                        current.objective = obj
                        improved = True
                        accepted = True
                        if obj < best.objective - 1e-12:
                            best = current.copy(name="order_best")
                        break
                if accepted:
                    break
        history.append({"iteration": iteration, "objective": current.objective})
        if not improved:
            break
    return best, {"objective_before": state.objective, "objective_after": best.objective, "history": history}


def optimize_boundaries_v2(
    state: PartitionState,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    targets: np.ndarray,
    weights: V2Weights,
    args: argparse.Namespace,
) -> tuple[PartitionState, dict[str, Any]]:
    current = state.copy(name="boundary_initial")
    evaluate_state(current, theta, original_edges, original_ranks, targets, weights, args)
    best = current.copy(name="boundary_best")
    history = [{"iteration": 0, "objective": current.objective}]
    min_cell_angle = math.radians(args.min_cell_angle_deg)
    for iteration in range(1, args.boundary_opt_iters + 1):
        improved = False
        for layer in sorted(current.boundaries):
            b = current.boundaries[layer].copy()
            if len(b) <= 2:
                continue
            for m in range(1, len(b) - 1):
                low = b[m - 1] + min_cell_angle
                high = b[m + 1] - min_cell_angle
                if high < low:
                    continue
                candidates = np.linspace(low, high, args.boundary_search_samples)
                best_value = b[m]
                best_obj = current.objective
                for value in candidates:
                    trial_boundaries = {k: v.copy() for k, v in current.boundaries.items()}
                    trial_boundaries[layer][m] = float(value)
                    trial = PartitionState(
                        "trial",
                        current.layer_ids.copy(),
                        {k: v.copy() for k, v in current.orders.items()},
                        trial_boundaries,
                        current.layer_edges.copy(),
                        math.inf,
                    )
                    obj, _details = evaluate_state(trial, theta, original_edges, original_ranks, targets, weights, args)
                    if obj < best_obj - 1e-12:
                        best_obj = obj
                        best_value = float(value)
                if abs(best_value - b[m]) > 1e-15:
                    b[m] = best_value
                    current.boundaries[layer] = b.copy()
                    current.objective = best_obj
                    improved = True
                    if best_obj < best.objective - 1e-12:
                        best = current.copy(name="boundary_best")
        history.append({"iteration": iteration, "objective": current.objective})
        if not improved:
            break
    return best, {"objective_before": state.objective, "objective_after": best.objective, "history": history}


def metrics_from_details(details: dict[str, Any]) -> V2Metrics:
    eff = details["effective_metrics"]
    raw = details["raw_metrics"]
    width = details["width"]
    aspect = details["aspect"]
    degree = details["degree"]
    angle = details["angle_errors"]
    return V2Metrics(
        original_edge_count=int(eff["original_edge_count"]),
        raw_final_edge_count=int(raw["final_edge_count"]),
        effective_final_edge_count=int(eff["final_edge_count"]),
        preserved_edge_count=int(eff["preserved_edge_count"]),
        lost_edge_count=int(eff["lost_edge_count"]),
        new_edge_count=int(eff["new_edge_count"]),
        adjacency_precision=float(eff["adjacency_precision"]),
        adjacency_recall=float(eff["adjacency_recall"]),
        adjacency_f1=float(eff["adjacency_f1"]),
        raw_adjacency_precision=float(raw["adjacency_precision"]),
        raw_adjacency_recall=float(raw["adjacency_recall"]),
        raw_adjacency_f1=float(raw["adjacency_f1"]),
        long_layer_conflict_count=int(details["long_layer_conflict_count"]),
        long_layer_conflict_rate=float(details["long_layer_conflict_rate"]),
        min_angular_width_deg=math.degrees(float(width["min_width"])),
        max_angular_width_deg=math.degrees(float(width["max_width"])),
        p05_angular_width_deg=math.degrees(float(width["p05_width"])),
        p95_angular_width_deg=math.degrees(float(width["p95_width"])),
        lower_bound_cell_count=int(width["lower_bound_count"]),
        exceed_max_target_count=int(width["exceed_max_target_count"]),
        mean_aspect_ratio=float(aspect["mean_aspect"]),
        min_aspect_ratio=float(aspect["min_aspect"]),
        max_aspect_ratio=float(aspect["max_aspect"]),
        aspect_ratio_violation_count=int(aspect["violation_count"]),
        mean_degree_inflation=float(degree["mean"]),
        max_degree_inflation=int(degree["max"]),
        cells_degree_inflation_ge_3=int(degree["ge3"]),
        cells_degree_inflation_ge_5=int(degree["ge5"]),
        mean_angular_error_deg=math.degrees(float(np.mean(angle))),
        median_angular_error_deg=math.degrees(float(np.median(angle))),
        p90_angular_error_deg=math.degrees(float(np.percentile(angle, 90))),
        p95_angular_error_deg=math.degrees(float(np.percentile(angle, 95))),
        max_angular_error_deg=math.degrees(float(np.max(angle))),
        angle_error_above_10_count=int(np.sum(angle > math.radians(10.0))),
        angle_error_above_15_count=int(np.sum(angle > math.radians(15.0))),
        capacity_error=float(details["capacity_error"]),
        max_capacity_deviation=float(np.max(np.abs(details["capacity_deviation"]))),
        nonedge_overlap_error=float(details["nonedge_overlap_error"]),
        objective=float(details["objective"]),
    )


def baseline_v2(
    theta: np.ndarray,
    initial_layer_ids: np.ndarray,
    layer_edges: np.ndarray,
    original_edges: set[tuple[int, int]],
    initial_ranks: dict[int, int],
    targets: np.ndarray,
    weights: V2Weights,
    args: argparse.Namespace,
) -> tuple[V2Metrics, set[tuple[int, int]], set[tuple[int, int]]]:
    orders_old, intervals_old = old_step4_intervals(theta, initial_layer_ids, len(layer_edges) - 1)
    boundaries: dict[int, np.ndarray] = {}
    for layer, order in orders_old.items():
        starts = [intervals_old[idx][0] for idx in order]
        ends = [intervals_old[idx][1] for idx in order]
        boundaries[layer] = np.asarray([starts[0], *ends], dtype=float)
        # Old Step 4 has per-layer seams; normalize only for metrics by preserving intervals.
    raw_edges = build_adjacency_from_explicit_intervals(
        initial_layer_ids, orders_old, intervals_old, 1e-10, raw=True, ratio=args.target_overlap_ratio
    )
    eff_edges = build_adjacency_from_explicit_intervals(
        initial_layer_ids, orders_old, intervals_old, 1e-10, raw=False, ratio=args.target_overlap_ratio
    )
    eff = compute_adjacency_metrics(original_edges, eff_edges)
    raw = compute_adjacency_metrics(original_edges, raw_edges)
    angle = angular_errors(theta, intervals_old)
    widths = np.asarray([right - left for left, right in intervals_old.values()], dtype=float)
    long_count = count_long_layer_conflicts(original_edges, initial_layer_ids)
    degree = degree_inflation_stats(len(theta), original_edges, eff_edges)
    metrics = V2Metrics(
        original_edge_count=len(original_edges),
        raw_final_edge_count=len(raw_edges),
        effective_final_edge_count=len(eff_edges),
        preserved_edge_count=int(eff["preserved_edge_count"]),
        lost_edge_count=int(eff["lost_edge_count"]),
        new_edge_count=int(eff["new_edge_count"]),
        adjacency_precision=float(eff["adjacency_precision"]),
        adjacency_recall=float(eff["adjacency_recall"]),
        adjacency_f1=float(eff["adjacency_f1"]),
        raw_adjacency_precision=float(raw["adjacency_precision"]),
        raw_adjacency_recall=float(raw["adjacency_recall"]),
        raw_adjacency_f1=float(raw["adjacency_f1"]),
        long_layer_conflict_count=long_count,
        long_layer_conflict_rate=long_count / max(1, len(original_edges)),
        min_angular_width_deg=math.degrees(float(np.min(widths))),
        max_angular_width_deg=math.degrees(float(np.max(widths))),
        p05_angular_width_deg=math.degrees(float(np.percentile(widths, 5))),
        p95_angular_width_deg=math.degrees(float(np.percentile(widths, 95))),
        lower_bound_cell_count=int(np.sum(widths <= math.radians(args.min_cell_angle_deg) + math.radians(0.02))),
        exceed_max_target_count=-1,
        mean_aspect_ratio=0.0,
        min_aspect_ratio=0.0,
        max_aspect_ratio=0.0,
        aspect_ratio_violation_count=-1,
        mean_degree_inflation=float(degree["mean"]),
        max_degree_inflation=int(degree["max"]),
        cells_degree_inflation_ge_3=int(degree["ge3"]),
        cells_degree_inflation_ge_5=int(degree["ge5"]),
        mean_angular_error_deg=math.degrees(float(np.mean(angle))),
        median_angular_error_deg=math.degrees(float(np.median(angle))),
        p90_angular_error_deg=math.degrees(float(np.percentile(angle, 90))),
        p95_angular_error_deg=math.degrees(float(np.percentile(angle, 95))),
        max_angular_error_deg=math.degrees(float(np.max(angle))),
        angle_error_above_10_count=int(np.sum(angle > math.radians(10.0))),
        angle_error_above_15_count=int(np.sum(angle > math.radians(15.0))),
        capacity_error=capacity_error(initial_layer_ids, targets),
        max_capacity_deviation=float(np.max(np.abs((layer_counts(initial_layer_ids, len(targets)) - targets) / np.maximum(1, targets)))),
        nonedge_overlap_error=0.0,
        objective=math.nan,
    )
    return metrics, raw_edges, eff_edges


def build_adjacency_from_explicit_intervals(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    intervals: dict[int, tuple[float, float]],
    contact_tolerance: float,
    *,
    raw: bool,
    ratio: float,
) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for _layer, order in orders.items():
        if len(order) <= 1:
            continue
        for pos, i in enumerate(order):
            edges.add(tuple(sorted((int(i), int(order[(pos + 1) % len(order)])))))
    for layer in sorted(orders):
        if layer + 1 not in orders:
            continue
        for i in orders[layer]:
            li, ri = intervals[i]
            wi = ri - li
            for j in orders[layer + 1]:
                lj, rj = intervals[j]
                wj = rj - lj
                ov = interval_overlap(li, ri, lj, rj)
                if raw:
                    if ov > contact_tolerance:
                        edges.add(tuple(sorted((int(i), int(j)))))
                elif ov / max(1e-12, min(wi, wj)) >= ratio:
                    edges.add(tuple(sorted((int(i), int(j)))))
    return edges


def polygons_from_state_v2(
    state: PartitionState, cell_ids: np.ndarray, max_arc_step_deg: float
) -> tuple[list[np.ndarray], list[int], dict[int, tuple[float, float, float, float]], np.ndarray]:
    intervals = intervals_by_index(state.layer_ids, state.orders, state.boundaries)
    polygons: list[np.ndarray] = []
    ids: list[int] = []
    geometry: dict[int, tuple[float, float, float, float]] = {}
    centroids = np.zeros((len(cell_ids), 2), dtype=float)
    for idx in range(len(cell_ids)):
        left, right = intervals[idx]
        layer = int(state.layer_ids[idx])
        r_inner = float(state.layer_edges[layer])
        r_outer = float(state.layer_edges[layer + 1])
        poly = adaptive_annular_sector_polygon(r_inner, r_outer, left, right, max_arc_step_deg)
        polygons.append(poly)
        ids.append(int(cell_ids[idx]))
        geometry[idx] = (r_inner, r_outer, left, right)
        centroids[idx] = polygon_centroid(poly)
    return polygons, ids, geometry, centroids


def add_direction_labels(ax: plt.Axes, radius: float) -> None:
    for angle, label, ha, va in [
        (math.pi / 2, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3 * math.pi / 2, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ]:
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=11, color="#444444")


def draw_edge_set(
    ax: plt.Axes,
    centers: np.ndarray,
    edges: list[tuple[int, int]],
    color: str,
    *,
    linestyle: str = "solid",
    linewidth: float = 0.6,
    alpha: float = 0.75,
    label: str | None = None,
) -> None:
    if not edges:
        return
    ax.add_collection(
        LineCollection(
            [(centers[i], centers[j]) for i, j in edges], colors=color, linestyles=linestyle, linewidths=linewidth, alpha=alpha, label=label
        )
    )


def save_scaffold(
    path: Path,
    centers: np.ndarray,
    layer_edges: np.ndarray,
    cell_ids: np.ndarray,
    original_edges: set[tuple[int, int]],
    layer_ids: np.ndarray,
    title: str,
) -> None:
    same, adjacent, long = classify_edges_by_layer(original_edges, layer_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=280)
    fig.patch.set_facecolor("white")
    circle = np.linspace(0.0, TAU, 721)
    for r in layer_edges:
        ax.plot(r * np.cos(circle), r * np.sin(circle), color="#dddddd", linewidth=0.45, zorder=1)
    draw_edge_set(ax, centers, same, "#2870b8", label=f"same-layer original ({len(same)})", linewidth=0.7, alpha=0.82)
    draw_edge_set(ax, centers, adjacent, "#c9483a", label=f"adjacent-layer original ({len(adjacent)})", linewidth=0.7, alpha=0.82)
    draw_edge_set(ax, centers, long, "#9e9e9e", linestyle="dashed", label=f"long-layer conflict ({len(long)})", linewidth=0.55, alpha=0.65)
    ax.scatter(centers[:, 0], centers[:, 1], s=7.0, color="#111111", linewidths=0, zorder=5)
    for idx, cell_id in enumerate(cell_ids):
        ax.text(
            float(centers[idx, 0]),
            float(centers[idx, 1]),
            str(int(cell_id)),
            ha="center",
            va="center",
            fontsize=2.9,
            color="#111111",
            zorder=6,
        )
    add_direction_labels(ax, float(np.max(layer_edges)) + 0.12)
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    setup_axis(ax, title, float(np.max(layer_edges)) + 0.22)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_partition(path: Path, state: PartitionState, polygons: list[np.ndarray], ids: list[int], max_arc_step_deg: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=300)
    fig.patch.set_facecolor("white")
    for pos, r in enumerate(state.layer_edges):
        pts = adaptive_circle_points(float(r), 0.0, TAU, max_arc_step_deg)
        is_outer = pos in (0, len(state.layer_edges) - 1)
        ax.plot(pts[:, 0], pts[:, 1], color="#222222" if is_outer else "#d7d7d7", linewidth=0.75 if is_outer else 0.45, zorder=1)
    ax.add_collection(PolyCollection(polygons, facecolors="none", edgecolors="#4b4b4b", linewidths=0.38, zorder=2))
    for poly, cell_id in zip(polygons, ids):
        xy = polygon_centroid(poly)
        ax.text(float(xy[0]), float(xy[1]), str(cell_id), ha="center", va="center", fontsize=3.0, color="#111111", zorder=3)
    add_direction_labels(ax, float(np.max(state.layer_edges)) + 0.12)
    setup_axis(ax, "Step 4 - Topology-aware Annular Partition v2", float(np.max(state.layer_edges)) + 0.22)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_final_topology_only(
    path: Path,
    centers: np.ndarray,
    cell_ids: np.ndarray,
    state: PartitionState,
    final_edges: set[tuple[int, int]],
    edge_kind: dict[tuple[int, int], str],
) -> None:
    same = [e for e in sorted(final_edges) if edge_kind.get(e) == "same_layer"]
    cross = [e for e in sorted(final_edges) if edge_kind.get(e) == "cross_layer"]
    fig, ax = plt.subplots(figsize=(10.0, 10.0), dpi=280)
    fig.patch.set_facecolor("white")
    draw_edge_set(ax, centers, same, "#355c9a", linewidth=0.58, alpha=0.72, label=f"same-layer final ({len(same)})")
    draw_edge_set(ax, centers, cross, "#8f5aa8", linewidth=0.52, alpha=0.58, label=f"cross-layer final ({len(cross)})")
    ax.scatter(centers[:, 0], centers[:, 1], s=8, color="#111111", zorder=4)
    for idx, cell_id in enumerate(cell_ids):
        ax.text(float(centers[idx, 0]), float(centers[idx, 1]), str(int(cell_id)), fontsize=2.8, ha="center", va="center", zorder=5)
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    setup_axis(ax, "Final Partition Adjacency", float(np.max(state.layer_edges)) + 0.22)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_adjacency_comparison(
    path: Path,
    centers: np.ndarray,
    cell_ids: np.ndarray,
    state: PartitionState,
    original_edges: set[tuple[int, int]],
    final_edges: set[tuple[int, int]],
) -> None:
    preserved = sorted(original_edges & final_edges)
    lost = sorted(original_edges - final_edges)
    new = sorted(final_edges - original_edges)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2), dpi=240)
    fig.patch.set_facecolor("white")
    panels = [
        ("Original adjacency", sorted(original_edges), "#555555", f"original edges = {len(original_edges)}"),
        ("Final adjacency", sorted(final_edges), "#355c9a", f"final edges = {len(final_edges)}"),
    ]
    for ax, (title, edges, color, note) in zip(axes[:2], panels):
        draw_edge_set(ax, centers, edges, color, linewidth=0.45, alpha=0.42)
        ax.scatter(centers[:, 0], centers[:, 1], s=5, color="#111111", zorder=3)
        ax.text(0.02, 0.98, note, transform=ax.transAxes, ha="left", va="top", fontsize=9)
        setup_axis(ax, title, float(np.max(state.layer_edges)) + 0.18)
    ax = axes[2]
    draw_edge_set(ax, centers, preserved, "#238443", linewidth=0.45, alpha=0.40, label=f"preserved ({len(preserved)})")
    draw_edge_set(ax, centers, lost, "#c9252d", linestyle="dashed", linewidth=0.5, alpha=0.38, label=f"lost ({len(lost)})")
    draw_edge_set(ax, centers, new, "#e08214", linestyle="dotted", linewidth=0.45, alpha=0.28, label=f"new ({len(new)})")
    ax.scatter(centers[:, 0], centers[:, 1], s=5, color="#111111", zorder=3)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8, frameon=False)
    setup_axis(ax, "Difference", float(np.max(state.layer_edges)) + 0.18)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_edges_only(
    path: Path,
    centers: np.ndarray,
    cell_ids: np.ndarray,
    state: PartitionState,
    edges: set[tuple[int, int]],
    title: str,
    color: str,
    linestyle: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 10.0), dpi=280)
    fig.patch.set_facecolor("white")
    draw_edge_set(ax, centers, sorted(edges), color, linestyle=linestyle, linewidth=0.64, alpha=0.74, label=f"{len(edges)} edges")
    ax.scatter(centers[:, 0], centers[:, 1], s=8, color="#111111", zorder=3)
    for idx, cell_id in enumerate(cell_ids):
        ax.text(float(centers[idx, 0]), float(centers[idx, 1]), str(int(cell_id)), fontsize=3.0, ha="center", va="center", zorder=4)
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    setup_axis(ax, title, float(np.max(state.layer_edges)) + 0.22)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_layer_capacity(path: Path, initial_counts: np.ndarray, optimized_counts: np.ndarray, targets: np.ndarray) -> None:
    x = np.arange(len(targets))
    fig, ax = plt.subplots(figsize=(9.4, 4.8), dpi=220)
    fig.patch.set_facecolor("white")
    ax.bar(x - 0.24, initial_counts, width=0.24, label="initial", color="#bdbdbd")
    ax.bar(x, optimized_counts, width=0.24, label="optimized", color="#4c78a8")
    ax.bar(x + 0.24, targets, width=0.24, label="target", color="#72b7b2")
    ax.set_xlabel("Layer")
    ax.set_ylabel("Cell count")
    ax.set_title("Layer Capacity Diagnostic", loc="left", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.legend(frameon=False, ncol=3)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.28)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_width_distribution(path: Path, widths: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(9.4, 4.8), dpi=220)
    fig.patch.set_facecolor("white")
    ax.hist(np.degrees(widths), bins=24, color="#4c78a8", alpha=0.86)
    ax.set_xlabel("Angular width (deg)")
    ax.set_ylabel("Cell count")
    ax.set_title("Angular Width Distribution", loc="left", fontsize=12, fontweight="bold")
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.28)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def edge_status_rows(original_edges: set[tuple[int, int]], final_edges: set[tuple[int, int]], cell_ids: np.ndarray) -> list[dict[str, Any]]:
    rows = []
    for i, j in sorted(original_edges | final_edges):
        in_original = (i, j) in original_edges
        in_final = (i, j) in final_edges
        status = "preserved" if in_original and in_final else "lost" if in_original else "new"
        rows.append(
            {
                "cell_id_a": int(cell_ids[i]),
                "cell_id_b": int(cell_ids[j]),
                "in_original": int(in_original),
                "in_final": int(in_final),
                "status": status,
            }
        )
    return rows


def assignment_rows(
    state: PartitionState,
    initial_layer_ids: np.ndarray,
    initial_orders: dict[int, list[int]],
    cell_ids: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    radius: np.ndarray,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    final_edges: set[tuple[int, int]],
    geometry: dict[int, tuple[float, float, float, float]],
) -> list[dict[str, Any]]:
    initial_slot = {int(idx): pos for order in initial_orders.values() for pos, idx in enumerate(order, start=1)}
    optimized_slot = {int(idx): pos for order in state.orders.values() for pos, idx in enumerate(order, start=1)}
    original_neighbors = build_neighbors(len(cell_ids), original_edges)
    preserved_neighbors = build_neighbors(len(cell_ids), original_edges & final_edges)
    lost_neighbors = build_neighbors(len(cell_ids), original_edges - final_edges)
    new_neighbors = build_neighbors(len(cell_ids), final_edges - original_edges)
    out = []
    for idx in range(len(cell_ids)):
        r_inner, r_outer, left, right = geometry[idx]
        center_theta = 0.5 * (left + right)
        out.append(
            {
                "cell_id": int(cell_ids[idx]),
                "row": int(rows[idx]),
                "col": int(cols[idx]),
                "initial_layer": int(initial_layer_ids[idx]),
                "optimized_layer": int(state.layer_ids[idx]),
                "layer_changed": int(initial_layer_ids[idx] != state.layer_ids[idx]),
                "initial_slot": int(initial_slot[idx]),
                "optimized_slot": int(optimized_slot[idx]),
                "source_radius": f"{float(radius[idx]):.10f}",
                "source_theta_deg": f"{math.degrees(float(theta[idx])):.6f}",
                "optimized_center_theta_deg": f"{math.degrees(float(normalize_angle(center_theta))):.6f}",
                "angular_displacement_deg": f"{math.degrees(circular_angle_distance(center_theta, float(theta[idx]))):.6f}",
                "theta_left_deg": f"{math.degrees(left):.6f}",
                "theta_right_deg": f"{math.degrees(right):.6f}",
                "angular_width_deg": f"{math.degrees(right - left):.6f}",
                "r_inner": f"{r_inner:.10f}",
                "r_outer": f"{r_outer:.10f}",
                "original_degree": int(len(original_neighbors[idx])),
                "preserved_degree": int(len(preserved_neighbors[idx])),
                "lost_degree": int(len(lost_neighbors[idx])),
                "new_degree": int(len(new_neighbors[idx])),
            }
        )
    return out


def geometry_validation(
    polygons: list[np.ndarray], geometry: dict[int, tuple[float, float, float, float]], state: PartitionState
) -> dict[str, Any]:
    max_outer_error = 0.0
    max_inner_error = 0.0
    negative = 0
    invalid = 0
    for idx, poly in enumerate(polygons):
        r_inner, r_outer, _left, _right = geometry[idx]
        radii = np.hypot(poly[:, 0], poly[:, 1])
        outer_count = len(poly) // 2
        outer_r = radii[:outer_count]
        inner_r = radii[outer_count:]
        max_outer_error = max(max_outer_error, float(np.max(np.abs(outer_r - r_outer))))
        max_inner_error = max(max_inner_error, float(np.max(np.abs(inner_r - r_inner))))
        area = polygon_area(poly)
        negative += int(area < -1e-12)
        invalid += int((not np.all(np.isfinite(poly))) or abs(area) < 1e-14)
    coverage_errors = {}
    for layer, b in state.boundaries.items():
        coverage_errors[str(layer)] = float(abs(np.sum(np.diff(b)) - TAU))
    return {
        "max_outer_arc_radial_error": max_outer_error,
        "max_inner_arc_radial_error": max_inner_error,
        "polygon_negative_area_count": negative,
        "invalid_polygon_count": invalid,
        "angular_coverage_error_per_layer": coverage_errors,
    }


def run_pipeline(args: argparse.Namespace, *, write_outputs: bool = True, sweep_label: str = "") -> dict[str, Any]:
    start = time.perf_counter()
    records = read_center_polar_cells(args.input_csv)
    cell_ids, rows, cols, radius, theta, original_centers = arrays_from_records(records)
    original_edges = read_topology_constraints(args.topology_csv, records)
    initial_layer_ids, layer_edges = assign_layers_by_radius(
        radius, args.layers, inner_radius=args.inner_radius, outer_radius=args.outer_radius
    )
    targets = target_layer_counts(len(cell_ids), layer_edges)
    seam_angle = float(
        normalize_angle(math.radians(args.seam_angle_deg) if args.seam_angle_deg is not None else choose_global_seam_angle(theta))
    )
    unwrapped_theta = unwrap_angles(theta, seam_angle)
    weights = V2Weights(
        w_lost=args.w_lost,
        w_new=args.w_new,
        w_layer_jump=args.w_layer_jump,
        w_capacity=args.w_capacity,
        w_angle=args.w_angle,
        w_angle_cap=args.w_angle_cap,
        w_rank=args.w_rank,
        w_width_balance=args.w_width_balance,
        w_aspect=args.w_aspect,
        w_degree_inflation=args.w_degree_inflation,
        w_nonedge_overlap=args.w_nonedge_overlap,
        w_radial_deviation=args.w_radial_deviation,
        w_radial_inversion=args.w_radial_inversion,
        w_empty_layer=args.w_empty_layer,
    )
    initial_orders = build_layer_orders(initial_layer_ids, unwrapped_theta, args.layers)
    initial_ranks = rank_maps(initial_orders)
    baseline_metrics_obj, baseline_raw, baseline_effective = baseline_v2(
        theta, initial_layer_ids, layer_edges, original_edges, initial_ranks, targets, weights, args
    )

    optimized_layer_ids, layer_summary = optimize_layer_assignment_v2(
        initial_layer_ids,
        radius,
        original_edges,
        targets,
        max_layer_shift=args.max_layer_shift,
        max_iters=args.layer_opt_iters,
        capacity_tolerance=args.capacity_tolerance,
        weights=weights,
    )
    initial_state = make_state_v2(
        "initial_v2", initial_layer_ids, layer_edges, unwrapped_theta, seam_angle, min_cell_angle=math.radians(args.min_cell_angle_deg)
    )
    evaluate_state(initial_state, theta, original_edges, initial_ranks, targets, weights, args)
    layer_state = make_state_v2(
        "layer_optimized_v2",
        optimized_layer_ids,
        layer_edges,
        unwrapped_theta,
        seam_angle,
        min_cell_angle=math.radians(args.min_cell_angle_deg),
    )
    evaluate_state(layer_state, theta, original_edges, initial_ranks, targets, weights, args)
    current_state = layer_state if layer_state.objective <= initial_state.objective + 1e-12 else initial_state
    states = [initial_state.copy(name="initial_v2"), current_state.copy(name="layer_selected_v2")]

    order_state, order_summary = optimize_orders_v2(
        current_state, unwrapped_theta, seam_angle, theta, original_edges, initial_ranks, targets, weights, args
    )
    if order_state.objective <= current_state.objective + 1e-12:
        current_state = order_state
    states.append(current_state.copy(name="order_selected_v2"))
    boundary_state, boundary_summary = optimize_boundaries_v2(current_state, theta, original_edges, initial_ranks, targets, weights, args)
    if boundary_state.objective <= current_state.objective + 1e-12:
        current_state = boundary_state
    states.append(current_state.copy(name="boundary_selected_v2"))

    scored = []
    stage_metrics: dict[str, Any] = {}
    for state in states:
        _obj, details = evaluate_state(state, theta, original_edges, initial_ranks, targets, weights, args)
        metrics = metrics_from_details(details)
        stage_metrics[state.name] = asdict(metrics)
        # Prefer F1, then fewer new edges, then precision, then objective.
        scored.append(
            (
                metrics.adjacency_f1,
                -metrics.new_edge_count,
                metrics.adjacency_precision,
                -metrics.objective,
                state.copy(name=state.name),
                details,
                metrics,
            )
        )
    scored.sort(reverse=True, key=lambda item: item[:4])
    selected_state = scored[0][4].copy(name=f"selected_{scored[0][4].name}")
    selected_details = scored[0][5]
    selected_metrics = scored[0][6]
    effective_edges = selected_details["effective_edges"]
    _raw2, _eff2, edge_kind = build_adjacencies_v2(
        selected_state.layer_ids, selected_state.orders, selected_state.boundaries, args.target_overlap_ratio
    )

    polygons, ids, geometry, optimized_centers = polygons_from_state_v2(selected_state, cell_ids, args.max_arc_step_deg)
    geom = geometry_validation(polygons, geometry, selected_state)
    runtime = time.perf_counter() - start

    if write_outputs:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        save_scaffold(
            args.output_dir / "step03a_initial_topology_scaffold.png",
            original_centers,
            layer_edges,
            cell_ids,
            original_edges,
            initial_layer_ids,
            "Step 3A - Initial Original-topology Scaffold",
        )
        save_scaffold(
            args.output_dir / "step03b_optimized_topology_scaffold.png",
            optimized_centers,
            layer_edges,
            cell_ids,
            original_edges,
            selected_state.layer_ids,
            "Step 3B - Optimized Original-topology Scaffold",
        )
        save_partition(args.output_dir / "step04_topology_aware_partition_v2.png", selected_state, polygons, ids, args.max_arc_step_deg)
        save_final_topology_only(
            args.output_dir / "step04_final_topology_only.png", optimized_centers, cell_ids, selected_state, effective_edges, edge_kind
        )
        save_adjacency_comparison(
            args.output_dir / "step04_adjacency_comparison.png",
            optimized_centers,
            cell_ids,
            selected_state,
            original_edges,
            effective_edges,
        )
        save_edges_only(
            args.output_dir / "step04_lost_edges_only.png",
            optimized_centers,
            cell_ids,
            selected_state,
            original_edges - effective_edges,
            "Lost Original Edges Only",
            "#c9252d",
            "dashed",
        )
        save_edges_only(
            args.output_dir / "step04_new_edges_only.png",
            optimized_centers,
            cell_ids,
            selected_state,
            effective_edges - original_edges,
            "New Final Edges Only",
            "#e08214",
            "dotted",
        )
        save_layer_capacity(
            args.output_dir / "layer_capacity_diagnostic.png",
            layer_counts(initial_layer_ids, args.layers),
            layer_counts(selected_state.layer_ids, args.layers),
            targets,
        )
        save_width_distribution(args.output_dir / "angular_width_distribution.png", selected_details["width"]["widths"])
        write_csv(
            args.output_dir / "step04_topology_edge_status_v2.csv",
            edge_status_rows(original_edges, effective_edges, cell_ids),
            ["cell_id_a", "cell_id_b", "in_original", "in_final", "status"],
        )
        write_csv(
            args.output_dir / "step04_topology_aware_assignment_v2.csv",
            assignment_rows(
                selected_state,
                initial_layer_ids,
                initial_orders,
                cell_ids,
                rows,
                cols,
                radius,
                theta,
                original_edges,
                effective_edges,
                geometry,
            ),
            [
                "cell_id",
                "row",
                "col",
                "initial_layer",
                "optimized_layer",
                "layer_changed",
                "initial_slot",
                "optimized_slot",
                "source_radius",
                "source_theta_deg",
                "optimized_center_theta_deg",
                "angular_displacement_deg",
                "theta_left_deg",
                "theta_right_deg",
                "angular_width_deg",
                "r_inner",
                "r_outer",
                "original_degree",
                "preserved_degree",
                "lost_degree",
                "new_degree",
            ],
        )

    degree = selected_details["degree"]
    max_inflation = int(np.max(degree["inflation"]))
    max_inflation_cells = [int(cell_ids[i]) for i in np.flatnonzero(degree["inflation"] == max_inflation)]
    payload = {
        "parameters": {k: portable_parameter(v) for k, v in vars(args).items()},
        "sweep_label": sweep_label,
        "runtime_seconds": runtime,
        "adjacency_metric_primary": "effective_adjacency",
        "effective_adjacency_definition": "same-layer cyclic contacts plus adjacent-layer contacts whose angular overlap / min(cell widths) >= target_overlap_ratio",
        "raw_contact_adjacency_definition": "same-layer cyclic contacts plus any adjacent-layer positive angular overlap above numerical tolerance",
        "seam_angle_deg": math.degrees(seam_angle),
        "initial_layer_counts": layer_counts(initial_layer_ids, args.layers).tolist(),
        "optimized_layer_counts": layer_counts(selected_state.layer_ids, args.layers).tolist(),
        "target_layer_counts": targets.tolist(),
        "layer_capacity_deviation": ((layer_counts(selected_state.layer_ids, args.layers) - targets) / np.maximum(1, targets)).tolist(),
        "baseline_effective_metrics": asdict(baseline_metrics_obj),
        "optimized_effective_metrics": asdict(selected_metrics),
        "stage_metrics": stage_metrics,
        "layer_optimization": layer_summary,
        "order_optimization": order_summary,
        "boundary_optimization": boundary_summary,
        "geometry_validation": geom,
        "lost_edges": [{"cell_id_a": int(cell_ids[i]), "cell_id_b": int(cell_ids[j])} for i, j in sorted(original_edges - effective_edges)],
        "new_edges": [{"cell_id_a": int(cell_ids[i]), "cell_id_b": int(cell_ids[j])} for i, j in sorted(effective_edges - original_edges)],
        "remaining_long_layer_conflicts": [
            {
                "cell_id_a": int(cell_ids[i]),
                "cell_id_b": int(cell_ids[j]),
                "layer_delta": int(abs(selected_state.layer_ids[i] - selected_state.layer_ids[j])),
            }
            for i, j in sorted(original_edges)
            if abs(int(selected_state.layer_ids[i]) - int(selected_state.layer_ids[j])) > 1
        ],
        "max_degree_inflation": max_inflation,
        "max_degree_inflation_cells": max_inflation_cells,
        "layer_changed_cells": [
            {"cell_id": int(cell_ids[i]), "initial_layer": int(initial_layer_ids[i]), "optimized_layer": int(selected_state.layer_ids[i])}
            for i in np.flatnonzero(initial_layer_ids != selected_state.layer_ids)
        ],
    }
    if write_outputs:
        (args.output_dir / "step04_topology_aware_metrics_v2.json").write_text(
            json.dumps(json_safe(payload), ensure_ascii=False, indent=2), encoding="utf-8"
        )
    return payload


def run_parameter_sweep(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    base = argparse.Namespace(**vars(args))
    base.layer_opt_iters = min(args.layer_opt_iters, 18)
    base.order_opt_iters = min(args.order_opt_iters, 10)
    base.boundary_opt_iters = min(args.boundary_opt_iters, 8)
    for w_new in [3.0, 4.0, 5.0]:
        for w_capacity in [2.0, 4.0, 6.0]:
            for w_degree in [1.0, 2.0, 3.0]:
                for max_factor in [2.0, 2.5, 3.0]:
                    trial = argparse.Namespace(**vars(base))
                    trial.w_new = w_new
                    trial.w_capacity = w_capacity
                    trial.w_degree_inflation = w_degree
                    trial.max_angle_factor = max_factor
                    payload = run_pipeline(
                        trial, write_outputs=False, sweep_label=f"new{w_new}_cap{w_capacity}_deg{w_degree}_max{max_factor}"
                    )
                    m = payload["optimized_effective_metrics"]
                    rows.append(
                        {
                            "w_new": w_new,
                            "w_capacity": w_capacity,
                            "w_degree_inflation": w_degree,
                            "max_angle_factor": max_factor,
                            "precision": m["adjacency_precision"],
                            "recall": m["adjacency_recall"],
                            "f1": m["adjacency_f1"],
                            "lost": m["lost_edge_count"],
                            "new": m["new_edge_count"],
                            "maximum_angular_width": m["max_angular_width_deg"],
                            "minimum_angular_width": m["min_angular_width_deg"],
                            "maximum_degree_inflation": m["max_degree_inflation"],
                            "capacity_deviation": m["max_capacity_deviation"],
                            "mean_angular_error": m["mean_angular_error_deg"],
                            "maximum_angular_error": m["max_angular_error_deg"],
                            "runtime": payload["runtime_seconds"],
                        }
                    )
    fieldnames = list(rows[0].keys()) if rows else []
    write_csv(args.output_dir / "parameter_sweep_results.csv", rows, fieldnames)
    if rows:
        fig, ax = plt.subplots(figsize=(8.0, 5.4), dpi=220)
        fig.patch.set_facecolor("white")
        sc = ax.scatter([r["new"] for r in rows], [r["f1"] for r in rows], c=[r["precision"] for r in rows], cmap="viridis", s=38)
        ax.set_xlabel("New edges")
        ax.set_ylabel("F1")
        ax.set_title("Parameter Sweep Pareto View", loc="left", fontsize=12, fontweight="bold")
        fig.colorbar(sc, ax=ax, label="Precision")
        ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.25)
        fig.savefig(args.output_dir / "parameter_sweep_pareto.png", facecolor="white", bbox_inches="tight")
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step 4 v2: capacity-balanced, geometry-aware topology-aware annular partition.")
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--topology-csv", type=Path, default=DEFAULT_TOPOLOGY_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--layers", type=int, default=6)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=1.85)
    parser.add_argument("--seam-angle-deg", type=float, default=None)
    parser.add_argument("--max-layer-shift", type=int, default=1)
    parser.add_argument("--layer-opt-iters", type=int, default=30)
    parser.add_argument("--capacity-tolerance", type=float, default=0.25)
    parser.add_argument("--w-lost", type=float, default=4.0)
    parser.add_argument("--w-new", type=float, default=4.5)
    parser.add_argument("--w-layer-jump", type=float, default=4.0)
    parser.add_argument("--w-capacity", type=float, default=4.0)
    parser.add_argument("--w-angle", type=float, default=1.2)
    parser.add_argument("--w-angle-cap", type=float, default=3.0)
    parser.add_argument("--w-rank", type=float, default=0.3)
    parser.add_argument("--w-width-balance", type=float, default=2.0)
    parser.add_argument("--w-aspect", type=float, default=1.5)
    parser.add_argument("--w-degree-inflation", type=float, default=2.0)
    parser.add_argument("--w-nonedge-overlap", type=float, default=3.0)
    parser.add_argument("--w-radial-deviation", type=float, default=1.0)
    parser.add_argument("--w-radial-inversion", type=float, default=1.0)
    parser.add_argument("--w-empty-layer", type=float, default=10000.0)
    parser.add_argument("--order-opt-iters", type=int, default=28)
    parser.add_argument("--max-rank-shift", type=int, default=2)
    parser.add_argument("--boundary-opt-iters", type=int, default=18)
    parser.add_argument("--boundary-search-samples", type=int, default=15)
    parser.add_argument("--min-cell-angle-deg", type=float, default=1.5)
    parser.add_argument("--max-angle-factor", type=float, default=2.5)
    parser.add_argument("--min-aspect", type=float, default=0.35)
    parser.add_argument("--max-aspect", type=float, default=4.0)
    parser.add_argument("--target-overlap-ratio", type=float, default=0.1)
    parser.add_argument("--max-angle-displacement-deg", type=float, default=15.0)
    parser.add_argument("--max-arc-step-deg", type=float, default=0.5)
    parser.add_argument("--run-parameter-sweep", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.outer_radius <= args.inner_radius:
        raise ValueError("--outer-radius must be greater than --inner-radius.")
    if args.max_arc_step_deg <= 0:
        raise ValueError("--max-arc-step-deg must be positive.")
    payload = run_pipeline(args, write_outputs=True)
    if args.run_parameter_sweep:
        run_parameter_sweep(args)
    m = payload["optimized_effective_metrics"]
    b = payload["baseline_effective_metrics"]
    print("\nBaseline effective adjacency:")
    print(f"- Precision: {b['adjacency_precision']:.4f}")
    print(f"- Recall: {b['adjacency_recall']:.4f}")
    print(f"- F1: {b['adjacency_f1']:.4f}")
    print(f"- Lost: {b['lost_edge_count']}")
    print(f"- New: {b['new_edge_count']}")
    print("\nTopology-aware v2 effective adjacency:")
    print(f"- Precision: {m['adjacency_precision']:.4f}")
    print(f"- Recall: {m['adjacency_recall']:.4f}")
    print(f"- F1: {m['adjacency_f1']:.4f}")
    print(f"- Lost: {m['lost_edge_count']}")
    print(f"- New: {m['new_edge_count']}")
    print(f"- Layer counts: {payload['optimized_layer_counts']} target={payload['target_layer_counts']}")
    print(f"- Angular width min/max: {m['min_angular_width_deg']:.3f} / {m['max_angular_width_deg']:.3f} deg")
    print(
        f"- Mean/median/max angular error: {m['mean_angular_error_deg']:.3f} / {m['median_angular_error_deg']:.3f} / {m['max_angular_error_deg']:.3f} deg"
    )
    print(f"- Max degree inflation: {m['max_degree_inflation']} cells={payload['max_degree_inflation_cells']}")
    print("\nOutputs:")
    for name in [
        "step03a_initial_topology_scaffold.png",
        "step03b_optimized_topology_scaffold.png",
        "step04_topology_aware_partition_v2.png",
        "step04_final_topology_only.png",
        "step04_adjacency_comparison.png",
        "step04_lost_edges_only.png",
        "step04_new_edges_only.png",
        "step04_topology_aware_assignment_v2.csv",
        "step04_topology_aware_metrics_v2.json",
        "step04_topology_edge_status_v2.csv",
        "layer_capacity_diagnostic.png",
        "angular_width_distribution.png",
    ]:
        print(f"- {args.output_dir / name}")


if __name__ == "__main__":
    main()
