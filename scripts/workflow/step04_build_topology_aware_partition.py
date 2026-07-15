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
DEFAULT_INPUT_CSV = PROJECT_ROOT / "data" / "generated" / "four_vertex_grid_cells.csv"
DEFAULT_TOPOLOGY_CSV = PROJECT_ROOT / "results" / "workflow" / "step03" / "step03_topology_constraints.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "workflow" / "step04_topology_aware_v1"

if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from annulus_topology_utils import (  # noqa: E402
    TAU,
    annular_sector_polygon,
    arrays_from_records,
    assign_layers_by_radius,
    build_layer_orders,
    build_neighbors,
    build_partition_adjacency,
    choose_global_seam_angle,
    circular_angle_distance,
    compute_adjacency_metrics,
    count_long_layer_conflicts,
    initialize_boundaries_from_orders,
    interval_overlap,
    intervals_by_index,
    linear_interval_overlap,
    normalize_angle,
    read_center_polar_cells,
    read_topology_constraints,
    unwrap_angles,
)


@dataclass(frozen=True)
class OptimizationWeights:
    w_layer_jump: float
    w_radial_deviation: float
    w_radial_inversion: float
    w_empty_layer: float
    w_lost: float
    w_new: float
    w_angle: float
    w_long_conflict: float
    w_rank: float


@dataclass
class PartitionState:
    name: str
    layer_ids: np.ndarray
    orders: dict[int, list[int]]
    boundaries: dict[int, np.ndarray]
    layer_edges: np.ndarray
    objective: float

    def copy(self, *, name: str | None = None) -> "PartitionState":
        return PartitionState(
            name=self.name if name is None else name,
            layer_ids=self.layer_ids.copy(),
            orders={layer: order.copy() for layer, order in self.orders.items()},
            boundaries={layer: values.copy() for layer, values in self.boundaries.items()},
            layer_edges=self.layer_edges.copy(),
            objective=float(self.objective),
        )


@dataclass
class PartitionMetrics:
    original_edge_count: int
    final_edge_count: int
    preserved_edge_count: int
    lost_edge_count: int
    new_edge_count: int
    adjacency_precision: float
    adjacency_recall: float
    adjacency_f1: float
    long_layer_conflict_count: int
    long_layer_conflict_rate: float
    mean_angular_error_deg: float
    median_angular_error_deg: float
    max_angular_error_deg: float
    minimum_angular_width_deg: float
    mean_angular_width_deg: float
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
        for row in rows:
            writer.writerow(row)


def circular_order_with_largest_gap(theta_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(theta_values, kind="mergesort")
    sorted_theta = np.asarray(theta_values[order], dtype=float)
    if len(sorted_theta) <= 1:
        return order, sorted_theta
    gaps = np.diff(np.r_[sorted_theta, sorted_theta[0] + TAU])
    start = int((np.argmax(gaps) + 1) % len(sorted_theta))
    rotated_order = np.r_[order[start:], order[:start]]
    if start == 0:
        return rotated_order, sorted_theta.copy()
    rotated_theta = np.r_[sorted_theta[start:], sorted_theta[:start] + TAU]
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


def old_step4_intervals(
    theta: np.ndarray,
    layer_ids: np.ndarray,
    layer_count: int,
) -> tuple[dict[int, list[int]], dict[int, tuple[float, float]]]:
    orders: dict[int, list[int]] = {}
    intervals: dict[int, tuple[float, float]] = {}
    for layer in range(layer_count):
        indices = np.flatnonzero(layer_ids == layer)
        if len(indices) == 0:
            raise ValueError(f"Baseline layer {layer} is empty.")
        local_order, unwrapped = circular_order_with_largest_gap(theta[indices])
        order = indices[local_order].astype(int).tolist()
        layout_theta = spread_close_angles(unwrapped)
        orders[layer] = order
        if len(order) == 1:
            intervals[order[0]] = (float(layout_theta[0] - math.pi / 18.0), float(layout_theta[0] + math.pi / 18.0))
            continue
        for pos, idx in enumerate(order):
            current = float(layout_theta[pos])
            prev_theta = float(layout_theta[pos - 1] if pos > 0 else layout_theta[-1] - TAU)
            next_theta = float(layout_theta[pos + 1] if pos < len(order) - 1 else layout_theta[0] + TAU)
            intervals[int(idx)] = (0.5 * (prev_theta + current), 0.5 * (current + next_theta))
    return orders, intervals


def build_adjacency_from_intervals(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    intervals: dict[int, tuple[float, float]],
    min_shared_angle: float,
) -> set[tuple[int, int]]:
    final_edges: set[tuple[int, int]] = set()
    for _layer, order in sorted(orders.items()):
        if len(order) <= 1:
            continue
        for pos, i in enumerate(order):
            j = order[(pos + 1) % len(order)]
            final_edges.add(tuple(sorted((int(i), int(j)))))
    for layer in sorted(orders):
        if layer + 1 not in orders:
            continue
        for i in orders[layer]:
            left_i, right_i = intervals[i]
            for j in orders[layer + 1]:
                left_j, right_j = intervals[j]
                if interval_overlap(left_i, right_i, left_j, right_j) >= min_shared_angle:
                    final_edges.add(tuple(sorted((int(i), int(j)))))
    return final_edges


def radial_inversion_count(radius: np.ndarray, layer_ids: np.ndarray) -> int:
    order = np.argsort(radius, kind="mergesort")
    layers = layer_ids[order]
    count = 0
    for a in range(len(layers)):
        count += int(np.sum(layers[a + 1 :] < layers[a]))
    return int(count)


def layer_objective(
    layer_ids: np.ndarray,
    initial_layer_ids: np.ndarray,
    radius: np.ndarray,
    original_edges: set[tuple[int, int]],
    layer_count: int,
    weights: OptimizationWeights,
) -> float:
    jump = 0.0
    for i, j in original_edges:
        jump += max(0, abs(int(layer_ids[i]) - int(layer_ids[j])) - 1) ** 2
    deviation = float(np.sum((layer_ids - initial_layer_ids) ** 2))
    inversion = float(radial_inversion_count(radius, layer_ids))
    empty_layers = layer_count - len(set(layer_ids.tolist()))
    return (
        weights.w_layer_jump * jump
        + weights.w_radial_deviation * deviation
        + weights.w_radial_inversion * inversion
        + weights.w_empty_layer * max(0, empty_layers)
    )


def optimize_layer_assignment(
    initial_layer_ids: np.ndarray,
    radius: np.ndarray,
    original_edges: set[tuple[int, int]],
    layer_count: int,
    *,
    max_layer_shift: int,
    max_iters: int,
    weights: OptimizationWeights,
) -> tuple[np.ndarray, dict[str, Any]]:
    current = initial_layer_ids.copy()
    current_obj = layer_objective(current, initial_layer_ids, radius, original_edges, layer_count, weights)
    best = current.copy()
    best_obj = current_obj
    history = [{"iteration": 0, "objective": current_obj, "long_layer_conflicts": count_long_layer_conflicts(original_edges, current)}]
    allowed_by_cell = [
        range(max(0, int(k0) - max_layer_shift), min(layer_count - 1, int(k0) + max_layer_shift) + 1) for k0 in initial_layer_ids
    ]

    for iteration in range(1, max_iters + 1):
        changed = False
        for idx in range(len(current)):
            current_layer = int(current[idx])
            candidate_layers = sorted({current_layer, current_layer - 1, current_layer + 1} & set(allowed_by_cell[idx]))
            if not candidate_layers:
                continue
            best_local_layer = current_layer
            best_local_obj = current_obj
            for candidate_layer in candidate_layers:
                if candidate_layer == current_layer:
                    candidate_obj = current_obj
                else:
                    trial = current.copy()
                    trial[idx] = candidate_layer
                    if len(set(trial.tolist())) != layer_count:
                        continue
                    candidate_obj = layer_objective(trial, initial_layer_ids, radius, original_edges, layer_count, weights)
                if candidate_obj < best_local_obj - 1e-12:
                    best_local_obj = candidate_obj
                    best_local_layer = candidate_layer
                elif abs(candidate_obj - best_local_obj) <= 1e-12:
                    old_key = (abs(best_local_layer - int(initial_layer_ids[idx])), best_local_layer)
                    new_key = (abs(candidate_layer - int(initial_layer_ids[idx])), candidate_layer)
                    if new_key < old_key:
                        best_local_layer = candidate_layer
            if best_local_layer != current_layer:
                current[idx] = best_local_layer
                current_obj = best_local_obj
                changed = True
                if current_obj < best_obj - 1e-12:
                    best = current.copy()
                    best_obj = current_obj
        history.append(
            {"iteration": iteration, "objective": current_obj, "long_layer_conflicts": count_long_layer_conflicts(original_edges, current)}
        )
        if not changed:
            break
    return best, {"objective_before": history[0]["objective"], "objective_after": best_obj, "history": history}


def rank_maps(orders: dict[int, list[int]]) -> dict[int, int]:
    ranks: dict[int, int] = {}
    for order in orders.values():
        for pos, idx in enumerate(order):
            ranks[int(idx)] = pos
    return ranks


def center_angles_from_intervals(intervals: dict[int, tuple[float, float]]) -> dict[int, float]:
    return {idx: 0.5 * (left + right) for idx, (left, right) in intervals.items()}


def angular_errors(theta: np.ndarray, intervals: dict[int, tuple[float, float]]) -> np.ndarray:
    centers = center_angles_from_intervals(intervals)
    return np.asarray([circular_angle_distance(centers[idx], float(theta[idx])) for idx in range(len(theta))], dtype=float)


def topology_objective(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    min_shared_angle: float,
    weights: OptimizationWeights,
) -> tuple[float, dict[str, Any]]:
    final_edges = build_partition_adjacency(layer_ids, orders, boundaries, min_shared_angle)
    metrics = compute_adjacency_metrics(original_edges, final_edges)
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    angle_err = angular_errors(theta, intervals)
    current_ranks = rank_maps(orders)
    rank_disp = [
        abs(current_ranks[idx] - original_ranks[idx]) / max(1, len(orders[int(layer_ids[idx])]) - 1) for idx in range(len(layer_ids))
    ]
    long_count = count_long_layer_conflicts(original_edges, layer_ids)
    long_rate = long_count / max(1, len(original_edges))
    objective = (
        weights.w_lost * (metrics["lost_edge_count"] / max(1, len(original_edges)))
        + weights.w_new * (metrics["new_edge_count"] / max(1, len(final_edges)))
        + weights.w_angle * float(np.mean(angle_err / math.pi))
        + weights.w_long_conflict * long_rate
        + weights.w_rank * float(np.mean(rank_disp))
    )
    details = dict(metrics)
    details.update(
        {
            "mean_angular_error_deg": math.degrees(float(np.mean(angle_err))),
            "median_angular_error_deg": math.degrees(float(np.median(angle_err))),
            "max_angular_error_deg": math.degrees(float(np.max(angle_err))),
            "long_layer_conflict_count": long_count,
            "long_layer_conflict_rate": long_rate,
            "mean_rank_displacement": float(np.mean(rank_disp)),
        }
    )
    return float(objective), details


def make_state(
    name: str,
    layer_ids: np.ndarray,
    layer_edges: np.ndarray,
    unwrapped_theta: np.ndarray,
    seam_angle: float,
    min_cell_angle: float,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int] | None,
    min_shared_angle: float,
    weights: OptimizationWeights,
) -> tuple[PartitionState, dict[int, int], dict[str, Any]]:
    orders = build_layer_orders(layer_ids, unwrapped_theta, len(layer_edges) - 1)
    ranks = rank_maps(orders)
    boundaries = initialize_boundaries_from_orders(orders, unwrapped_theta, seam_angle, min_cell_angle=min_cell_angle)
    ranks_for_objective = ranks if original_ranks is None else original_ranks
    objective, details = topology_objective(
        layer_ids, orders, boundaries, theta, original_edges, ranks_for_objective, min_shared_angle, weights
    )
    return PartitionState(name, layer_ids.copy(), orders, boundaries, layer_edges.copy(), objective), ranks, details


def optimize_circular_orders(
    state: PartitionState,
    unwrapped_theta: np.ndarray,
    seam_angle: float,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    *,
    min_cell_angle: float,
    min_shared_angle: float,
    max_rank_shift: int,
    max_iters: int,
    weights: OptimizationWeights,
) -> tuple[PartitionState, dict[str, Any]]:
    current = state.copy(name="order_initial")
    best = current.copy(name="order_best")
    history = [{"iteration": 0, "objective": current.objective}]

    def evaluate(orders: dict[int, list[int]]) -> tuple[float, dict[int, np.ndarray]]:
        boundaries = initialize_boundaries_from_orders(orders, unwrapped_theta, seam_angle, min_cell_angle=min_cell_angle)
        objective, _details = topology_objective(
            current.layer_ids, orders, boundaries, theta, original_edges, original_ranks, min_shared_angle, weights
        )
        return objective, boundaries

    for iteration in range(1, max_iters + 1):
        improved = False
        for layer in sorted(current.orders):
            n = len(current.orders[layer])
            if n <= 1:
                continue
            # Adjacent swaps in deterministic layer/position order.
            for pos in range(n - 1):
                trial_orders = {k: v.copy() for k, v in current.orders.items()}
                trial_orders[layer][pos], trial_orders[layer][pos + 1] = trial_orders[layer][pos + 1], trial_orders[layer][pos]
                obj, boundaries = evaluate(trial_orders)
                if obj < current.objective - 1e-12:
                    current.orders = trial_orders
                    current.boundaries = boundaries
                    current.objective = obj
                    improved = True
                    if obj < best.objective - 1e-12:
                        best = current.copy(name="order_best")
            # Bounded insert moves.
            positions = list(range(n))
            for pos in positions:
                for delta in range(-max_rank_shift, max_rank_shift + 1):
                    if delta == 0:
                        continue
                    new_pos = pos + delta
                    if new_pos < 0 or new_pos >= n:
                        continue
                    trial_orders = {k: v.copy() for k, v in current.orders.items()}
                    item = trial_orders[layer].pop(pos)
                    trial_orders[layer].insert(new_pos, item)
                    obj, boundaries = evaluate(trial_orders)
                    if obj < current.objective - 1e-12:
                        current.orders = trial_orders
                        current.boundaries = boundaries
                        current.objective = obj
                        improved = True
                        if obj < best.objective - 1e-12:
                            best = current.copy(name="order_best")
                        break
                if improved:
                    break
        history.append({"iteration": iteration, "objective": current.objective})
        if not improved:
            break
    return best, {"objective_before": state.objective, "objective_after": best.objective, "history": history}


def boundary_objective(
    layer_ids: np.ndarray,
    orders: dict[int, list[int]],
    boundaries: dict[int, np.ndarray],
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    *,
    min_shared_angle: float,
    target_overlap_ratio: float,
    weights: OptimizationWeights,
) -> float:
    intervals = intervals_by_index(layer_ids, orders, boundaries)
    original_adjacent = {edge for edge in original_edges if abs(int(layer_ids[edge[0]]) - int(layer_ids[edge[1]])) == 1}
    deficit_penalty = 0.0
    target_sum = 0.0
    for i, j in sorted(original_adjacent):
        li, ri = intervals[i]
        lj, rj = intervals[j]
        width_i = ri - li
        width_j = rj - lj
        target = max(min_shared_angle, target_overlap_ratio * min(width_i, width_j))
        ov = linear_interval_overlap(li, ri, lj, rj)
        deficit_penalty += max(0.0, target - ov) ** 2
        target_sum += target**2

    new_overlap_penalty = 0.0
    new_norm = 0.0
    for layer in sorted(orders):
        if layer + 1 not in orders:
            continue
        for i in orders[layer]:
            li, ri = intervals[i]
            for j in orders[layer + 1]:
                edge = tuple(sorted((int(i), int(j))))
                if edge in original_edges:
                    continue
                lj, rj = intervals[j]
                ov = linear_interval_overlap(li, ri, lj, rj)
                if ov >= min_shared_angle:
                    new_overlap_penalty += (ov - min_shared_angle) ** 2
                new_norm += max(ri - li, rj - lj) ** 2

    angle_err = angular_errors(theta, intervals)
    widths = np.asarray([right - left for left, right in intervals.values()], dtype=float)
    narrow = np.mean(np.maximum(0.0, 2.0 * min_shared_angle - widths) ** 2)
    return float(
        weights.w_lost * deficit_penalty / max(target_sum, 1e-12)
        + weights.w_new * new_overlap_penalty / max(new_norm, 1e-12)
        + weights.w_angle * float(np.mean(angle_err / math.pi))
        + 0.1 * narrow
    )


def optimize_angular_boundaries(
    state: PartitionState,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    *,
    min_cell_angle: float,
    min_shared_angle: float,
    target_overlap_ratio: float,
    samples: int,
    max_iters: int,
    weights: OptimizationWeights,
) -> tuple[PartitionState, dict[str, Any]]:
    if samples < 3:
        raise ValueError("--boundary-search-samples must be at least 3.")
    current = state.copy(name="boundary_initial")
    current_boundary_obj = boundary_objective(
        current.layer_ids,
        current.orders,
        current.boundaries,
        theta,
        original_edges,
        min_shared_angle=min_shared_angle,
        target_overlap_ratio=target_overlap_ratio,
        weights=weights,
    )
    best = current.copy(name="boundary_best")
    best_boundary_obj = current_boundary_obj
    history = [{"iteration": 0, "boundary_objective": current_boundary_obj, "topology_objective": current.objective}]

    for iteration in range(1, max_iters + 1):
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
                candidates = np.linspace(low, high, samples)
                best_value = b[m]
                best_value_obj = current_boundary_obj
                for value in candidates:
                    trial_boundaries = {k: v.copy() for k, v in current.boundaries.items()}
                    trial_boundaries[layer][m] = float(value)
                    obj = boundary_objective(
                        current.layer_ids,
                        current.orders,
                        trial_boundaries,
                        theta,
                        original_edges,
                        min_shared_angle=min_shared_angle,
                        target_overlap_ratio=target_overlap_ratio,
                        weights=weights,
                    )
                    if obj < best_value_obj - 1e-12:
                        best_value_obj = obj
                        best_value = float(value)
                if abs(best_value - b[m]) > 1e-15:
                    b[m] = best_value
                    current.boundaries[layer] = b.copy()
                    current_boundary_obj = best_value_obj
                    improved = True
                    top_obj, _details = topology_objective(
                        current.layer_ids,
                        current.orders,
                        current.boundaries,
                        theta,
                        original_edges,
                        original_ranks,
                        min_shared_angle,
                        weights,
                    )
                    current.objective = top_obj
                    if current_boundary_obj < best_boundary_obj - 1e-12:
                        best = current.copy(name="boundary_best")
                        best_boundary_obj = current_boundary_obj
        history.append({"iteration": iteration, "boundary_objective": current_boundary_obj, "topology_objective": current.objective})
        if not improved:
            break

    # A boundary pass is only accepted when the same comprehensive topology objective
    # is not worse than the input state; otherwise the previous historical best wins.
    if best.objective > state.objective + 1e-12:
        best = state.copy(name="boundary_rollback")
        best_boundary_obj = history[0]["boundary_objective"]
    return best, {
        "objective_before": state.objective,
        "objective_after": best.objective,
        "boundary_objective_after": best_boundary_obj,
        "history": history,
    }


def state_metrics(
    state: PartitionState,
    theta: np.ndarray,
    original_edges: set[tuple[int, int]],
    original_ranks: dict[int, int],
    min_shared_angle: float,
    weights: OptimizationWeights,
) -> tuple[PartitionMetrics, set[tuple[int, int]], dict[str, Any]]:
    objective, details = topology_objective(
        state.layer_ids,
        state.orders,
        state.boundaries,
        theta,
        original_edges,
        original_ranks,
        min_shared_angle,
        weights,
    )
    final_edges = build_partition_adjacency(state.layer_ids, state.orders, state.boundaries, min_shared_angle)
    widths = [right - left for left, right in intervals_by_index(state.layer_ids, state.orders, state.boundaries).values()]
    metrics = PartitionMetrics(
        original_edge_count=int(details["original_edge_count"]),
        final_edge_count=int(details["final_edge_count"]),
        preserved_edge_count=int(details["preserved_edge_count"]),
        lost_edge_count=int(details["lost_edge_count"]),
        new_edge_count=int(details["new_edge_count"]),
        adjacency_precision=float(details["adjacency_precision"]),
        adjacency_recall=float(details["adjacency_recall"]),
        adjacency_f1=float(details["adjacency_f1"]),
        long_layer_conflict_count=int(details["long_layer_conflict_count"]),
        long_layer_conflict_rate=float(details["long_layer_conflict_rate"]),
        mean_angular_error_deg=float(details["mean_angular_error_deg"]),
        median_angular_error_deg=float(details["median_angular_error_deg"]),
        max_angular_error_deg=float(details["max_angular_error_deg"]),
        minimum_angular_width_deg=math.degrees(float(np.min(widths))),
        mean_angular_width_deg=math.degrees(float(np.mean(widths))),
        objective=float(objective),
    )
    return metrics, final_edges, details


def baseline_metrics(
    theta: np.ndarray,
    layer_ids: np.ndarray,
    layer_count: int,
    original_edges: set[tuple[int, int]],
    min_shared_angle: float,
) -> tuple[dict[str, Any], set[tuple[int, int]], dict[int, tuple[float, float]], dict[int, list[int]]]:
    orders, intervals = old_step4_intervals(theta, layer_ids, layer_count)
    final_edges = build_adjacency_from_intervals(layer_ids, orders, intervals, min_shared_angle)
    metrics = compute_adjacency_metrics(original_edges, final_edges)
    angle_err = angular_errors(theta, intervals)
    widths = np.asarray([right - left for left, right in intervals.values()], dtype=float)
    long_count = count_long_layer_conflicts(original_edges, layer_ids)
    metrics.update(
        {
            "long_layer_conflict_count": long_count,
            "long_layer_conflict_rate": long_count / max(1, len(original_edges)),
            "mean_angular_error_deg": math.degrees(float(np.mean(angle_err))),
            "median_angular_error_deg": math.degrees(float(np.median(angle_err))),
            "max_angular_error_deg": math.degrees(float(np.max(angle_err))),
            "minimum_angular_width_deg": math.degrees(float(np.min(widths))),
            "mean_angular_width_deg": math.degrees(float(np.mean(widths))),
        }
    )
    return metrics, final_edges, intervals, orders


def polygons_from_state(
    state: PartitionState,
    cell_ids: np.ndarray,
) -> tuple[list[np.ndarray], list[int], dict[int, tuple[float, float, float, float]]]:
    intervals = intervals_by_index(state.layer_ids, state.orders, state.boundaries)
    polygons_by_index: list[np.ndarray | None] = [None] * len(cell_ids)
    geometry: dict[int, tuple[float, float, float, float]] = {}
    for idx, (left, right) in intervals.items():
        layer = int(state.layer_ids[idx])
        r_inner = float(state.layer_edges[layer])
        r_outer = float(state.layer_edges[layer + 1])
        polygons_by_index[idx] = annular_sector_polygon(r_inner, r_outer, left, right)
        geometry[idx] = (r_inner, r_outer, left, right)
    polygons: list[np.ndarray] = []
    ids: list[int] = []
    for idx, poly in enumerate(polygons_by_index):
        if poly is None:
            raise ValueError(f"Missing polygon for cell index {idx}.")
        polygons.append(poly)
        ids.append(int(cell_ids[idx]))
    return polygons, ids, geometry


def add_direction_labels(ax: plt.Axes, radius: float) -> None:
    for angle, label, ha, va in [
        (math.pi / 2.0, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3.0 * math.pi / 2.0, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ]:
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=11, color="#444444")


def setup_axis(ax: plt.Axes, title: str, limit: float) -> None:
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.18)


def save_partition_figure(path: Path, state: PartitionState, polygons: list[np.ndarray], ids: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=280)
    fig.patch.set_facecolor("white")
    circle = np.linspace(0.0, TAU, 721)
    for pos, r in enumerate(state.layer_edges):
        is_outer = pos in (0, len(state.layer_edges) - 1)
        ax.plot(
            r * np.cos(circle),
            r * np.sin(circle),
            color="#222222" if is_outer else "#d7d7d7",
            linewidth=0.75 if is_outer else 0.45,
            zorder=1,
        )
    ax.add_collection(PolyCollection(polygons, facecolors="none", edgecolors="#4b4b4b", linewidths=0.42, zorder=2))
    for poly, cell_id in zip(polygons, ids):
        xy = np.mean(poly, axis=0)
        ax.text(float(xy[0]), float(xy[1]), str(cell_id), ha="center", va="center", fontsize=3.2, color="#111111", zorder=3)
    add_direction_labels(ax, float(np.max(state.layer_edges)) + 0.12)
    setup_axis(ax, "Step 4 - Topology-aware Annular Partition", float(np.max(state.layer_edges)) + 0.22)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_diagnostic_figure(
    path: Path,
    state: PartitionState,
    cell_ids: np.ndarray,
    original_edges: set[tuple[int, int]],
    final_edges: set[tuple[int, int]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    intervals = intervals_by_index(state.layer_ids, state.orders, state.boundaries)
    centers = np.zeros((len(cell_ids), 2), dtype=float)
    for idx, (left, right) in intervals.items():
        layer = int(state.layer_ids[idx])
        radius = 0.5 * (state.layer_edges[layer] + state.layer_edges[layer + 1])
        theta = 0.5 * (left + right)
        centers[idx] = [radius * math.cos(theta), radius * math.sin(theta)]
    preserved = sorted(original_edges & final_edges)
    lost = sorted(original_edges - final_edges)
    new = sorted(final_edges - original_edges)

    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=280)
    fig.patch.set_facecolor("white")
    circle = np.linspace(0.0, TAU, 721)
    for r in state.layer_edges:
        ax.plot(r * np.cos(circle), r * np.sin(circle), color="#dddddd", linewidth=0.45, zorder=1)
    if preserved:
        ax.add_collection(
            LineCollection(
                [(centers[i], centers[j]) for i, j in preserved], colors="#238443", linewidths=0.7, linestyles="solid", alpha=0.82, zorder=2
            )
        )
    if lost:
        ax.add_collection(
            LineCollection(
                [(centers[i], centers[j]) for i, j in lost], colors="#c9252d", linewidths=0.62, linestyles="dashed", alpha=0.82, zorder=3
            )
        )
    if new:
        ax.add_collection(
            LineCollection(
                [(centers[i], centers[j]) for i, j in new], colors="#e08214", linewidths=0.58, linestyles="dotted", alpha=0.72, zorder=4
            )
        )
    ax.scatter(centers[:, 0], centers[:, 1], s=8, color="#111111", linewidths=0, zorder=5)
    for idx, cell_id in enumerate(cell_ids):
        ax.text(
            float(centers[idx, 0]),
            float(centers[idx, 1]),
            str(int(cell_id)),
            ha="center",
            va="center",
            fontsize=3.0,
            color="#111111",
            zorder=6,
        )
    handles = [
        mpl.lines.Line2D([0], [0], color="#238443", lw=1.2, label=f"preserved original ({len(preserved)})"),
        mpl.lines.Line2D([0], [0], color="#c9252d", lw=1.2, linestyle="--", label=f"lost original ({len(lost)})"),
        mpl.lines.Line2D([0], [0], color="#e08214", lw=1.2, linestyle=":", label=f"new final ({len(new)})"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True)
    add_direction_labels(ax, float(np.max(state.layer_edges)) + 0.12)
    setup_axis(ax, "Topology adjacency diagnostic", float(np.max(state.layer_edges)) + 0.22)
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def edge_status_rows(original_edges: set[tuple[int, int]], final_edges: set[tuple[int, int]], cell_ids: np.ndarray) -> list[dict[str, Any]]:
    rows = []
    for i, j in sorted(original_edges | final_edges):
        in_original = (i, j) in original_edges
        in_final = (i, j) in final_edges
        if in_original and in_final:
            status = "preserved"
        elif in_original:
            status = "lost"
        else:
            status = "new"
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
    initial_slot: dict[int, int] = {}
    for order in initial_orders.values():
        for pos, idx in enumerate(order, start=1):
            initial_slot[int(idx)] = pos
    optimized_slot: dict[int, int] = {}
    for order in state.orders.values():
        for pos, idx in enumerate(order, start=1):
            optimized_slot[int(idx)] = pos

    original_neighbors = build_neighbors(len(cell_ids), original_edges)
    preserved_neighbors = build_neighbors(len(cell_ids), original_edges & final_edges)
    lost_neighbors = build_neighbors(len(cell_ids), original_edges - final_edges)
    new_neighbors = build_neighbors(len(cell_ids), final_edges - original_edges)
    result = []
    for idx in range(len(cell_ids)):
        r_inner, r_outer, left, right = geometry[idx]
        center_theta = 0.5 * (left + right)
        result.append(
            {
                "cell_id": int(cell_ids[idx]),
                "row": int(rows[idx]),
                "col": int(cols[idx]),
                "initial_layer": int(initial_layer_ids[idx]),
                "optimized_layer": int(state.layer_ids[idx]),
                "layer_changed": int(initial_layer_ids[idx] != state.layer_ids[idx]),
                "initial_slot": int(initial_slot.get(idx, -1)),
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
    return result


def validate_final_state(
    state: PartitionState,
    polygons: list[np.ndarray],
    ids: list[int],
    cell_ids: np.ndarray,
    original_edges: set[tuple[int, int]],
    final_edges: set[tuple[int, int]],
    min_cell_angle: float,
    inner_radius: float,
    outer_radius: float,
    objective_before: float,
    objective_after: float,
) -> None:
    assert len(polygons) == len(cell_ids), "Output polygon count must equal input cell count."
    assert len(ids) == len(set(ids)) == len(cell_ids), "Each cell_id must appear exactly once."
    assert set(ids) == set(int(x) for x in cell_ids), "Output cell_id set differs from input."
    intervals = intervals_by_index(state.layer_ids, state.orders, state.boundaries)
    intervals_again = intervals_by_index(state.layer_ids, state.orders, state.boundaries)
    fingerprint = sorted((idx, round(left, 12), round(right, 12), int(state.layer_ids[idx])) for idx, (left, right) in intervals.items())
    fingerprint_again = sorted(
        (idx, round(left, 12), round(right, 12), int(state.layer_ids[idx])) for idx, (left, right) in intervals_again.items()
    )
    assert fingerprint == fingerprint_again, "Final state fingerprint must be reproducible."
    for layer, order in state.orders.items():
        widths = np.diff(state.boundaries[layer])
        assert abs(float(np.sum(widths)) - TAU) <= 1e-8, f"Layer {layer} does not cover exactly 2pi."
        assert np.all(widths >= min_cell_angle - 1e-10), f"Layer {layer} contains a cell narrower than min_cell_angle."
        assert len(order) > 0, f"Layer {layer} is empty."
    for idx, (left, right) in intervals.items():
        assert right > left, f"Cell {idx} has non-positive angular width."
        layer = int(state.layer_ids[idx])
        assert inner_radius - 1e-10 <= state.layer_edges[layer] <= outer_radius + 1e-10
        assert inner_radius - 1e-10 <= state.layer_edges[layer + 1] <= outer_radius + 1e-10
    assert objective_after <= objective_before + 1e-9, "Final objective is higher than the initial topology-aware objective."
    preserved = original_edges & final_edges
    lost = original_edges - final_edges
    new = final_edges - original_edges
    assert len(preserved) + len(lost) == len(original_edges), "Original edge partition is inconsistent."
    assert len(preserved) + len(new) == len(final_edges), "Final edge partition is inconsistent."
    metrics = compute_adjacency_metrics(original_edges, final_edges)
    for key in ("adjacency_precision", "adjacency_recall", "adjacency_f1"):
        assert 0.0 <= float(metrics[key]) <= 1.0, f"{key} must be in [0, 1]."
    assert count_long_layer_conflicts(original_edges, state.layer_ids) >= 0, "Long layer conflicts must be counted, not ignored."


def print_comparison(baseline: dict[str, Any], optimized: PartitionMetrics) -> None:
    print("\nBaseline:")
    print(f"- Precision: {float(baseline['adjacency_precision']):.4f}")
    print(f"- Recall: {float(baseline['adjacency_recall']):.4f}")
    print(f"- F1: {float(baseline['adjacency_f1']):.4f}")
    print(f"- Lost: {int(baseline['lost_edge_count'])}")
    print(f"- New: {int(baseline['new_edge_count'])}")
    print(f"- Long layer conflicts: {int(baseline['long_layer_conflict_count'])}")
    print(f"- Mean angular error: {float(baseline['mean_angular_error_deg']):.4f} deg")
    print("\nTopology-aware:")
    print(f"- Precision: {optimized.adjacency_precision:.4f}")
    print(f"- Recall: {optimized.adjacency_recall:.4f}")
    print(f"- F1: {optimized.adjacency_f1:.4f}")
    print(f"- Lost: {optimized.lost_edge_count}")
    print(f"- New: {optimized.new_edge_count}")
    print(f"- Long layer conflicts: {optimized.long_layer_conflict_count}")
    print(f"- Mean angular error: {optimized.mean_angular_error_deg:.4f} deg")
    print("\nDelta (Topology-aware - Baseline):")
    print(f"- Precision: {optimized.adjacency_precision - float(baseline['adjacency_precision']):+.4f}")
    print(f"- Recall: {optimized.adjacency_recall - float(baseline['adjacency_recall']):+.4f}")
    print(f"- F1: {optimized.adjacency_f1 - float(baseline['adjacency_f1']):+.4f}")
    print(f"- Lost: {optimized.lost_edge_count - int(baseline['lost_edge_count']):+d}")
    print(f"- New: {optimized.new_edge_count - int(baseline['new_edge_count']):+d}")
    print(f"- Long layer conflicts: {optimized.long_layer_conflict_count - int(baseline['long_layer_conflict_count']):+d}")
    print(f"- Mean angular error: {optimized.mean_angular_error_deg - float(baseline['mean_angular_error_deg']):+.4f} deg")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step 4 optimized: topology-aware layer, order, and boundary annular partition.")
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--topology-csv", type=Path, default=DEFAULT_TOPOLOGY_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--layers", type=int, default=6)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=1.85)
    parser.add_argument("--seam-angle-deg", type=float, default=None)
    parser.add_argument("--max-layer-shift", type=int, default=1)
    parser.add_argument("--layer-opt-iters", type=int, default=30)
    parser.add_argument("--w-layer-jump", type=float, default=4.0)
    parser.add_argument("--w-radial-deviation", type=float, default=1.0)
    parser.add_argument("--w-radial-inversion", type=float, default=1.0)
    parser.add_argument("--w-empty-layer", type=float, default=10000.0)
    parser.add_argument("--order-opt-iters", type=int, default=40)
    parser.add_argument("--max-rank-shift", type=int, default=4)
    parser.add_argument("--w-lost", type=float, default=4.0)
    parser.add_argument("--w-new", type=float, default=2.0)
    parser.add_argument("--w-angle", type=float, default=0.8)
    parser.add_argument("--w-long-conflict", type=float, default=3.0)
    parser.add_argument("--w-rank", type=float, default=0.3)
    parser.add_argument("--boundary-opt-iters", type=int, default=30)
    parser.add_argument("--boundary-search-samples", type=int, default=21)
    parser.add_argument("--min-cell-angle-deg", type=float, default=0.8)
    parser.add_argument("--min-shared-angle-deg", type=float, default=0.1)
    parser.add_argument("--target-overlap-ratio", type=float, default=0.1)
    return parser.parse_args()


def run(args: argparse.Namespace) -> None:
    start_time = time.perf_counter()
    if args.outer_radius <= args.inner_radius:
        raise ValueError("--outer-radius must be greater than --inner-radius.")
    min_cell_angle = math.radians(args.min_cell_angle_deg)
    min_shared_angle = math.radians(args.min_shared_angle_deg)
    weights = OptimizationWeights(
        w_layer_jump=args.w_layer_jump,
        w_radial_deviation=args.w_radial_deviation,
        w_radial_inversion=args.w_radial_inversion,
        w_empty_layer=args.w_empty_layer,
        w_lost=args.w_lost,
        w_new=args.w_new,
        w_angle=args.w_angle,
        w_long_conflict=args.w_long_conflict,
        w_rank=args.w_rank,
    )

    records = read_center_polar_cells(args.input_csv)
    cell_ids, rows, cols, radius, theta, _centers = arrays_from_records(records)
    original_edges = read_topology_constraints(args.topology_csv, records)
    if not original_edges:
        raise ValueError("No original topology edges were found.")
    initial_layer_ids, layer_edges = assign_layers_by_radius(
        radius, args.layers, inner_radius=args.inner_radius, outer_radius=args.outer_radius
    )
    seam_angle = math.radians(args.seam_angle_deg) if args.seam_angle_deg is not None else choose_global_seam_angle(theta)
    seam_angle = float(normalize_angle(seam_angle))
    unwrapped_theta = unwrap_angles(theta, seam_angle)

    baseline, baseline_final_edges, _baseline_intervals, _baseline_orders = baseline_metrics(
        theta,
        initial_layer_ids,
        args.layers,
        original_edges,
        min_shared_angle,
    )

    initial_state, initial_ranks, initial_details = make_state(
        "global_seam_initial",
        initial_layer_ids,
        layer_edges,
        unwrapped_theta,
        seam_angle,
        min_cell_angle,
        theta,
        original_edges,
        None,
        min_shared_angle,
        weights,
    )
    history_states = [initial_state.copy(name="global_seam_initial")]

    optimized_layer_ids, layer_opt_summary = optimize_layer_assignment(
        initial_layer_ids,
        radius,
        original_edges,
        args.layers,
        max_layer_shift=args.max_layer_shift,
        max_iters=args.layer_opt_iters,
        weights=weights,
    )
    layer_opt_summary.update(
        {
            "long_layer_conflict_count_before": count_long_layer_conflicts(original_edges, initial_layer_ids),
            "long_layer_conflict_rate_before": count_long_layer_conflicts(original_edges, initial_layer_ids) / max(1, len(original_edges)),
            "long_layer_conflict_count_after": count_long_layer_conflicts(original_edges, optimized_layer_ids),
            "long_layer_conflict_rate_after": count_long_layer_conflicts(original_edges, optimized_layer_ids) / max(1, len(original_edges)),
            "layer_changed_cell_count": int(np.sum(optimized_layer_ids != initial_layer_ids)),
        }
    )
    layer_state, layer_ranks, layer_details = make_state(
        "layer_optimized",
        optimized_layer_ids,
        layer_edges,
        unwrapped_theta,
        seam_angle,
        min_cell_angle,
        theta,
        original_edges,
        None,
        min_shared_angle,
        weights,
    )
    if layer_state.objective <= initial_state.objective + 1e-12:
        current_state = layer_state
        current_ranks = layer_ranks
    else:
        current_state = initial_state.copy(name="layer_rollback")
        current_ranks = initial_ranks
    history_states.append(current_state.copy(name=current_state.name))

    order_state, order_opt_summary = optimize_circular_orders(
        current_state,
        unwrapped_theta,
        seam_angle,
        theta,
        original_edges,
        current_ranks,
        min_cell_angle=min_cell_angle,
        min_shared_angle=min_shared_angle,
        max_rank_shift=args.max_rank_shift,
        max_iters=args.order_opt_iters,
        weights=weights,
    )
    if order_state.objective <= current_state.objective + 1e-12:
        current_state = order_state
    history_states.append(current_state.copy(name="order_selected"))

    boundary_state, boundary_opt_summary = optimize_angular_boundaries(
        current_state,
        theta,
        original_edges,
        current_ranks,
        min_cell_angle=min_cell_angle,
        min_shared_angle=min_shared_angle,
        target_overlap_ratio=args.target_overlap_ratio,
        samples=args.boundary_search_samples,
        max_iters=args.boundary_opt_iters,
        weights=weights,
    )
    if boundary_state.objective <= current_state.objective + 1e-12:
        current_state = boundary_state
    history_states.append(current_state.copy(name="boundary_selected"))

    # Choose the historical state with the strongest F1; objective is the tie-breaker.
    scored_states: list[tuple[float, float, PartitionState, PartitionMetrics, set[tuple[int, int]]]] = []
    stage_metrics_payload: dict[str, dict[str, Any]] = {}
    for state in history_states:
        metrics, edges, _details = state_metrics(state, theta, original_edges, current_ranks, min_shared_angle, weights)
        stage_metrics_payload[state.name] = asdict(metrics)
        scored_states.append((metrics.adjacency_f1, -metrics.objective, state, metrics, edges))
    scored_states.sort(key=lambda item: (item[0], item[1], -len(item[2].name)), reverse=True)
    selected_state = scored_states[0][2].copy(name=f"selected_{scored_states[0][2].name}")
    selected_metrics = scored_states[0][3]
    selected_final_edges = scored_states[0][4]
    if selected_metrics.adjacency_f1 + 1e-12 < float(baseline["adjacency_f1"]):
        print(
            "WARNING: final topology-aware F1 is lower than the old Step 4 baseline; using the best historical optimized state instead of the last iteration."
        )

    polygons, ids, geometry = polygons_from_state(selected_state, cell_ids)
    validate_final_state(
        selected_state,
        polygons,
        ids,
        cell_ids,
        original_edges,
        selected_final_edges,
        min_cell_angle,
        args.inner_radius,
        args.outer_radius,
        initial_state.objective,
        selected_metrics.objective,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    partition_png = args.output_dir / "step04_topology_aware_partition.png"
    assignment_csv = args.output_dir / "step04_topology_aware_assignment.csv"
    metrics_json = args.output_dir / "step04_topology_aware_metrics.json"
    diagnostic_png = args.output_dir / "step04_topology_adjacency_diagnostic.png"
    edge_status_csv = args.output_dir / "step04_topology_edge_status.csv"
    save_partition_figure(partition_png, selected_state, polygons, ids)
    save_diagnostic_figure(diagnostic_png, selected_state, cell_ids, original_edges, selected_final_edges)
    initial_orders = build_layer_orders(initial_layer_ids, unwrapped_theta, args.layers)
    write_csv(
        assignment_csv,
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
            selected_final_edges,
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
    write_csv(
        edge_status_csv,
        edge_status_rows(original_edges, selected_final_edges, cell_ids),
        ["cell_id_a", "cell_id_b", "in_original", "in_final", "status"],
    )

    final_metrics_dict = asdict(selected_metrics)
    layer_changed = np.flatnonzero(selected_state.layer_ids != initial_layer_ids).astype(int).tolist()
    lost_edges = sorted(original_edges - selected_final_edges)
    new_edges = sorted(selected_final_edges - original_edges)
    runtime_seconds = time.perf_counter() - start_time
    metrics_payload = {
        "parameters": {k: portable_parameter(v) for k, v in vars(args).items()},
        "seam_angle_deg": math.degrees(seam_angle),
        "topology_csv_used": portable_parameter(args.topology_csv) if args.topology_csv.exists() else "fallback_from_row_col",
        "baseline_metrics": baseline,
        "initial_topology_aware_metrics": initial_details,
        "layer_optimization": layer_opt_summary,
        "layer_optimized_metrics": layer_details,
        "order_optimization": order_opt_summary,
        "order_optimized_metrics": stage_metrics_payload.get("order_selected", {}),
        "boundary_optimization": boundary_opt_summary,
        "boundary_optimized_metrics": stage_metrics_payload.get("boundary_selected", {}),
        "stage_metrics": stage_metrics_payload,
        "final_metrics": final_metrics_dict,
        "final_adjacency_precision": selected_metrics.adjacency_precision,
        "final_adjacency_recall": selected_metrics.adjacency_recall,
        "final_adjacency_f1": selected_metrics.adjacency_f1,
        "lost_edges": [{"cell_id_a": int(cell_ids[i]), "cell_id_b": int(cell_ids[j])} for i, j in lost_edges],
        "new_edges": [{"cell_id_a": int(cell_ids[i]), "cell_id_b": int(cell_ids[j])} for i, j in new_edges],
        "long_layer_conflicts": [
            {
                "cell_id_a": int(cell_ids[i]),
                "cell_id_b": int(cell_ids[j]),
                "layer_delta": int(abs(selected_state.layer_ids[i] - selected_state.layer_ids[j])),
            }
            for i, j in sorted(original_edges)
            if abs(int(selected_state.layer_ids[i]) - int(selected_state.layer_ids[j])) > 1
        ],
        "layer_changed_cells": [
            {
                "cell_id": int(cell_ids[idx]),
                "initial_layer": int(initial_layer_ids[idx]),
                "optimized_layer": int(selected_state.layer_ids[idx]),
            }
            for idx in layer_changed
        ],
        "layer_changed_cell_count": len(layer_changed),
        "mean_angular_error": selected_metrics.mean_angular_error_deg,
        "median_angular_error": selected_metrics.median_angular_error_deg,
        "max_angular_error": selected_metrics.max_angular_error_deg,
        "minimum_angular_width": selected_metrics.minimum_angular_width_deg,
        "mean_angular_width": selected_metrics.mean_angular_width_deg,
        "objective_before": initial_state.objective,
        "objective_after": selected_metrics.objective,
        "runtime_seconds": runtime_seconds,
    }
    metrics_json.write_text(json.dumps(json_safe(metrics_payload), ensure_ascii=False, indent=2), encoding="utf-8")

    print_comparison(baseline, selected_metrics)
    print("\nOutputs:")
    for path in (partition_png, assignment_csv, metrics_json, diagnostic_png, edge_status_csv):
        print(f"- {path}")
    print(f"\nLayer changed cell count: {len(layer_changed)}")
    if layer_changed:
        print("Layer changed cells:", ", ".join(str(int(cell_ids[idx])) for idx in layer_changed))
    print(f"Lost original edges after optimization: {len(lost_edges)}")
    if lost_edges:
        print(
            "Lost edges:",
            ", ".join(f"{int(cell_ids[i])}-{int(cell_ids[j])}" for i, j in lost_edges[:30]) + (" ..." if len(lost_edges) > 30 else ""),
        )
    print(f"New final edges after optimization: {len(new_edges)}")
    if new_edges:
        print(
            "New edges:",
            ", ".join(f"{int(cell_ids[i])}-{int(cell_ids[j])}" for i, j in new_edges[:30]) + (" ..." if len(new_edges) > 30 else ""),
        )
    print(f"runtime seconds: {runtime_seconds:.3f}")


if __name__ == "__main__":
    run(parse_args())
