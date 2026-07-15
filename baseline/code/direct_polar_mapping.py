#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import math
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
STEP01_CSV = PROJECT_ROOT / "data" / "input" / "beijing_grid_cells.csv"
BOUNDARY_RAYS_CSV = PROJECT_ROOT / "data" / "input" / "beijing_boundary_rays.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "baselines" / "polar"
INNER_RADIUS = 1.0
OUTER_RADIUS = 1.85
RING_AREA = math.pi * (OUTER_RADIUS**2 - INNER_RADIUS**2)


@dataclass(frozen=True)
class SourceCell:
    cell_id: int
    row: int
    col: int
    center_lon: float
    center_lat: float


@dataclass
class BaselineResult:
    cells: list[SourceCell]
    source_polygons: list[np.ndarray]
    annulus_polygons: list[np.ndarray]
    centers: np.ndarray
    corner_theta: list[np.ndarray]
    corner_rho: list[np.ndarray]


def read_step01_cells(path: Path) -> list[SourceCell]:
    with path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    cells = [
        SourceCell(
            cell_id=int(row["cell_id"]),
            row=int(row["row"]),
            col=int(row["col"]),
            center_lon=float(row["center_lon"]),
            center_lat=float(row["center_lat"]),
        )
        for row in rows
    ]
    if not cells:
        raise ValueError(f"No cells found in {path}")
    return cells


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def project_lonlat(points: np.ndarray, center_lon: float, center_lat: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return np.column_stack(
        [
            (points[:, 0] - center_lon) * math.cos(math.radians(center_lat)),
            points[:, 1] - center_lat,
        ]
    )


def read_boundary_rays(path: Path = BOUNDARY_RAYS_CSV) -> dict[tuple[float, float], float]:
    """Read precomputed boundary distances for every source-grid corner.

    The compact ray table replaces the original machine-local dependency on a
    full administrative-boundary GeoJSON while preserving the baseline output.
    """
    lookup: dict[tuple[float, float], float] = {}
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (round(float(row["source_lon"]), 10), round(float(row["source_lat"]), 10))
            lookup[key] = float(row["boundary_radius"])
    return lookup


def infer_grid_steps(cells: list[SourceCell]) -> tuple[float, float]:
    lons = np.asarray(sorted({round(cell.center_lon, 6) for cell in cells}), dtype=float)
    lats = np.asarray(sorted({round(cell.center_lat, 6) for cell in cells}), dtype=float)
    lon_diffs = np.diff(lons)
    lat_diffs = np.diff(lats)
    lon_step = float(np.nanmedian(lon_diffs[lon_diffs > 1e-8]))
    lat_step = float(np.nanmedian(lat_diffs[lat_diffs > 1e-8]))
    return lon_step, lat_step


def source_cell_polygon(cell: SourceCell, lon_step: float, lat_step: float) -> np.ndarray:
    half_lon = lon_step / 2.0
    half_lat = lat_step / 2.0
    # Four original grid corners: lower-left, lower-right, upper-right, upper-left.
    return np.asarray(
        [
            [cell.center_lon - half_lon, cell.center_lat - half_lat],
            [cell.center_lon + half_lon, cell.center_lat - half_lat],
            [cell.center_lon + half_lon, cell.center_lat + half_lat],
            [cell.center_lon - half_lon, cell.center_lat + half_lat],
        ],
        dtype=float,
    )


def choose_origin_from_center_cells(cells: list[SourceCell]) -> tuple[float, float]:
    center_cells = [cell for cell in cells if cell.cell_id in {96, 97, 110, 111}]
    if len(center_cells) == 4:
        return (
            float(np.mean([cell.center_lon for cell in center_cells])),
            float(np.mean([cell.center_lat for cell in center_cells])),
        )
    return (
        float(np.mean([cell.center_lon for cell in cells])),
        float(np.mean([cell.center_lat for cell in cells])),
    )


def map_polygon_to_annulus(
    polygon_lonlat: np.ndarray,
    *,
    center_lon: float,
    center_lat: float,
    boundary_rays: dict[tuple[float, float], float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    projected = project_lonlat(polygon_lonlat, center_lon, center_lat)
    theta = np.mod(np.arctan2(projected[:, 1], projected[:, 0]), 2.0 * math.pi)
    distance = np.hypot(projected[:, 0], projected[:, 1])
    keys = [(round(float(lon), 10), round(float(lat), 10)) for lon, lat in polygon_lonlat]
    try:
        boundary_radius = np.asarray([boundary_rays[key] for key in keys], dtype=float)
    except KeyError as exc:
        raise KeyError(f"Missing precomputed boundary ray for source-grid corner {exc.args[0]}") from exc
    rho = distance / np.maximum(boundary_radius, 1e-12)
    radius = INNER_RADIUS + rho * (OUTER_RADIUS - INNER_RADIUS)
    mapped = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])
    return mapped, theta, rho


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
    return np.asarray(
        [
            float(np.sum((x + np.roll(x, -1)) * cross) / (6.0 * area)),
            float(np.sum((y + np.roll(y, -1)) * cross) / (6.0 * area)),
        ],
        dtype=float,
    )


def edge_lengths(poly: np.ndarray) -> np.ndarray:
    return np.linalg.norm(np.roll(poly, -1, axis=0) - poly, axis=1)


def topology_edges(cells: list[SourceCell]) -> list[tuple[int, int]]:
    index_by_grid = {(cell.row, cell.col): i for i, cell in enumerate(cells)}
    edges: list[tuple[int, int]] = []
    for i, cell in enumerate(cells):
        for key in ((cell.row, cell.col + 1), (cell.row + 1, cell.col)):
            j = index_by_grid.get(key)
            if j is not None:
                edges.append((i, j))
    return edges


def connected_components(node_count: int, edges: list[tuple[int, int]]) -> int:
    graph = [[] for _ in range(node_count)]
    for i, j in edges:
        graph[i].append(j)
        graph[j].append(i)
    seen = np.zeros(node_count, dtype=bool)
    components = 0
    for start in range(node_count):
        if seen[start]:
            continue
        components += 1
        queue: deque[int] = deque([start])
        seen[start] = True
        while queue:
            cur = queue.popleft()
            for nxt in graph[cur]:
                if not seen[nxt]:
                    seen[nxt] = True
                    queue.append(nxt)
    return components


def ensure_ccw(poly: np.ndarray) -> np.ndarray:
    return poly if polygon_area(poly) >= 0 else poly[::-1].copy()


def clip_polygon(subject: np.ndarray, clipper: np.ndarray) -> np.ndarray:
    subject = ensure_ccw(subject)
    clipper = ensure_ccw(clipper)

    def cross2(u: np.ndarray, v: np.ndarray) -> float:
        return float(u[0] * v[1] - u[1] * v[0])

    def inside(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> bool:
        return cross2(b - a, p - a) >= -1e-10

    def intersect(p1: np.ndarray, p2: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        d1 = p2 - p1
        d2 = b - a
        denom = cross2(d1, d2)
        if abs(denom) < 1e-12:
            return p2
        t = cross2(a - p1, d2) / denom
        return p1 + t * d1

    output = subject.copy()
    for a, b in zip(clipper, np.roll(clipper, -1, axis=0)):
        if len(output) == 0:
            break
        current = output
        result: list[np.ndarray] = []
        s = current[-1]
        for e in current:
            if inside(e, a, b):
                if not inside(s, a, b):
                    result.append(intersect(s, e, a, b))
                result.append(e)
            elif inside(s, a, b):
                result.append(intersect(s, e, a, b))
            s = e
        output = np.asarray(result, dtype=float) if result else np.empty((0, 2), dtype=float)
    return output


def polygon_overlap_area(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or len(b) < 3:
        return 0.0
    amin, amax = np.min(a, axis=0), np.max(a, axis=0)
    bmin, bmax = np.min(b, axis=0), np.max(b, axis=0)
    if np.any(amax < bmin) or np.any(bmax < amin):
        return 0.0
    inter = clip_polygon(a, b)
    if len(inter) < 3:
        return 0.0
    return abs(polygon_area(inter))


def non_adjacent_overlap(polygons: list[np.ndarray], edges: list[tuple[int, int]]) -> tuple[int, float]:
    adjacent = {tuple(sorted(edge)) for edge in edges}
    pairs = 0
    total_area = 0.0
    for i in range(len(polygons)):
        for j in range(i + 1, len(polygons)):
            if (i, j) in adjacent:
                continue
            area = polygon_overlap_area(polygons[i], polygons[j])
            if area > 1e-7:
                pairs += 1
                total_area += area
    return pairs, total_area


def classify_edges_by_layers(centers: np.ndarray, edges: list[tuple[int, int]], layer_count: int = 6) -> dict[str, int | float]:
    radius = np.hypot(centers[:, 0], centers[:, 1])
    quantiles = np.quantile(radius, np.linspace(0.0, 1.0, layer_count + 1))
    quantiles[0] -= 1e-9
    quantiles[-1] += 1e-9
    layers = np.clip(np.searchsorted(quantiles, radius, side="right") - 1, 0, layer_count - 1)
    same = adjacent = long = 0
    for i, j in edges:
        delta = abs(int(layers[i]) - int(layers[j]))
        if delta == 0:
            same += 1
        elif delta == 1:
            adjacent += 1
        else:
            long += 1
    total = max(len(edges), 1)
    return {
        "same_layer_edges": same,
        "adjacent_layer_edges": adjacent,
        "long_diagonal_edges": long,
        "long_diagonal_ratio": long / total,
    }


def build_result() -> BaselineResult:
    cells = read_step01_cells(STEP01_CSV)
    center_lon, center_lat = choose_origin_from_center_cells(cells)
    boundary_rays = read_boundary_rays()
    lon_step, lat_step = infer_grid_steps(cells)

    source_polygons: list[np.ndarray] = []
    annulus_polygons: list[np.ndarray] = []
    corner_theta: list[np.ndarray] = []
    corner_rho: list[np.ndarray] = []
    for cell in cells:
        source_poly = source_cell_polygon(cell, lon_step, lat_step)
        mapped, theta, rho = map_polygon_to_annulus(
            source_poly,
            center_lon=center_lon,
            center_lat=center_lat,
            boundary_rays=boundary_rays,
        )
        source_polygons.append(source_poly)
        annulus_polygons.append(mapped)
        corner_theta.append(theta)
        corner_rho.append(rho)

    centers = np.asarray([polygon_centroid(poly) for poly in annulus_polygons], dtype=float)
    return BaselineResult(cells, source_polygons, annulus_polygons, centers, corner_theta, corner_rho)


def direction_labels(ax, radius: float) -> None:
    for angle, label, ha, va in [
        (math.pi / 2.0, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3.0 * math.pi / 2.0, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ]:
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=10, color="#444444")


def save_plot(result: BaselineResult, output_path: Path, *, draw_edges: bool = False) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.6, 9.6), dpi=260)
    fig.patch.set_facecolor("white")
    theta = np.linspace(0.0, 2.0 * math.pi, 721)
    ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.add_collection(PolyCollection(result.annulus_polygons, facecolors="none", edgecolors="#5d5d5d", linewidths=0.34, zorder=2))
    if draw_edges:
        segments = [(result.centers[i], result.centers[j]) for i, j in topology_edges(result.cells)]
        ax.add_collection(LineCollection(segments, colors="#c9483a", linewidths=0.48, alpha=0.82, zorder=3))
    ax.scatter(result.centers[:, 0], result.centers[:, 1], s=4.2, color="#111111", linewidths=0, zorder=4)
    for xy, cell in zip(result.centers, result.cells):
        ax.text(float(xy[0]), float(xy[1]), str(cell.cell_id), ha="center", va="center", fontsize=2.7, color="#111111", zorder=5)
    direction_labels(ax, OUTER_RADIUS + 0.13)
    ax.set_title("Baseline 1 - Four-Corner Direct Polar Mapping", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    coords = np.vstack(result.annulus_polygons)
    limit = max(float(np.abs(coords).max()), OUTER_RADIUS) + 0.22
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.20)
    fig.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def metrics(result: BaselineResult) -> dict[str, object]:
    edges = topology_edges(result.cells)
    areas_signed = np.asarray([polygon_area(poly) for poly in result.annulus_polygons], dtype=float)
    areas = np.abs(areas_signed)
    perimeters = np.asarray([float(np.sum(edge_lengths(poly))) for poly in result.annulus_polygons], dtype=float)
    aspect = []
    for poly in result.annulus_polygons:
        lengths = edge_lengths(poly)
        positive = lengths[lengths > 1e-12]
        aspect.append(float(np.max(positive) / np.min(positive)) if len(positive) else np.nan)
    aspect = np.asarray(aspect, dtype=float)
    vertices = np.vstack(result.annulus_polygons)
    vertex_radius = np.hypot(vertices[:, 0], vertices[:, 1])
    center_radius = np.hypot(result.centers[:, 0], result.centers[:, 1])
    edge_lengths_mapped = np.asarray([float(np.linalg.norm(result.centers[i] - result.centers[j])) for i, j in edges], dtype=float)
    overlap_pairs, overlap_area = non_adjacent_overlap(result.annulus_polygons, edges)
    area_sum = float(np.sum(areas))
    return {
        "experiment": "baseline_polar_mapping",
        "title": "Baseline 1 - four-corner direct polar mapping",
        "kind": "four_corner_direct_polar",
        "cell_count": len(result.cells),
        "topology_edges": len(edges),
        "connected_components": connected_components(len(result.cells), edges),
        "zero_area_cells": int(np.sum(areas <= 1e-10)),
        "clockwise_orientation_cells": int(np.sum(areas_signed < -1e-10)),
        "area_sum": area_sum,
        "area_mean": float(np.mean(areas)),
        "area_cv": float(np.std(areas) / max(np.mean(areas), 1e-12)),
        "area_min": float(np.min(areas)),
        "area_max": float(np.max(areas)),
        "perimeter_mean": float(np.mean(perimeters)),
        "aspect_ratio_median": float(np.nanmedian(aspect)),
        "aspect_ratio_max": float(np.nanmax(aspect)),
        "vertices_outside_outer": int(np.sum(vertex_radius > OUTER_RADIUS + 1e-8)),
        "vertices_inside_inner": int(np.sum(vertex_radius < INNER_RADIUS - 1e-8)),
        "centers_outside_outer": int(np.sum(center_radius > OUTER_RADIUS + 1e-8)),
        "centers_inside_inner": int(np.sum(center_radius < INNER_RADIUS - 1e-8)),
        "non_adjacent_overlap_pairs": overlap_pairs,
        "non_adjacent_overlap_area": overlap_area,
        "ring_gap_area_proxy": max(RING_AREA - area_sum + overlap_area, 0.0),
        "ring_area_coverage_proxy": min(area_sum / RING_AREA, 999.0),
        "topology_edge_length_mean": float(np.mean(edge_lengths_mapped)) if len(edge_lengths_mapped) else 0.0,
        "topology_edge_length_max": float(np.max(edge_lengths_mapped)) if len(edge_lengths_mapped) else 0.0,
        "topology_edge_length_cv": float(np.std(edge_lengths_mapped) / max(np.mean(edge_lengths_mapped), 1e-12))
        if len(edge_lengths_mapped)
        else 0.0,
        **classify_edges_by_layers(result.centers, edges),
    }


def save_data(result: BaselineResult) -> None:
    corner_rows: list[dict[str, object]] = []
    for cell, source_poly, mapped_poly, theta, rho in zip(
        result.cells,
        result.source_polygons,
        result.annulus_polygons,
        result.corner_theta,
        result.corner_rho,
    ):
        for corner_index, (src, xy, th, rr) in enumerate(zip(source_poly, mapped_poly, theta, rho), start=1):
            corner_rows.append(
                {
                    "cell_id": cell.cell_id,
                    "row": cell.row,
                    "col": cell.col,
                    "corner_index": corner_index,
                    "source_lon": f"{float(src[0]):.10f}",
                    "source_lat": f"{float(src[1]):.10f}",
                    "theta_rad": f"{float(th):.10f}",
                    "theta_deg": f"{math.degrees(float(th)):.6f}",
                    "rho_unclipped": f"{float(rr):.10f}",
                    "R": f"{float(INNER_RADIUS + rr * (OUTER_RADIUS - INNER_RADIUS)):.10f}",
                    "annulus_x": f"{float(xy[0]):.10f}",
                    "annulus_y": f"{float(xy[1]):.10f}",
                }
            )
    write_csv(
        OUTPUT_DIR / "baseline_polar_mapping_corners.csv",
        corner_rows,
        [
            "cell_id",
            "row",
            "col",
            "corner_index",
            "source_lon",
            "source_lat",
            "theta_rad",
            "theta_deg",
            "rho_unclipped",
            "R",
            "annulus_x",
            "annulus_y",
        ],
    )

    cell_rows: list[dict[str, object]] = []
    for cell, center, poly in zip(result.cells, result.centers, result.annulus_polygons):
        area_signed = polygon_area(poly)
        cell_rows.append(
            {
                "cell_id": cell.cell_id,
                "row": cell.row,
                "col": cell.col,
                "center_x": f"{float(center[0]):.10f}",
                "center_y": f"{float(center[1]):.10f}",
                "vertex_count": len(poly),
                "signed_area": f"{area_signed:.10f}",
                "area": f"{abs(area_signed):.10f}",
            }
        )
    write_csv(
        OUTPUT_DIR / "baseline_polar_mapping_cells.csv",
        cell_rows,
        ["cell_id", "row", "col", "center_x", "center_y", "vertex_count", "signed_area", "area"],
    )


def write_readme(metric_row: dict[str, object]) -> None:
    text = f"""# Baseline 1 - Polar Mapping

这个 baseline 将 `data/input/beijing_grid_cells.csv` 中的 176 个北京网格 cell 作为输入，并使用 `data/input/beijing_boundary_rays.csv` 中的预计算边界射线距离。

做法：

1. 用每个 cell 的中心经纬度和网格步长重建原始四个角点。
2. 以 cell 96、97、110、111 的共享中心顶点作为极坐标原点。
3. 每个角点单独计算 `theta` 和相对北京边界射线距离的 `rho_unclipped`。
4. 直接映射到圆环半径 `R = {INNER_RADIUS} + rho_unclipped * ({OUTER_RADIUS} - {INNER_RADIUS})`。

该方法不做拓扑修正、不做层列重排、不裁剪越界角点，因此用于观察“直接四角点极坐标映射”的原始几何问题。

## 输出文件

- `baseline_polar_mapping_result.png`: 四角点直接映射后的圆环网格。
- `baseline_polar_mapping_topology_overlay.png`: 叠加原始上下左右邻接边的结果图。
- `baseline_polar_mapping_corners.csv`: 每个 cell 四个角点的原始经纬度、theta、rho、圆环坐标。
- `baseline_polar_mapping_cells.csv`: 每个 cell 的圆环中心、面积和顶点数。
- `baseline_polar_mapping_metrics.csv`: 单独指标表。
- `code/direct_polar_mapping.py`: 复现脚本。

## 当前指标摘要

- cell_count: {metric_row["cell_count"]}
- topology_edges: {metric_row["topology_edges"]}
- vertices_outside_outer: {metric_row["vertices_outside_outer"]}
- vertices_inside_inner: {metric_row["vertices_inside_inner"]}
- non_adjacent_overlap_pairs: {metric_row["non_adjacent_overlap_pairs"]}
- non_adjacent_overlap_area: {float(metric_row["non_adjacent_overlap_area"]):.10f}
- area_cv: {float(metric_row["area_cv"]):.6f}
- ring_area_coverage_proxy: {float(metric_row["ring_area_coverage_proxy"]):.6f}
"""
    (OUTPUT_DIR / "README.md").write_text(text, encoding="utf-8")


def run() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = build_result()
    save_plot(result, OUTPUT_DIR / "baseline_polar_mapping_result.png")
    save_plot(result, OUTPUT_DIR / "baseline_polar_mapping_topology_overlay.png", draw_edges=True)
    save_data(result)
    metric_row = metrics(result)
    write_csv(OUTPUT_DIR / "baseline_polar_mapping_metrics.csv", [metric_row], list(metric_row.keys()))
    write_readme(metric_row)
    print(f"cells: {len(result.cells)}")
    print(f"wrote: {OUTPUT_DIR}")


if __name__ == "__main__":
    run()
