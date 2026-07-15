#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import math
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
STEP01_CSV = PROJECT_ROOT / "data" / "input" / "beijing_grid_cells.csv"
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
class AnnulusMesh:
    cells: list[SourceCell]
    polygons: list[list[tuple]]
    source_node_coords: dict[tuple, np.ndarray]
    center_lon: float
    center_lat: float
    inner_nodes: list[tuple]
    outer_nodes: list[tuple]
    boundary_nodes: set[tuple]


@dataclass
class BaselineOutput:
    key: str
    title: str
    kind: str
    cells: list[SourceCell]
    polygons: list[np.ndarray]
    centers: np.ndarray
    node_coords: dict[tuple, np.ndarray]
    mesh: AnnulusMesh


def read_step01_cells(path: Path = STEP01_CSV) -> list[SourceCell]:
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
        writer.writerows(rows)


def infer_grid_steps(cells: list[SourceCell]) -> tuple[float, float]:
    lons = np.asarray(sorted({round(cell.center_lon, 6) for cell in cells}), dtype=float)
    lats = np.asarray(sorted({round(cell.center_lat, 6) for cell in cells}), dtype=float)
    lon_diffs = np.diff(lons)
    lat_diffs = np.diff(lats)
    return float(np.nanmedian(lon_diffs[lon_diffs > 1e-8])), float(np.nanmedian(lat_diffs[lat_diffs > 1e-8]))


def center_vertex_from_cells(cells: list[SourceCell]) -> tuple[int, int]:
    center_cells = [cell for cell in cells if cell.cell_id in {96, 97, 110, 111}]
    if len(center_cells) == 4:
        shared = None
        for cell in center_cells:
            vertices = {
                (cell.row, cell.col),
                (cell.row, cell.col + 1),
                (cell.row + 1, cell.col + 1),
                (cell.row + 1, cell.col),
            }
            shared = vertices if shared is None else shared & vertices
        if shared:
            return sorted(shared)[0]
    return (
        int(round(np.mean([cell.row for cell in cells]))),
        int(round(np.mean([cell.col for cell in cells]))),
    )


def project_lonlat(points: np.ndarray, center_lon: float, center_lat: float) -> np.ndarray:
    points = np.atleast_2d(np.asarray(points, dtype=float))
    return np.column_stack(
        [
            (points[:, 0] - center_lon) * math.cos(math.radians(center_lat)),
            points[:, 1] - center_lat,
        ]
    )


def build_punctured_mesh(cells: list[SourceCell] | None = None, *, hole_fraction: float = 0.20) -> AnnulusMesh:
    cells = cells or read_step01_cells()
    lon_step, lat_step = infer_grid_steps(cells)
    center_vertex = center_vertex_from_cells(cells)

    vertex_coords_lonlat: dict[tuple[int, int], np.ndarray] = {}
    raw_polygons: list[list[tuple]] = []
    for cell in cells:
        vertices = [
            (cell.row, cell.col),
            (cell.row, cell.col + 1),
            (cell.row + 1, cell.col + 1),
            (cell.row + 1, cell.col),
        ]
        coords = [
            np.asarray([cell.center_lon - lon_step / 2.0, cell.center_lat - lat_step / 2.0], dtype=float),
            np.asarray([cell.center_lon + lon_step / 2.0, cell.center_lat - lat_step / 2.0], dtype=float),
            np.asarray([cell.center_lon + lon_step / 2.0, cell.center_lat + lat_step / 2.0], dtype=float),
            np.asarray([cell.center_lon - lon_step / 2.0, cell.center_lat + lat_step / 2.0], dtype=float),
        ]
        for vertex, coord in zip(vertices, coords):
            vertex_coords_lonlat.setdefault(vertex, coord)
        raw_polygons.append([("v", r, c) for r, c in vertices])

    center_lon, center_lat = vertex_coords_lonlat[center_vertex]
    hole_radius = hole_fraction * min(abs(lon_step) * math.cos(math.radians(center_lat)), abs(lat_step))
    hole_lonlat = {
        ("h", "E"): np.asarray([center_lon + hole_radius / math.cos(math.radians(center_lat)), center_lat], dtype=float),
        ("h", "N"): np.asarray([center_lon, center_lat + hole_radius], dtype=float),
        ("h", "W"): np.asarray([center_lon - hole_radius / math.cos(math.radians(center_lat)), center_lat], dtype=float),
        ("h", "S"): np.asarray([center_lon, center_lat - hole_radius], dtype=float),
    }

    polygons: list[list[tuple]] = []
    center_node = ("v", center_vertex[0], center_vertex[1])
    replacement_by_position = {
        0: [("h", "E"), ("h", "N")],
        1: [("h", "S"), ("h", "E")],
        2: [("h", "W"), ("h", "S")],
        3: [("h", "N"), ("h", "W")],
    }
    for polygon in raw_polygons:
        if center_node not in polygon:
            polygons.append(polygon)
            continue
        pos = polygon.index(center_node)
        new_polygon: list[tuple] = []
        for idx, node in enumerate(polygon):
            if idx == pos:
                new_polygon.extend(replacement_by_position[pos])
            else:
                new_polygon.append(node)
        polygons.append(new_polygon)

    source_node_coords: dict[tuple, np.ndarray] = {}
    for node in {node for polygon in polygons for node in polygon}:
        if node[0] == "h":
            source_node_coords[node] = project_lonlat(hole_lonlat[node], float(center_lon), float(center_lat))[0]
        else:
            source_node_coords[node] = project_lonlat(
                vertex_coords_lonlat[(int(node[1]), int(node[2]))], float(center_lon), float(center_lat)
            )[0]

    edge_count: dict[tuple[tuple, tuple], int] = defaultdict(int)
    for polygon in polygons:
        for a, b in zip(polygon, polygon[1:] + polygon[:1]):
            edge_count[tuple(sorted((a, b)))] += 1
    inner_nodes = [("h", "E"), ("h", "N"), ("h", "W"), ("h", "S")]
    inner_set = set(inner_nodes)
    outer_nodes = sorted(
        {node for edge, count in edge_count.items() for node in edge if count == 1 and node not in inner_set},
        key=lambda node: math.atan2(float(source_node_coords[node][1]), float(source_node_coords[node][0])),
    )
    boundary_nodes = set(outer_nodes) | inner_set
    return AnnulusMesh(cells, polygons, source_node_coords, float(center_lon), float(center_lat), inner_nodes, outer_nodes, boundary_nodes)


def boundary_positions(mesh: AnnulusMesh) -> dict[tuple, np.ndarray]:
    fixed: dict[tuple, np.ndarray] = {}
    for node in mesh.outer_nodes:
        theta = math.atan2(float(mesh.source_node_coords[node][1]), float(mesh.source_node_coords[node][0]))
        fixed[node] = np.asarray([OUTER_RADIUS * math.cos(theta), OUTER_RADIUS * math.sin(theta)], dtype=float)
    inner_angles = {("h", "E"): 0.0, ("h", "N"): math.pi / 2.0, ("h", "W"): math.pi, ("h", "S"): 3.0 * math.pi / 2.0}
    for node, theta in inner_angles.items():
        fixed[node] = np.asarray([INNER_RADIUS * math.cos(theta), INNER_RADIUS * math.sin(theta)], dtype=float)
    return fixed


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


def mesh_edges(polygons: list[list[tuple]]) -> set[tuple[tuple, tuple]]:
    edges: set[tuple[tuple, tuple]] = set()
    for polygon in polygons:
        for a, b in zip(polygon, polygon[1:] + polygon[:1]):
            edges.add(tuple(sorted((a, b))))
    return edges


def uniform_weights(mesh: AnnulusMesh) -> dict[tuple, dict[tuple, float]]:
    weights: dict[tuple, dict[tuple, float]] = defaultdict(dict)
    for a, b in mesh_edges(mesh.polygons):
        weights[a][b] = weights[a].get(b, 0.0) + 1.0
        weights[b][a] = weights[b].get(a, 0.0) + 1.0
    return weights


def triangulate_polygon(polygon: list[tuple]) -> list[tuple[tuple, tuple, tuple]]:
    return [(polygon[0], polygon[i], polygon[i + 1]) for i in range(1, len(polygon) - 1)]


def cotangent(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u = b - a
    v = c - a
    denom = float(u[0] * v[1] - u[1] * v[0])
    if abs(denom) < 1e-14:
        return 0.0
    return float(np.dot(u, v) / abs(denom))


def cotangent_weights(mesh: AnnulusMesh) -> dict[tuple, dict[tuple, float]]:
    weights: dict[tuple, dict[tuple, float]] = defaultdict(dict)

    def add(a: tuple, b: tuple, value: float) -> None:
        value = max(1e-8, 0.5 * value)
        weights[a][b] = weights[a].get(b, 0.0) + value
        weights[b][a] = weights[b].get(a, 0.0) + value

    for polygon in mesh.polygons:
        for a, b, c in triangulate_polygon(polygon):
            pa = mesh.source_node_coords[a]
            pb = mesh.source_node_coords[b]
            pc = mesh.source_node_coords[c]
            add(b, c, cotangent(pa, pb, pc))
            add(c, a, cotangent(pb, pc, pa))
            add(a, b, cotangent(pc, pa, pb))
    return weights


def solve_harmonic(mesh: AnnulusMesh, weights: dict[tuple, dict[tuple, float]]) -> dict[tuple, np.ndarray]:
    fixed = boundary_positions(mesh)
    nodes = sorted(mesh.source_node_coords)
    free_nodes = [node for node in nodes if node not in fixed]
    index = {node: i for i, node in enumerate(free_nodes)}
    matrix = lil_matrix((len(free_nodes), len(free_nodes)), dtype=float)
    rhs = np.zeros((len(free_nodes), 2), dtype=float)
    for node in free_nodes:
        row = index[node]
        total = 0.0
        for other, weight in weights[node].items():
            total += weight
            if other in fixed:
                rhs[row] += weight * fixed[other]
            else:
                matrix[row, index[other]] -= weight
        matrix[row, row] += max(total, 1e-12)
    solved = spsolve(matrix.tocsr(), rhs) if free_nodes else np.empty((0, 2), dtype=float)
    coords = dict(fixed)
    for node, xy in zip(free_nodes, np.asarray(solved, dtype=float)):
        coords[node] = xy
    return coords


def polygons_from_node_coords(mesh: AnnulusMesh, node_coords: dict[tuple, np.ndarray]) -> list[np.ndarray]:
    return [np.asarray([node_coords[node] for node in polygon], dtype=float) for polygon in mesh.polygons]


def output_from_node_coords(key: str, title: str, kind: str, mesh: AnnulusMesh, node_coords: dict[tuple, np.ndarray]) -> BaselineOutput:
    polygons = polygons_from_node_coords(mesh, node_coords)
    centers = np.asarray([polygon_centroid(poly) for poly in polygons], dtype=float)
    return BaselineOutput(key, title, kind, mesh.cells, polygons, centers, node_coords, mesh)


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


def experiment_metrics(output: BaselineOutput) -> dict[str, object]:
    edges = topology_edges(output.cells)
    areas_signed = np.asarray([polygon_area(poly) for poly in output.polygons], dtype=float)
    areas = np.abs(areas_signed)
    perimeters = np.asarray([float(np.sum(edge_lengths(poly))) for poly in output.polygons], dtype=float)
    aspect = []
    for poly in output.polygons:
        lengths = edge_lengths(poly)
        positive = lengths[lengths > 1e-12]
        aspect.append(float(np.max(positive) / np.min(positive)) if len(positive) else np.nan)
    aspect = np.asarray(aspect, dtype=float)
    vertices = np.vstack(output.polygons)
    vertex_radius = np.hypot(vertices[:, 0], vertices[:, 1])
    center_radius = np.hypot(output.centers[:, 0], output.centers[:, 1])
    edge_lengths_mapped = np.asarray([float(np.linalg.norm(output.centers[i] - output.centers[j])) for i, j in edges], dtype=float)
    overlap_pairs, overlap_area = non_adjacent_overlap(output.polygons, edges)
    area_sum = float(np.sum(areas))
    return {
        "experiment": output.key,
        "title": output.title,
        "kind": output.kind,
        "cell_count": len(output.cells),
        "topology_edges": len(edges),
        "connected_components": connected_components(len(output.cells), edges),
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
        **classify_edges_by_layers(output.centers, edges),
    }


def direction_labels(ax, radius: float) -> None:
    for angle, label, ha, va in [
        (math.pi / 2.0, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3.0 * math.pi / 2.0, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ]:
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=10, color="#444444")


def save_plot(output: BaselineOutput, output_path: Path, *, draw_edges: bool = False) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.6, 9.6), dpi=260)
    fig.patch.set_facecolor("white")
    theta = np.linspace(0.0, 2.0 * math.pi, 721)
    ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.add_collection(PolyCollection(output.polygons, facecolors="none", edgecolors="#5d5d5d", linewidths=0.34, zorder=2))
    if draw_edges:
        segments = [(output.centers[i], output.centers[j]) for i, j in topology_edges(output.cells)]
        ax.add_collection(LineCollection(segments, colors="#c9483a", linewidths=0.48, alpha=0.82, zorder=3))
    ax.scatter(output.centers[:, 0], output.centers[:, 1], s=4.2, color="#111111", linewidths=0, zorder=4)
    for xy, cell in zip(output.centers, output.cells):
        ax.text(float(xy[0]), float(xy[1]), str(cell.cell_id), ha="center", va="center", fontsize=2.7, color="#111111", zorder=5)
    direction_labels(ax, OUTER_RADIUS + 0.13)
    ax.set_title(output.title, loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    coords = np.vstack(output.polygons)
    limit = max(float(np.abs(coords).max()), OUTER_RADIUS) + 0.22
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.20)
    fig.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def save_cells_csv(output: BaselineOutput, path: Path) -> None:
    rows: list[dict[str, object]] = []
    for cell, center, poly in zip(output.cells, output.centers, output.polygons):
        area_signed = polygon_area(poly)
        rows.append(
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
    write_csv(path, rows, ["cell_id", "row", "col", "center_x", "center_y", "vertex_count", "signed_area", "area"])


def save_vertices_csv(output: BaselineOutput, path: Path) -> None:
    rows: list[dict[str, object]] = []
    for node in sorted(output.node_coords):
        xy = output.node_coords[node]
        rows.append(
            {
                "node": repr(node),
                "node_type": node[0],
                "is_boundary": int(node in output.mesh.boundary_nodes),
                "x": f"{float(xy[0]):.10f}",
                "y": f"{float(xy[1]):.10f}",
                "r": f"{float(np.hypot(xy[0], xy[1])):.10f}",
                "theta_rad": f"{float(np.mod(math.atan2(float(xy[1]), float(xy[0])), 2.0 * math.pi)):.10f}",
            }
        )
    write_csv(path, rows, ["node", "node_type", "is_boundary", "x", "y", "r", "theta_rad"])


def save_standard_outputs(output: BaselineOutput, directory: Path, stem: str, readme_text: str) -> dict[str, object]:
    directory.mkdir(parents=True, exist_ok=True)
    save_plot(output, directory / f"{stem}_result.png")
    save_plot(output, directory / f"{stem}_topology_overlay.png", draw_edges=True)
    save_cells_csv(output, directory / f"{stem}_cells.csv")
    save_vertices_csv(output, directory / f"{stem}_vertices.csv")
    metric_row = experiment_metrics(output)
    write_csv(directory / f"{stem}_metrics.csv", [metric_row], list(metric_row.keys()))
    (directory / "README.md").write_text(readme_text.format(**metric_row), encoding="utf-8")
    return metric_row


def optimize_area_preserving(mesh: AnnulusMesh, initial_coords: dict[tuple, np.ndarray]) -> dict[tuple, np.ndarray]:
    fixed = boundary_positions(mesh)
    free_nodes = [node for node in sorted(mesh.source_node_coords) if node not in fixed]
    index = {node: i for i, node in enumerate(free_nodes)}
    edges = sorted(mesh_edges(mesh.polygons))
    x0 = np.concatenate([initial_coords[node] for node in free_nodes])
    initial_polygons = polygons_from_node_coords(mesh, initial_coords)
    target_area = float(np.mean([abs(polygon_area(poly)) for poly in initial_polygons]))
    initial_edge_vectors = {edge: initial_coords[edge[0]] - initial_coords[edge[1]] for edge in edges}

    def unpack(x: np.ndarray) -> dict[tuple, np.ndarray]:
        coords = dict(fixed)
        for node, i in index.items():
            coords[node] = x[2 * i : 2 * i + 2]
        return coords

    def residuals(x: np.ndarray) -> np.ndarray:
        coords = unpack(x)
        residual: list[float] = []
        for polygon in mesh.polygons:
            poly = np.asarray([coords[node] for node in polygon], dtype=float)
            area = max(abs(polygon_area(poly)), 1e-10)
            residual.append(1.25 * math.log(area / target_area))
        for a, b in edges:
            if a in fixed and b in fixed:
                continue
            current = coords[a] - coords[b]
            residual.extend((0.55 * (current - initial_edge_vectors[(a, b)])).tolist())
        for node in free_nodes:
            residual.extend((0.18 * (coords[node] - initial_coords[node])).tolist())
        return np.asarray(residual, dtype=float)

    result = least_squares(residuals, x0, max_nfev=400, ftol=1e-7, xtol=1e-7, gtol=1e-7, verbose=0)
    return unpack(result.x)
