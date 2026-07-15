#!/usr/bin/env python3
"""Render workflow step 01 with direct four-vertex polar mapping.

The original workflow mapped each source cell center to the annulus and then
rendered a fixed square around that mapped center.  This variant reconstructs
the four original grid corners for each source cell, maps every corner through
the same polar boundary-ray normalization, and renders the mapped quadrilateral
directly.  It also writes a generated cell-level CSV whose center fields are
the centroids of those mapped four-vertex polygons so later workflow stages can
advance from the four-vertex geometry instead of the old center mapping.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_INPUT_CSV = PROJECT_ROOT / "data" / "input" / "beijing_grid_cells.csv"
DEFAULT_BOUNDARY_RAYS_CSV = PROJECT_ROOT / "data" / "input" / "beijing_boundary_rays.csv"
DEFAULT_OUTPUT_PNG = PROJECT_ROOT / "results" / "workflow" / "step01_center_polar.png"
DEFAULT_OUTPUT_CELLS_CSV = PROJECT_ROOT / "data" / "generated" / "four_vertex_grid_cells.csv"
DEFAULT_OUTPUT_CORNERS_CSV = PROJECT_ROOT / "data" / "generated" / "four_vertex_grid_corners.csv"
DEFAULT_OUTPUT_METRICS_CSV = PROJECT_ROOT / "results" / "workflow" / "step01_four_vertex_metrics.csv"


@dataclass(frozen=True)
class SourceCell:
    cell_id: int
    row: int
    col: int
    center_lon: float
    center_lat: float
    source_theta: float
    source_rho: float
    source_radius: float
    source_x: float
    source_y: float
    source_side: float


@dataclass
class FourVertexResult:
    cells: list[SourceCell]
    source_polygons: list[np.ndarray]
    annulus_polygons: list[np.ndarray]
    centers: np.ndarray
    center_theta: np.ndarray
    center_rho: np.ndarray
    center_radius: np.ndarray
    corner_theta: list[np.ndarray]
    corner_rho: list[np.ndarray]
    corner_radius: list[np.ndarray]


def read_cells(path: Path) -> list[SourceCell]:
    records: list[SourceCell] = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            records.append(
                SourceCell(
                    cell_id=int(row["cell_id"]),
                    row=int(row["row"]),
                    col=int(row["col"]),
                    center_lon=float(row["center_lon"]),
                    center_lat=float(row["center_lat"]),
                    source_theta=float(row["theta_rad"]),
                    source_rho=float(row["rho"]),
                    source_radius=float(row["R"]),
                    source_x=float(row["X_center"]),
                    source_y=float(row["Y_center"]),
                    source_side=float(row["fixed_cell_side"]),
                )
            )
    if not records:
        raise ValueError(f"No cell records found in {path}")
    return records


def read_boundary_rays(path: Path) -> dict[tuple[float, float], float]:
    lookup: dict[tuple[float, float], float] = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (round(float(row["source_lon"]), 10), round(float(row["source_lat"]), 10))
            lookup[key] = float(row["boundary_radius"])
    if not lookup:
        raise ValueError(f"No boundary ray records found in {path}")
    return lookup


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def infer_grid_steps(cells: list[SourceCell]) -> tuple[float, float]:
    lons = np.asarray(sorted({round(cell.center_lon, 6) for cell in cells}), dtype=float)
    lats = np.asarray(sorted({round(cell.center_lat, 6) for cell in cells}), dtype=float)
    lon_diffs = np.diff(lons)
    lat_diffs = np.diff(lats)
    return float(np.nanmedian(lon_diffs[lon_diffs > 1e-8])), float(np.nanmedian(lat_diffs[lat_diffs > 1e-8]))


def source_cell_polygon(cell: SourceCell, lon_step: float, lat_step: float) -> np.ndarray:
    half_lon = lon_step / 2.0
    half_lat = lat_step / 2.0
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
    return float(np.mean([cell.center_lon for cell in cells])), float(np.mean([cell.center_lat for cell in cells]))


def project_lonlat(points: np.ndarray, center_lon: float, center_lat: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return np.column_stack(
        [
            (points[:, 0] - center_lon) * math.cos(math.radians(center_lat)),
            points[:, 1] - center_lat,
        ]
    )


def map_polygon_to_annulus(
    polygon_lonlat: np.ndarray,
    *,
    center_lon: float,
    center_lat: float,
    boundary_rays: dict[tuple[float, float], float],
    inner_radius: float,
    outer_radius: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    projected = project_lonlat(polygon_lonlat, center_lon, center_lat)
    theta = np.mod(np.arctan2(projected[:, 1], projected[:, 0]), 2.0 * math.pi)
    distance = np.hypot(projected[:, 0], projected[:, 1])
    keys = [(round(float(lon), 10), round(float(lat), 10)) for lon, lat in polygon_lonlat]
    try:
        boundary_radius = np.asarray([boundary_rays[key] for key in keys], dtype=float)
    except KeyError as exc:
        raise KeyError(f"Missing precomputed boundary ray for source-grid corner {exc.args[0]}") from exc
    rho = distance / np.maximum(boundary_radius, 1e-12)
    radius = inner_radius + rho * (outer_radius - inner_radius)
    mapped = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])
    return mapped, theta, rho, radius


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


def build_result(
    input_csv: Path,
    boundary_rays_csv: Path,
    *,
    inner_radius: float,
    outer_radius: float,
) -> FourVertexResult:
    cells = read_cells(input_csv)
    boundary_rays = read_boundary_rays(boundary_rays_csv)
    center_lon, center_lat = choose_origin_from_center_cells(cells)
    lon_step, lat_step = infer_grid_steps(cells)

    source_polygons: list[np.ndarray] = []
    annulus_polygons: list[np.ndarray] = []
    corner_theta: list[np.ndarray] = []
    corner_rho: list[np.ndarray] = []
    corner_radius: list[np.ndarray] = []
    for cell in cells:
        source_poly = source_cell_polygon(cell, lon_step, lat_step)
        mapped, theta, rho, radius = map_polygon_to_annulus(
            source_poly,
            center_lon=center_lon,
            center_lat=center_lat,
            boundary_rays=boundary_rays,
            inner_radius=inner_radius,
            outer_radius=outer_radius,
        )
        source_polygons.append(source_poly)
        annulus_polygons.append(mapped)
        corner_theta.append(theta)
        corner_rho.append(rho)
        corner_radius.append(radius)

    centers = np.asarray([polygon_centroid(poly) for poly in annulus_polygons], dtype=float)
    center_radius = np.hypot(centers[:, 0], centers[:, 1])
    center_theta = np.mod(np.arctan2(centers[:, 1], centers[:, 0]), 2.0 * math.pi)
    center_rho = (center_radius - inner_radius) / max(outer_radius - inner_radius, 1e-12)
    return FourVertexResult(cells, source_polygons, annulus_polygons, centers, center_theta, center_rho, center_radius, corner_theta, corner_rho, corner_radius)


def add_direction_labels(ax: plt.Axes, radius: float) -> None:
    for angle, label, ha, va in (
        (math.pi / 2.0, "N", "center", "bottom"),
        (0.0, "E", "left", "center"),
        (3.0 * math.pi / 2.0, "S", "center", "top"),
        (math.pi, "W", "right", "center"),
    ):
        ax.text(radius * math.cos(angle), radius * math.sin(angle), label, ha=ha, va=va, fontsize=10, color="#444444")


def render(result: FourVertexResult, output_png: Path, *, inner_radius: float, outer_radius: float) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.2, 10.2), dpi=280)
    theta = np.linspace(0.0, 2.0 * math.pi, 721)
    for radius in (inner_radius, outer_radius):
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="#222222", linewidth=0.9, zorder=1)
    ax.add_collection(PolyCollection(result.annulus_polygons, facecolors="none", edgecolors="#696969", linewidths=0.38, zorder=2))
    ax.scatter(result.centers[:, 0], result.centers[:, 1], s=4.5, color="#111111", linewidths=0, zorder=3)
    for cell, (x, y) in zip(result.cells, result.centers):
        ax.text(float(x), float(y), str(cell.cell_id), ha="center", va="center", fontsize=3.35, color="#111111", zorder=4)
    add_direction_labels(ax, outer_radius + 0.13)
    ax.set_title("Beijing Grid Four-Vertex Polar Mapping to Annulus", loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal", adjustable="box")
    limit = max(float(np.abs(np.vstack(result.annulus_polygons)).max()), outer_radius) + 0.35
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle="--", linewidth=0.35, alpha=0.22)
    fig.savefig(output_png, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def write_generated_data(
    result: FourVertexResult,
    *,
    output_cells_csv: Path,
    output_corners_csv: Path,
    output_metrics_csv: Path,
    inner_radius: float,
    outer_radius: float,
) -> None:
    corner_rows: list[dict[str, object]] = []
    for cell, source_poly, mapped_poly, theta, rho, radius in zip(
        result.cells,
        result.source_polygons,
        result.annulus_polygons,
        result.corner_theta,
        result.corner_rho,
        result.corner_radius,
    ):
        for corner_index, (source, xy, th, rr, rad) in enumerate(zip(source_poly, mapped_poly, theta, rho, radius), start=1):
            corner_rows.append(
                {
                    "cell_id": cell.cell_id,
                    "row": cell.row,
                    "col": cell.col,
                    "corner_index": corner_index,
                    "source_lon": f"{float(source[0]):.10f}",
                    "source_lat": f"{float(source[1]):.10f}",
                    "theta_rad": f"{float(th):.10f}",
                    "theta_deg": f"{math.degrees(float(th)):.6f}",
                    "rho": f"{float(rr):.10f}",
                    "R": f"{float(rad):.10f}",
                    "annulus_x": f"{float(xy[0]):.10f}",
                    "annulus_y": f"{float(xy[1]):.10f}",
                }
            )
    write_csv(
        output_corners_csv,
        corner_rows,
        ["cell_id", "row", "col", "corner_index", "source_lon", "source_lat", "theta_rad", "theta_deg", "rho", "R", "annulus_x", "annulus_y"],
    )

    cell_rows: list[dict[str, object]] = []
    signed_areas = []
    for idx, (cell, center, poly) in enumerate(zip(result.cells, result.centers, result.annulus_polygons)):
        signed_area = polygon_area(poly)
        signed_areas.append(signed_area)
        cell_rows.append(
            {
                "cell_id": cell.cell_id,
                "row": cell.row,
                "col": cell.col,
                "center_lon": f"{cell.center_lon:.10f}",
                "center_lat": f"{cell.center_lat:.10f}",
                "theta_rad": f"{float(result.center_theta[idx]):.10f}",
                "theta_deg": f"{math.degrees(float(result.center_theta[idx])):.6f}",
                "rho": f"{float(result.center_rho[idx]):.10f}",
                "R": f"{float(result.center_radius[idx]):.10f}",
                "X_center": f"{float(center[0]):.10f}",
                "Y_center": f"{float(center[1]):.10f}",
                "fixed_cell_side": f"{math.sqrt(abs(signed_area)):.10f}",
                "four_vertex_area": f"{abs(signed_area):.10f}",
                "four_vertex_signed_area": f"{signed_area:.10f}",
                "source_center_theta_rad": f"{cell.source_theta:.10f}",
                "source_center_R": f"{cell.source_radius:.10f}",
                "source_center_X": f"{cell.source_x:.10f}",
                "source_center_Y": f"{cell.source_y:.10f}",
            }
        )
    write_csv(
        output_cells_csv,
        cell_rows,
        [
            "cell_id",
            "row",
            "col",
            "center_lon",
            "center_lat",
            "theta_rad",
            "theta_deg",
            "rho",
            "R",
            "X_center",
            "Y_center",
            "fixed_cell_side",
            "four_vertex_area",
            "four_vertex_signed_area",
            "source_center_theta_rad",
            "source_center_R",
            "source_center_X",
            "source_center_Y",
        ],
    )

    areas = np.abs(np.asarray(signed_areas, dtype=float))
    vertices = np.vstack(result.annulus_polygons)
    vertex_radius = np.hypot(vertices[:, 0], vertices[:, 1])
    center_radius = np.hypot(result.centers[:, 0], result.centers[:, 1])
    metrics = {
        "method": "four_vertex_direct_polar",
        "cell_count": len(result.cells),
        "zero_area_cells": int(np.sum(areas <= 1e-10)),
        "clockwise_orientation_cells": int(np.sum(np.asarray(signed_areas) < -1e-10)),
        "vertices_outside_outer": int(np.sum(vertex_radius > outer_radius + 1e-8)),
        "vertices_inside_inner": int(np.sum(vertex_radius < inner_radius - 1e-8)),
        "centers_outside_outer": int(np.sum(center_radius > outer_radius + 1e-8)),
        "centers_inside_inner": int(np.sum(center_radius < inner_radius - 1e-8)),
        "area_sum": float(np.sum(areas)),
        "area_mean": float(np.mean(areas)),
        "area_cv": float(np.std(areas) / max(np.mean(areas), 1e-12)),
        "area_min": float(np.min(areas)),
        "area_max": float(np.max(areas)),
    }
    write_csv(output_metrics_csv, [metrics], list(metrics.keys()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--boundary-rays-csv", type=Path, default=DEFAULT_BOUNDARY_RAYS_CSV)
    parser.add_argument("--output-png", type=Path, default=DEFAULT_OUTPUT_PNG)
    parser.add_argument("--output-cells-csv", type=Path, default=DEFAULT_OUTPUT_CELLS_CSV)
    parser.add_argument("--output-corners-csv", type=Path, default=DEFAULT_OUTPUT_CORNERS_CSV)
    parser.add_argument("--output-metrics-csv", type=Path, default=DEFAULT_OUTPUT_METRICS_CSV)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=1.85)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.outer_radius <= args.inner_radius:
        raise ValueError("--outer-radius must be greater than --inner-radius")
    result = build_result(args.input_csv, args.boundary_rays_csv, inner_radius=args.inner_radius, outer_radius=args.outer_radius)
    render(result, args.output_png, inner_radius=args.inner_radius, outer_radius=args.outer_radius)
    write_generated_data(
        result,
        output_cells_csv=args.output_cells_csv,
        output_corners_csv=args.output_corners_csv,
        output_metrics_csv=args.output_metrics_csv,
        inner_radius=args.inner_radius,
        outer_radius=args.outer_radius,
    )
    print(f"cells: {len(result.cells)}")
    print(f"wrote: {args.output_png}")
    print(f"wrote: {args.output_cells_csv}")
    print(f"wrote: {args.output_corners_csv}")
