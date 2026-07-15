#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "baselines" / "harmonic"
BASELINE_ROOT = THIS_DIR
if str(BASELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(BASELINE_ROOT))

from baseline_common import (  # noqa: E402
    INNER_RADIUS,
    OUTER_RADIUS,
    build_punctured_mesh,
    output_from_node_coords,
    save_standard_outputs,
    solve_harmonic,
    uniform_weights,
)

STEM = "baseline_harmonic_annulus"


def run() -> None:
    mesh = build_punctured_mesh()
    node_coords = solve_harmonic(mesh, uniform_weights(mesh))
    output = output_from_node_coords(
        STEM,
        "Baseline 2 - Direct Harmonic Annulus Parameterization",
        "uniform_laplacian_harmonic",
        mesh,
        node_coords,
    )
    readme = f"""# Baseline 2 - 直接圆环调和参数化

该 baseline 对原始北京共享顶点网格构造内边界和外边界，然后固定双边界并求解均匀权重图 Laplacian。

流程：

1. 从 `data/input/beijing_grid_cells.csv` 重建 176 个共享顶点 cell。
2. 以 cell 96、97、110、111 的共享顶点构造一个小内边界。
3. 将外边界顶点按原始角度顺序固定到半径 `{OUTER_RADIUS}` 的圆上。
4. 将内边界顶点按 E/N/W/S 顺序固定到半径 `{INNER_RADIUS}` 的圆上。
5. 对所有内部顶点求解均匀权重调和方程 `Lu = 0`。

该方法对应“直接圆环调和参数化”，重点是固定双边界并让内部顶点平滑插值，不加入保角或保面积优化。

## 输出文件

- `{STEM}_result.png`
- `{STEM}_topology_overlay.png`
- `{STEM}_cells.csv`
- `{STEM}_vertices.csv`
- `{STEM}_metrics.csv`
- `code/uniform_harmonic_annulus.py`

## 当前指标摘要

- cell_count: {{cell_count}}
- topology_edges: {{topology_edges}}
- vertices_outside_outer: {{vertices_outside_outer}}
- vertices_inside_inner: {{vertices_inside_inner}}
- non_adjacent_overlap_pairs: {{non_adjacent_overlap_pairs}}
- non_adjacent_overlap_area: {{non_adjacent_overlap_area:.10f}}
- area_cv: {{area_cv:.6f}}
- ring_area_coverage_proxy: {{ring_area_coverage_proxy:.6f}}
"""
    metric_row = save_standard_outputs(output, OUTPUT_DIR, STEM, readme)
    print(f"cells: {metric_row['cell_count']}")
    print(f"wrote: {OUTPUT_DIR}")


if __name__ == "__main__":
    run()
