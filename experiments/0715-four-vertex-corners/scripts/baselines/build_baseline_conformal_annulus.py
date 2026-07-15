#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "baselines" / "conformal"
BASELINE_ROOT = THIS_DIR
if str(BASELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(BASELINE_ROOT))

from baseline_common import (  # noqa: E402
    INNER_RADIUS,
    OUTER_RADIUS,
    build_punctured_mesh,
    cotangent_weights,
    output_from_node_coords,
    save_standard_outputs,
    solve_harmonic,
)

STEM = "baseline_conformal_annulus"


def run() -> None:
    mesh = build_punctured_mesh()
    node_coords = solve_harmonic(mesh, cotangent_weights(mesh))
    output = output_from_node_coords(
        STEM,
        "Baseline 3 - Conformal Annulus Parameterization",
        "cotangent_conformal_harmonic",
        mesh,
        node_coords,
    )
    readme = f"""# Baseline 3 - 保角圆环参数化

该 baseline 使用 cotangent 权重替代均匀权重，在固定内外圆边界后求解离散调和映射。

流程：

1. 从 176 个北京网格 cell 重建共享顶点网格，并在中心共享顶点处构造内边界。
2. 外边界按原始角度顺序固定到半径 `{OUTER_RADIUS}` 的圆上。
3. 内边界按 E/N/W/S 固定到半径 `{INNER_RADIUS}` 的圆上。
4. 将四边形/五边形 cell 扇形三角化，计算 cotangent 权重。
5. 对内部顶点求解 cotangent Laplacian `Lu = 0`，作为离散保角近似。

该方法对应“保角圆环参数化” baseline：目标是减少局部角度畸变，但不额外控制面积。

## 输出文件

- `{STEM}_result.png`
- `{STEM}_topology_overlay.png`
- `{STEM}_cells.csv`
- `{STEM}_vertices.csv`
- `{STEM}_metrics.csv`
- `scripts/baselines/build_baseline_conformal_annulus.py`

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
