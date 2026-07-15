#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "baselines" / "area_preserving"
BASELINE_ROOT = THIS_DIR
if str(BASELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(BASELINE_ROOT))

from baseline_common import (  # noqa: E402
    INNER_RADIUS,
    OUTER_RADIUS,
    build_punctured_mesh,
    cotangent_weights,
    optimize_area_preserving,
    output_from_node_coords,
    save_standard_outputs,
    solve_harmonic,
)

STEM = "baseline_area_preserving_annulus"


def run() -> None:
    mesh = build_punctured_mesh()
    conformal_coords = solve_harmonic(mesh, cotangent_weights(mesh))
    node_coords = optimize_area_preserving(mesh, conformal_coords)
    output = output_from_node_coords(
        STEM,
        "Baseline 4 - Area-Preserving Annulus Parameterization",
        "area_preserving_optimized",
        mesh,
        node_coords,
    )
    readme = f"""# Baseline 4 - 保面积圆环参数化

该 baseline 先得到保角圆环参数化，再以其为初始化做单元面积均衡优化。

流程：

1. 重建北京共享顶点网格并构造中心内边界。
2. 外边界固定到半径 `{OUTER_RADIUS}` 的圆上，内边界固定到半径 `{INNER_RADIUS}` 的圆上。
3. 先用 cotangent Laplacian 得到保角初始化。
4. 以保角初始布局中单元绝对面积的均值作为 `A_target`，最小化 `log(A_i / A_target)` 的面积残差。
5. 同时加入弱边长保持和弱初始位置残差。当前目标不含显式翻转、重叠或内孔障碍项，因此面积均衡不代表几何可行。

该方法对应“保面积圆环参数化” baseline：目标是降低单元面积变化，代价是可能牺牲局部角度和形状。

## 输出文件

- `{STEM}_result.png`
- `{STEM}_topology_overlay.png`
- `{STEM}_cells.csv`
- `{STEM}_vertices.csv`
- `{STEM}_metrics.csv`
- `scripts/baselines/build_baseline_area_preserving_annulus.py`

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
