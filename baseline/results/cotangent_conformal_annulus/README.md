# Baseline 3 - 保角圆环参数化

该 baseline 使用 cotangent 权重替代均匀权重，在固定内外圆边界后求解离散调和映射。

流程：

1. 从 176 个北京网格 cell 重建共享顶点网格，并在中心共享顶点处构造内边界。
2. 外边界按原始角度顺序固定到半径 `1.85` 的圆上。
3. 内边界按 E/N/W/S 固定到半径 `1.0` 的圆上。
4. 将四边形/五边形 cell 扇形三角化，计算 cotangent 权重。
5. 对内部顶点求解 cotangent Laplacian `Lu = 0`，作为离散保角近似。

该方法对应“保角圆环参数化” baseline：目标是减少局部角度畸变，但不额外控制面积。

## 输出文件

- `baseline_conformal_annulus_result.png`
- `baseline_conformal_annulus_topology_overlay.png`
- `baseline_conformal_annulus_cells.csv`
- `baseline_conformal_annulus_vertices.csv`
- `baseline_conformal_annulus_metrics.csv`
- `code/cotangent_conformal_annulus.py`

## 当前指标摘要

- cell_count: 176
- topology_edges: 312
- vertices_outside_outer: 0
- vertices_inside_inner: 128
- non_adjacent_overlap_pairs: 46
- non_adjacent_overlap_area: 1.6064807102
- area_cv: 1.465563
- ring_area_coverage_proxy: 1.672716
