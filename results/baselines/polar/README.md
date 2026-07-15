# Baseline 1 - Polar Mapping

这个 baseline 将 `data/input/beijing_grid_cells.csv` 中的 176 个北京网格 cell 作为输入，并使用 `data/input/beijing_boundary_rays.csv` 中的预计算边界射线距离。

做法：

1. 用每个 cell 的中心经纬度和网格步长重建原始四个角点。
2. 以 cell 96、97、110、111 的共享中心顶点作为极坐标原点。
3. 每个角点单独计算 `theta` 和相对北京边界射线距离的 `rho_unclipped`。
4. 直接映射到圆环半径 `R = 1.0 + rho_unclipped * (1.85 - 1.0)`。

该方法不做拓扑修正、不做层列重排、不裁剪越界角点，因此用于观察“直接四角点极坐标映射”的原始几何问题。

## 输出文件

- `baseline_polar_mapping_result.png`: 四角点直接映射后的圆环网格。
- `baseline_polar_mapping_topology_overlay.png`: 叠加原始上下左右邻接边的结果图。
- `baseline_polar_mapping_corners.csv`: 每个 cell 四个角点的原始经纬度、theta、rho、圆环坐标。
- `baseline_polar_mapping_cells.csv`: 每个 cell 的圆环中心、面积和顶点数。
- `baseline_polar_mapping_metrics.csv`: 单独指标表。
- `scripts/baselines/build_baseline_polar_mapping.py`: 复现脚本。

## 当前指标摘要

- cell_count: 176
- topology_edges: 312
- vertices_outside_outer: 64
- vertices_inside_inner: 0
- non_adjacent_overlap_pairs: 3
- non_adjacent_overlap_area: 0.0012516596
- area_cv: 1.211691
- ring_area_coverage_proxy: 1.128473
