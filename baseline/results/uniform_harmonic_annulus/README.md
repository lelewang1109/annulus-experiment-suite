# Baseline 2 - 直接圆环调和参数化

该 baseline 对原始北京共享顶点网格构造内边界和外边界，然后固定双边界并求解均匀权重图 Laplacian。

流程：

1. 从 `data/input/beijing_grid_cells.csv` 重建 176 个共享顶点 cell。
2. 以 cell 96、97、110、111 的共享顶点构造一个小内边界。
3. 将外边界顶点按原始角度顺序固定到半径 `1.85` 的圆上。
4. 将内边界顶点按 E/N/W/S 顺序固定到半径 `1.0` 的圆上。
5. 对所有内部顶点求解均匀权重调和方程 `Lu = 0`。

该方法对应“直接圆环调和参数化”，重点是固定双边界并让内部顶点平滑插值，不加入保角或保面积优化。

## 输出文件

- `baseline_harmonic_annulus_result.png`
- `baseline_harmonic_annulus_topology_overlay.png`
- `baseline_harmonic_annulus_cells.csv`
- `baseline_harmonic_annulus_vertices.csv`
- `baseline_harmonic_annulus_metrics.csv`
- `code/uniform_harmonic_annulus.py`

## 当前指标摘要

- cell_count: 176
- topology_edges: 312
- vertices_outside_outer: 0
- vertices_inside_inner: 136
- non_adjacent_overlap_pairs: 46
- non_adjacent_overlap_area: 1.7224847577
- area_cv: 1.435614
- ring_area_coverage_proxy: 1.672514
