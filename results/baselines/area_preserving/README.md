# Baseline 4 - 保面积圆环参数化

该 baseline 先得到保角圆环参数化，再以其为初始化做单元面积均衡优化。

流程：

1. 重建北京共享顶点网格并构造中心内边界。
2. 外边界固定到半径 `1.85` 的圆上，内边界固定到半径 `1.0` 的圆上。
3. 先用 cotangent Laplacian 得到保角初始化。
4. 以保角初始布局中单元绝对面积的均值作为 `A_target`，最小化 `log(A_i / A_target)` 的面积残差。
5. 同时加入弱边长保持和弱初始位置残差。当前目标不含显式翻转、重叠或内孔障碍项，因此面积均衡不代表几何可行。

该方法对应“保面积圆环参数化” baseline：目标是降低单元面积变化，代价是可能牺牲局部角度和形状。

## 输出文件

- `baseline_area_preserving_annulus_result.png`
- `baseline_area_preserving_annulus_topology_overlay.png`
- `baseline_area_preserving_annulus_cells.csv`
- `baseline_area_preserving_annulus_vertices.csv`
- `baseline_area_preserving_annulus_metrics.csv`
- `scripts/baselines/build_baseline_area_preserving_annulus.py`

## 当前指标摘要

- cell_count: 176
- topology_edges: 312
- vertices_outside_outer: 0
- vertices_inside_inner: 268
- non_adjacent_overlap_pairs: 250
- non_adjacent_overlap_area: 4.7177744476
- area_cv: 0.210418
- ring_area_coverage_proxy: 1.794589
