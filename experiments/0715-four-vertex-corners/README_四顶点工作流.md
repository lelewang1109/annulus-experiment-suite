# 0715-角点 workflow

这个副本把 `scripts/workflow` 的第一步改成了直接四顶点极坐标映射。

## 改动

- `scripts/workflow/step01_build_center_polar_annulus.py`
  - 不再把 cell 中心点映射到圆环后画固定正方形。
  - 改为用原始中心经纬度和网格步长重建每个 cell 的四个角点。
  - 四个角点分别按 `theta / rho / R` 映射到圆环。
  - `results/workflow/step01_center_polar.png` 现在显示四顶点映射后的多边形。
  - 额外写出：
    - `data/generated/four_vertex_grid_cells.csv`
    - `data/generated/four_vertex_grid_corners.csv`
    - `results/workflow/step01_four_vertex_metrics.csv`

- `scripts/workflow/step02_build_center_topology.py`
  - 默认读取 `data/generated/four_vertex_grid_cells.csv`。
  - 优先用 `data/generated/four_vertex_grid_corners.csv` 画真实四顶点多边形。

- `scripts/workflow/step03_build_line_md_scaffold.py`
- `scripts/workflow/step04_build_line_md_layer_column_partition.py`
- `scripts/workflow/step04_build_topology_aware_partition.py`
- `scripts/workflow/step04_build_topology_aware_partition_v2.py`
  - 默认输入改为 `data/generated/four_vertex_grid_cells.csv`，即用四顶点多边形质心的 `X_center / Y_center / theta / R` 继续推进。

## 已运行结果

- step01: `results/workflow/step01_center_polar.png`
- step02: `results/workflow/step02_center_topology.png`
- step03: `results/workflow/step03/step03_line_md_scaffold.png`
- step04 initial: `results/workflow/step04_initial/step04_line_md_layer_column_partition.png`
- step04 topology-aware v2: `results/workflow/step04_topology_aware_v2/step04_topology_aware_partition_v2.png`

## 当前 v2 指标

- baseline effective F1: `0.5755968169761273`
- optimized effective F1: `0.711484593837535`
- optimized preserved/lost/new: `254 / 58 / 148`
- remaining long-layer conflicts: `11`
