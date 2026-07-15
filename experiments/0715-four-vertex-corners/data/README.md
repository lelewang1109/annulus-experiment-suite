# Data

## `input/beijing_grid_cells.csv`

规范输入，共 176 行。关键字段：

- `cell_id`: 稳定单元编号。
- `row`, `col`: 原始规则网格索引，用于重建四邻接。
- `center_lon`, `center_lat`: 单元中心经纬度。
- `theta_rad`, `theta_deg`, `rho`, `R`: 中心极坐标参考。
- `X_center`, `Y_center`: 环形参考中心。
- `fixed_cell_side`: Step 1/2 诊断图的固定小方格边长。

## `input/beijing_boundary_rays.csv`

包含 475 个原始网格角点在参考中心下的方向和北京边界射线距离。该表是为消除本机 `china.json` 绝对路径而保存的最小派生输入，能精确复现已提交的直接极坐标 baseline。

## `reference/`

参考图和锚点表只用于人工核对，默认脚本不依赖这些文件。
