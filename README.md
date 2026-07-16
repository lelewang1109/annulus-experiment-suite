# Annulus Experiment

面向时空污染比较的地理网格环形抽象实验。仓库以北京市 176 个编号网格为案例，比较四种连续参数化 baseline，并整理三个拓扑感知 workflow 版本。

该设计的目标不是用圆环替代地图，而是建立稳定的比较坐标：同一单元在不同月份和污染阶段中保持几何位置，只更换数据编码；不同区域则共享方向、中心到外围层次和数值编码语法。

## Repository layout

```text
.
├── baseline/       # 四个独立 baseline 实验，按方法名组织
├── workflow-v1/    # workflow v1：第一版拓扑感知优化
├── workflow-v2/    # workflow v2：中心点输入的完整结果快照
└── workflow-v3/    # workflow v3：0715 四角点/四顶点补充实验
```

每个实验目录都自包含：

```text
<experiment>/
├── data/       # 该实验需要的输入或生成数据
├── code/       # 该实验对应代码
├── results/    # 该实验对应结果
└── README.md   # 该实验内部说明
```


## Baseline

Baseline experiments 用来说明直接连续变形为什么不足以稳定保留拓扑、方向和面积关系。

| Method | Code | Results | Key metrics |
|---|---|---|---|
| Direct four-corner polar mapping | `baseline/code/direct_polar_mapping.py` | `baseline/results/direct_polar_mapping/` | Area CV 1.212; overlaps 3; inner violations 0; max edge ratio 13.50 |
| Uniform harmonic annulus | `baseline/code/uniform_harmonic_annulus.py` | `baseline/results/uniform_harmonic_annulus/` | Area CV 1.436; overlaps 46; inner violations 136; max edge ratio 15.31 |
| Cotangent/conformal annulus | `baseline/code/cotangent_conformal_annulus.py` | `baseline/results/cotangent_conformal_annulus/` | Area CV 1.466; overlaps 46; inner violations 128; max edge ratio 8.67 |
| Area-preserving optimization | `baseline/code/area_preserving_annulus.py` | `baseline/results/area_preserving_annulus/` | Area CV 0.210; overlaps 250; inner violations 268; max edge ratio 33.69 |

## Workflow v1

Workflow v1 保留第一版 topology-aware optimizer。它使用中心极坐标 cell 输入和 Step 3 拓扑约束，再优化 layer assignment、layer order 和 angular boundaries。

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Center-polar reference | `workflow-v1/code/step01_build_center_polar_annulus.py` | `workflow-v1/results/step01_center_polar/` |
| 02 | Center topology overlay | `workflow-v1/code/step02_build_center_topology.py` | `workflow-v1/results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `workflow-v1/code/step03_build_line_md_scaffold.py` | `workflow-v1/results/step03_line_md_scaffold/` |
| 04 | Topology-aware v1 optimizer | `workflow-v1/code/step04_build_topology_aware_v1.py` | `workflow-v1/results/step04_topology_aware_v1/` |

`workflow-v1/results/step04_initial_partition/` 保留早期 Step 04 initial partition 快照。Step 04 v1 topology-aware optimizer 的完整输出快照未提交；对应结果目录中保留了重新生成说明。

## Workflow v2

Workflow v2 是中心点输入下的 topology-aware 完整结果快照。

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Center-polar reference | `workflow-v2/code/step01_build_center_polar_annulus.py` | `workflow-v2/results/step01_center_polar/` |
| 02 | Center topology overlay | `workflow-v2/code/step02_build_center_topology.py` | `workflow-v2/results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `workflow-v2/code/step03_build_line_md_scaffold.py` | `workflow-v2/results/step03_line_md_scaffold/` |
| 04 | Topology-aware optimizer | `workflow-v2/code/step04_build_topology_aware.py` | `workflow-v2/results/step04_topology_aware/` |

| Metric | Value |
|---|---:|
| Precision | 0.595 |
| Recall | 0.763 |
| F1 | 0.669 |
| Preserved original edges | 238 |
| Lost original edges | 74 |
| New edges | 162 |
| Long-layer conflicts | 14 |
| Layer capacity error | 0.00675 |

## Workflow v3

Workflow v3 是 0715 四角点/四顶点补充实验。它改变 Step 2 之前的几何入口：每个 grid cell 先由四个经纬度角点重建，每个角点投影到圆环，再使用投影后 polygon centroid 进入后续 topology-aware workflow。

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Four-vertex polar mapping | `workflow-v3/code/step01_build_four_vertex_polar.py` | `workflow-v3/results/step01_four_vertex_polar/` |
| 02 | Center topology overlay from four-vertex centroids | `workflow-v3/code/step02_build_center_topology.py` | `workflow-v3/results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `workflow-v3/code/step03_build_line_md_scaffold.py` | `workflow-v3/results/step03_line_md_scaffold/` |
| 04 | Topology-aware optimizer on v3 input | `workflow-v3/code/step04_build_topology_aware.py` | `workflow-v3/results/step04_topology_aware/` |

Generated v3 data:

- `workflow-v3/data/generated/four_vertex_grid_cells.csv`
- `workflow-v3/data/generated/four_vertex_grid_corners.csv`
- `workflow-v3/results/step01_four_vertex_polar/step01_four_vertex_metrics.csv`

| Metric | Value |
|---|---:|
| Precision | 0.632 |
| Recall | 0.814 |
| F1 | 0.711 |
| Preserved original edges | 254 |
| Lost original edges | 58 |
| New edges | 148 |
| Long-layer conflicts | 11 |
| Layer capacity error | 0.00619 |

## Environment

- Python 3.10+
- NumPy
- Matplotlib
- SciPy

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Data boundary

各实验目录下的 `data/input/beijing_grid_cells.csv` 是对应实验的规范数据接口，包含 176 个单元的稳定 ID、行列索引、经纬度中心和中心极坐标。`beijing_boundary_rays.csv` 只保留直接极坐标 baseline 所需的 475 条角点射线距离，用于取代原先的本机 GeoJSON 绝对路径。原始 PM2.5/GeoJSON 数据准备工程不包含在此便携实验仓库中。
