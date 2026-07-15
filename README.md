# Annulus Experiment Suite

面向时空污染比较的地理网格环形抽象实验。仓库以北京市 176 个编号网格为案例，比较四种连续参数化 baseline，并实现一套从中心极坐标参考、原始邻接提取、径向分层到拓扑感知环形分区的可复现流程。

该设计的目标不是用圆环替代地图，而是建立稳定的比较坐标：同一单元在不同月份和污染阶段中保持几何位置，只更换数据编码；不同区域则共享方向、中心—外围层次和数值编码语法。


## Repository layout

```text
.
├── baseline/                  # baseline 四个独立实验：按方法名整理
├── workflow-v1/               # workflow v1：第一版拓扑感知优化
├── workflow-v2/               # workflow v2：中心点输入的完整结果快照
├── workflow-v3/               # workflow v3：0715 四角点/四顶点补充实验
└── EXPERIMENT_INDEX.md        # baseline 与 workflow-v1/v2/v3 的代码和结果索引
```

## Experiment map

本仓库包含两类实验：

- **Baseline experiments**：四种连续参数化或几何优化方法，整理在 `baseline/`。
- **Workflow experiments**：三个离散拓扑感知 workflow 版本，分别整理在 `workflow-v1/`、`workflow-v2/`、`workflow-v3/`。

先读 `EXPERIMENT_INDEX.md` 可以看到每个实验对应的代码入口、输入数据、结果目录和关键指标。`WORKFLOW_VERSIONS.md` 保留为 workflow 三个版本的简要说明。

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

## Reproduce

每个实验目录都自包含 `data/`、`code/` 和 `results/`。进入对应目录后运行该目录 `README.md` 中列出的代码入口即可复现相应结果。例如：

```bash
cd workflow-v3
python code/step01_build_four_vertex_polar.py
```

baseline 是四个独立方法，不按 step 命名；workflow-v1、workflow-v2、workflow-v3 则按 step 顺序组织代码和结果。根目录不再保留旧版 `data/`、`results/`、`scripts/`，避免和四个实验文件夹的最终结构混淆。

## Validate

GitHub Actions 会检查四个实验目录是否完整、关键说明文件是否存在，并阻止 `.DS_Store`、`.pyc` 和本机绝对路径进入公开仓库。

## Current result snapshot

| Method | Area CV | Non-adjacent overlaps | Inner-radius violations | Maximum edge ratio |
|---|---:|---:|---:|---:|
| Direct polar mapping | 1.212 | 3 | 0 | 13.50 |
| Uniform harmonic | 1.436 | 46 | 136 | 15.31 |
| Cotangent/conformal approximation | 1.466 | 46 | 128 | 8.67 |
| Area-preserving optimization | 0.210 | 250 | 268 | 33.69 |

拓扑感知 v2 的有效邻接 Precision / Recall / F1 为 `0.595 / 0.763 / 0.669`；保留 238 条原始边，丢失 74 条，新增 162 条。这些结果是多目标折中，不表示拓扑或方向已被完全保持。

## Data and version boundary

各实验目录下的 `data/input/beijing_grid_cells.csv` 是对应实验的规范数据接口，包含 176 个单元的稳定 ID、行列索引、经纬度中心和中心极坐标。`beijing_boundary_rays.csv` 只保留直接极坐标 baseline 所需的 475 条角点射线距离，用于取代原先的本机 GeoJSON 绝对路径。原始 PM2.5/GeoJSON 数据准备工程不包含在此便携实验仓库中。
