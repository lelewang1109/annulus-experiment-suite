# Annulus Experiment Suite

面向时空污染比较的地理网格环形抽象实验。仓库以北京市 176 个编号网格为案例，比较四种连续参数化 baseline，并实现一套从中心极坐标参考、原始邻接提取、径向分层到拓扑感知环形分区的可复现流程。

该设计的目标不是用圆环替代地图，而是建立稳定的比较坐标：同一单元在不同月份和污染阶段中保持几何位置，只更换数据编码；不同区域则共享方向、中心—外围层次和数值编码语法。


## Repository layout

```text
.
├── data/
│   ├── input/                 # 便携、可版本化的实验输入
│   └── reference/             # 只用于解释和人工核对的参考文件
├── baseline/                  # baseline 四个独立实验：按方法名整理
├── workflow-v1/               # workflow v1：第一版拓扑感知优化
├── workflow-v2/               # workflow v2：中心点输入的完整结果快照
├── workflow-v3/               # workflow v3：0715 四角点/四顶点补充实验
├── scripts/
│   ├── baselines/             # 原始 baseline 代码
│   ├── workflow/              # 原始中心点 workflow 代码
│   ├── run_all.py             # 按依赖顺序复现实验
│   └── validate_repository.py # 仓库与结果完整性检查
├── results/
│   ├── baselines/             # 原始 baseline 结果
│   └── workflow/              # 原始中心点 workflow 结果
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

先执行快速流程，生成 baseline 1–3 和 workflow Step 1–4 的非慢速结果：

```bash
python scripts/run_all.py --skip-slow
```

执行全部实验（面积优化和拓扑感知 v2 耗时更长）：

```bash
python scripts/run_all.py
```

只复现一个阶段：

```bash
python scripts/run_all.py --only step04-v2
```

`run_all.py` 对 v2 显式锁定 `w_new=4.5`、`w_angle=1.2`、`w_angle_cap=3.0` 和 `max_rank_shift=2`，与已提交的 JSON 结果保持一致，不受研究过程中脚本默认值变化的影响。脚本自身的默认值也已同步到该可复现配置。

## Validate

```bash
python scripts/validate_repository.py
```

该检查会验证输入单元数量、关键结果文件、v2 指标字段，并阻止 `.DS_Store`、`.pyc` 和本机绝对路径进入公开仓库。

## Current result snapshot

| Method | Area CV | Non-adjacent overlaps | Inner-radius violations | Maximum edge ratio |
|---|---:|---:|---:|---:|
| Direct polar mapping | 1.212 | 3 | 0 | 13.50 |
| Uniform harmonic | 1.436 | 46 | 136 | 15.31 |
| Cotangent/conformal approximation | 1.466 | 46 | 128 | 8.67 |
| Area-preserving optimization | 0.210 | 250 | 268 | 33.69 |

拓扑感知 v2 的有效邻接 Precision / Recall / F1 为 `0.595 / 0.763 / 0.669`；保留 238 条原始边，丢失 74 条，新增 162 条。这些结果是多目标折中，不表示拓扑或方向已被完全保持。

## Data and version boundary

`data/input/beijing_grid_cells.csv` 是本仓库的规范数据接口，包含 176 个单元的稳定 ID、行列索引、经纬度中心和中心极坐标。`beijing_boundary_rays.csv` 只保留直接极坐标 baseline 所需的 475 条角点射线距离，用于取代原先的本机 GeoJSON 绝对路径。原始 PM2.5/GeoJSON 数据准备工程不包含在此便携实验仓库中。
