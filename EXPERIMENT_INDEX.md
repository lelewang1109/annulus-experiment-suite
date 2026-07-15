# Experiment index

This project contains two experiment families:

1. **Baseline experiments**: four continuous annulus parameterization methods.
2. **Workflow experiments**: three versions of the discrete topology-aware workflow.

The root `scripts/`, `data/`, and `results/` folders keep the original public
experiment suite. The `v1/`, `v2/`, and `v3/` folders are version snapshots for
reading and comparing workflow variants without overwriting the root workflow.

## Directory map

```text
.
├── scripts/
│   ├── baselines/             # Canonical baseline code
│   └── workflow/              # Canonical center-point workflow code
├── results/
│   ├── baselines/             # Canonical baseline results
│   └── workflow/              # Canonical center-point workflow results
├── v1/                        # Workflow v1: first topology-aware optimizer
├── v2/                        # Workflow v2: center-point snapshot with results
└── v3/                        # Workflow v3: 0715 four-vertex/corner experiment
```

## Baseline experiments

The baseline group asks: what happens if the original grid geometry is mapped to
an annulus by continuous parameterization or geometry optimization?

| Baseline | Method | Code | Main result directory | Key metrics |
|---|---|---|---|---|
| B1 | Four-corner direct polar mapping | `scripts/baselines/build_baseline_polar_mapping.py` | `results/baselines/polar/` | Area CV 1.212; overlaps 3; inner violations 0; max edge ratio 13.50 |
| B2 | Uniform harmonic annulus parameterization | `scripts/baselines/build_baseline_harmonic_annulus.py` | `results/baselines/harmonic/` | Area CV 1.436; overlaps 46; inner violations 136; max edge ratio 15.31 |
| B3 | Cotangent/conformal harmonic approximation | `scripts/baselines/build_baseline_conformal_annulus.py` | `results/baselines/conformal/` | Area CV 1.466; overlaps 46; inner violations 128; max edge ratio 8.67 |
| B4 | Area-preserving optimization | `scripts/baselines/build_baseline_area_preserving_annulus.py` | `results/baselines/area_preserving/` | Area CV 0.210; overlaps 250; inner violations 268; max edge ratio 33.69 |

Each baseline result directory contains:

- a rendered result image,
- a topology overlay image,
- cell/vertex/corner CSV files when applicable,
- one `*_metrics.csv` file.

The baseline results are diagnostic. They are not the final recommended annular
layout because each method fails at least one core requirement: boundary
validity, non-overlap, readable cell shape, or topology preservation.

## Workflow experiments

The workflow group asks: can the annular layout be built as a discrete,
topology-aware visual indexing space rather than a continuous map deformation?

| Version | Experiment role | Geometry entry | Code | Result directory | Result status |
|---|---|---|---|---|---|
| v1 | First topology-aware optimizer | Center-polar cell table | `v1/scripts/workflow/step04_build_topology_aware_partition.py` | `v1/results/` | Code preserved; full output snapshot not committed |
| v2 | Current center-point workflow snapshot | Center-polar cell table | `v2/scripts/workflow/step04_build_topology_aware_partition_v2.py` | `v2/results/workflow/step04_topology_aware_v2/` | Full result snapshot committed |
| v3 | 0715 four-vertex/corner supplement | Four-vertex polar polygons, then generated centroids | `v3/scripts/workflow/step01_build_center_polar_annulus.py` and `v3/scripts/workflow/step04_build_topology_aware_partition_v2.py` | `v3/results/workflow/step04_topology_aware_v2/` | Full result snapshot committed |

### v1

v1 keeps the earliest topology-aware optimizer. It optimizes layer assignment,
layer order, and angular boundaries using center-polar cells and the Step 3
topology constraints.

Committed files:

- code: `v1/scripts/workflow/`
- input: `v1/data/input/`
- topology constraints: `v1/results/workflow/step03/step03_topology_constraints.csv`

No full v1 output snapshot is committed. To regenerate it, run from `v1/`:

```bash
python3 scripts/workflow/step04_build_topology_aware_partition.py
```

Default output:

- `v1/results/workflow/step04_topology_aware_v1/`

### v2

v2 is the center-point workflow snapshot. Step 1 uses the original center-polar
cell table, then Step 2-4 build topology constraints, an initial layer-slot
partition, and the topology-aware v2 optimized partition.

Important result files:

- `v2/results/workflow/step01_center_polar.png`
- `v2/results/workflow/step02_center_topology.png`
- `v2/results/workflow/step03/step03_topology_constraints.csv`
- `v2/results/workflow/step04_initial/step04_line_md_layer_column_partition.png`
- `v2/results/workflow/step04_topology_aware_v2/step04_topology_aware_partition_v2.png`
- `v2/results/workflow/step04_topology_aware_v2/step04_topology_aware_metrics_v2.json`

Optimized effective adjacency metrics:

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

### v3

v3 is the 0715 four-vertex/corner experiment. It keeps the later topology-aware
optimization logic, but changes the geometry entry before Step 2:

1. reconstruct each grid cell from four original lon/lat corners,
2. project each corner with the polar boundary-ray normalization,
3. write generated four-vertex cells and corners,
4. feed the four-vertex polygon centroids into the Step 2-4 workflow.

Important generated files:

- `v3/data/generated/four_vertex_grid_cells.csv`
- `v3/data/generated/four_vertex_grid_corners.csv`
- `v3/results/workflow/step01_four_vertex_metrics.csv`

Important result files:

- `v3/results/workflow/step01_center_polar.png`
- `v3/results/workflow/step02_center_topology.png`
- `v3/results/workflow/step03/step03_topology_constraints.csv`
- `v3/results/workflow/step04_initial/step04_line_md_layer_column_partition.png`
- `v3/results/workflow/step04_topology_aware_v2/step04_topology_aware_partition_v2.png`
- `v3/results/workflow/step04_topology_aware_v2/step04_topology_aware_metrics_v2.json`

Optimized effective adjacency metrics:

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

## Recommended reading order

1. `README.md`: project purpose and quick reproduction commands.
2. `EXPERIMENT_INDEX.md`: where each experiment's code and results live.
3. `results/README.md`: canonical result file list.
4. `v1/README.md`, `v2/README.md`, `v3/WORKFLOW_VERSION.md`: version-specific notes.
