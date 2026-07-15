# Experiment index

The project is organized into four top-level experiment folders:

```text
.
├── baseline/       # Four continuous/geometric baseline experiments
├── workflow-v1/    # Workflow experiment v1: first topology-aware optimizer
├── workflow-v2/    # Workflow experiment v2: center-point topology-aware snapshot
└── workflow-v3/    # Workflow experiment v3: 0715 four-vertex/corner supplement
```

Each experiment folder follows this pattern:

```text
<experiment>/
├── data/       # Inputs or generated data required by this experiment
├── code/       # Code for this experiment family
├── results/    # Result snapshots for this experiment family
└── README.md   # Local explanation and code/result map
```

## Baseline

Baseline experiments diagnose why direct continuous deformation is not enough.

| Method | Code | Results | Key metrics |
|---|---|---|---|
| Direct four-corner polar mapping | `baseline/code/direct_polar_mapping.py` | `baseline/results/direct_polar_mapping/` | Area CV 1.212; overlaps 3; inner violations 0; max edge ratio 13.50 |
| Uniform harmonic annulus | `baseline/code/uniform_harmonic_annulus.py` | `baseline/results/uniform_harmonic_annulus/` | Area CV 1.436; overlaps 46; inner violations 136; max edge ratio 15.31 |
| Cotangent/conformal annulus | `baseline/code/cotangent_conformal_annulus.py` | `baseline/results/cotangent_conformal_annulus/` | Area CV 1.466; overlaps 46; inner violations 128; max edge ratio 8.67 |
| Area-preserving optimization | `baseline/code/area_preserving_annulus.py` | `baseline/results/area_preserving_annulus/` | Area CV 0.210; overlaps 250; inner violations 268; max edge ratio 33.69 |

## Workflow v1

Workflow v1 preserves the first topology-aware optimizer. It uses center-polar
cell inputs and the Step 3 topology constraints, then optimizes layer assignment,
layer order, and angular boundaries.

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Center-polar reference | `workflow-v1/code/step01_build_center_polar_annulus.py` | `workflow-v1/results/step01_center_polar/` |
| 02 | Center topology overlay | `workflow-v1/code/step02_build_center_topology.py` | `workflow-v1/results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `workflow-v1/code/step03_build_line_md_scaffold.py` | `workflow-v1/results/step03_line_md_scaffold/` |
| 04 | Topology-aware v1 optimizer | `workflow-v1/code/step04_build_topology_aware_v1.py` | `workflow-v1/results/step04_topology_aware_v1/` |

The Step 04 v1 full output snapshot is not committed; the directory contains a
README explaining how to regenerate it.

## Workflow v2

Workflow v2 is the center-point topology-aware result snapshot.

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Center-polar reference | `workflow-v2/code/step01_build_center_polar_annulus.py` | `workflow-v2/results/step01_center_polar/` |
| 02 | Center topology overlay | `workflow-v2/code/step02_build_center_topology.py` | `workflow-v2/results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `workflow-v2/code/step03_build_line_md_scaffold.py` | `workflow-v2/results/step03_line_md_scaffold/` |
| 04 | Initial layer-column partition | `workflow-v2/code/step04_build_initial_partition.py` | `workflow-v2/results/step04_initial_partition/` |
| 05 | Topology-aware v2 optimizer | `workflow-v2/code/step05_build_topology_aware_v2.py` | `workflow-v2/results/step05_topology_aware_v2/` |

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

## Workflow v3

Workflow v3 is the 0715 four-vertex/corner experiment. It changes the geometry
entry before Step 2: each grid cell is reconstructed from four lon/lat corners,
each corner is projected to the annulus, and the resulting polygon centroid is
used in the later topology-aware workflow.

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Four-vertex polar mapping | `workflow-v3/code/step01_build_four_vertex_polar.py` | `workflow-v3/results/step01_four_vertex_polar/` |
| 02 | Center topology overlay from four-vertex centroids | `workflow-v3/code/step02_build_center_topology.py` | `workflow-v3/results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `workflow-v3/code/step03_build_line_md_scaffold.py` | `workflow-v3/results/step03_line_md_scaffold/` |
| 04 | Initial layer-column partition | `workflow-v3/code/step04_build_initial_partition.py` | `workflow-v3/results/step04_initial_partition/` |
| 05 | Topology-aware v2 optimizer on v3 input | `workflow-v3/code/step05_build_topology_aware_v2.py` | `workflow-v3/results/step05_topology_aware_v2/` |

Generated v3 data:

- `workflow-v3/data/generated/four_vertex_grid_cells.csv`
- `workflow-v3/data/generated/four_vertex_grid_corners.csv`
- `workflow-v3/results/step01_four_vertex_polar/step01_four_vertex_metrics.csv`

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
