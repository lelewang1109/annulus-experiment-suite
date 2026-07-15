# Workflow v3

Workflow v3 is the 0715 four-vertex/corner experiment. It keeps the later
topology-aware optimization logic, but changes Step 01 from center-point mapping
to four-vertex polar mapping.

## Code and results

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Four-vertex polar mapping | `code/step01_build_four_vertex_polar.py` | `results/step01_four_vertex_polar/` |
| 02 | Original topology over generated centroids | `code/step02_build_center_topology.py` | `results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `code/step03_build_line_md_scaffold.py` | `results/step03_line_md_scaffold/` |
| 04 | Initial layer-column partition | `code/step04_build_initial_partition.py` | `results/step04_initial_partition/` |
| 05 | Topology-aware v2 optimizer on v3 input | `code/step05_build_topology_aware_v2.py` | `results/step05_topology_aware_v2/` |

`code/step04_build_topology_aware_partition.py` is kept with its original file
name because Step 05 imports helper functions from it.

## Generated data

- `data/generated/four_vertex_grid_cells.csv`
- `data/generated/four_vertex_grid_corners.csv`
- `results/step01_four_vertex_polar/step01_four_vertex_metrics.csv`

## Optimized metrics

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
