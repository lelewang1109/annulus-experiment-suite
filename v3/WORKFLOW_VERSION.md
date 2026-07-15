# v3 four-vertex workflow snapshot

v3 is the 0715 four-corner experiment. It is a supplement to the original
center-point workflow, not a replacement for it.

## What changed

Step 1 reconstructs each selected source grid cell from four original
longitude/latitude corners. Each corner is projected through the polar
boundary-ray normalization. The resulting four-vertex polygon is written as
generated data, and its centroid becomes the input for the later topology-aware
workflow.

## Code

- `scripts/workflow/step01_build_center_polar_annulus.py`
- `scripts/workflow/step02_build_center_topology.py`
- `scripts/workflow/step03_build_line_md_scaffold.py`
- `scripts/workflow/step04_build_line_md_layer_column_partition.py`
- `scripts/workflow/step04_build_topology_aware_partition.py`
- `scripts/workflow/step04_build_topology_aware_partition_v2.py`

## Generated data

- `data/generated/four_vertex_grid_cells.csv`
- `data/generated/four_vertex_grid_corners.csv`
- `results/workflow/step01_four_vertex_metrics.csv`

## Results

Primary result directory:

- `results/workflow/step04_topology_aware_v2/`

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

The metrics are stored in
`results/workflow/step04_topology_aware_v2/step04_topology_aware_metrics_v2.json`.
