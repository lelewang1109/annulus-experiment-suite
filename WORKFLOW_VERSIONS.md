# Workflow versions

This repository keeps the original center-point workflow at the project root and
adds three explicit workflow version folders for clearer comparison.

## Version map

| Version | Geometry entry | Main code | Result snapshot | Status |
|---|---|---|---|---|
| v1 | Center-polar cells | `v1/scripts/workflow/step04_build_topology_aware_partition.py` | Not committed as a full result snapshot | First topology-aware optimizer; code preserved for reruns. |
| v2 | Center-polar cells | `v2/scripts/workflow/step04_build_topology_aware_partition_v2.py` | `v2/results/workflow/step04_topology_aware_v2/` | Current center-point workflow result. |
| v3 | Four-vertex polar cells | `v3/scripts/workflow/step01_build_center_polar_annulus.py` and `v3/scripts/workflow/step04_build_topology_aware_partition_v2.py` | `v3/results/workflow/step04_topology_aware_v2/` | 0715 four-corner experiment supplement. |

## v1

v1 is the first topology-aware annular partition optimizer. It keeps the
center-polar input table from `data/input/beijing_grid_cells.csv` and optimizes
layer assignment, layer order, and angular boundaries with the earlier objective
and metric payload.

The v1 code is preserved in `v1/scripts/workflow/`. The current repository does
not include a full committed v1 result snapshot because the public result set
was standardized around v2. To regenerate v1:

```bash
cd v1
python3 scripts/workflow/step04_build_topology_aware_partition.py
```

## v2

v2 is the current center-point workflow. It keeps Step 1 as a center-polar
mapping and writes the final topology-aware result under
`v2/results/workflow/step04_topology_aware_v2/`.

Current optimized effective adjacency metrics:

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

## v3

v3 is the 0715 four-corner experiment. It does not replace the center-point v2
workflow. Instead, it changes the geometry entry of Step 1: each source cell is
reconstructed from four original grid corners, each corner is projected by the
polar boundary-ray normalization, and the generated four-vertex polygon centroid
is used by the later topology-aware stages.

Generated v3 inputs:

- `v3/data/generated/four_vertex_grid_cells.csv`
- `v3/data/generated/four_vertex_grid_corners.csv`
- `v3/results/workflow/step01_four_vertex_metrics.csv`

Current optimized effective adjacency metrics:

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

## Root workflow

The root `scripts/`, `data/`, and `results/` directories are retained as the
original center-point experiment suite. The version folders are explicit
snapshots for comparing workflow variants without overwriting the original
project layout.
