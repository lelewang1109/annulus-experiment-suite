# Workflow v2

Workflow v2 is the center-point topology-aware workflow snapshot. It keeps the
original center-polar geometry entry and records the complete Step 1-5 result
chain.

## Code and results

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Center-polar reference | `code/step01_build_center_polar_annulus.py` | `results/step01_center_polar/` |
| 02 | Original topology over center-polar points | `code/step02_build_center_topology.py` | `results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `code/step03_build_line_md_scaffold.py` | `results/step03_line_md_scaffold/` |
| 04 | Initial layer-column partition | `code/step04_build_initial_partition.py` | `results/step04_initial_partition/` |
| 05 | Topology-aware v2 optimizer | `code/step05_build_topology_aware_v2.py` | `results/step05_topology_aware_v2/` |

`code/step04_build_topology_aware_partition.py` is kept with its original file
name because Step 05 imports helper functions from it.

## Optimized metrics

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
