# Workflow v1

Workflow v1 is the first topology-aware annular partition experiment. It uses
the center-polar input and the Step 3 topology scaffold, then runs the earlier
topology-aware optimizer.

## Code and results

| Step | Purpose | Code | Results |
|---|---|---|---|
| 01 | Center-polar reference | `code/step01_build_center_polar_annulus.py` | `results/step01_center_polar/` |
| 02 | Original topology over center-polar points | `code/step02_build_center_topology.py` | `results/step02_center_topology/` |
| 03 | Layer/edge scaffold | `code/step03_build_line_md_scaffold.py` | `results/step03_line_md_scaffold/` |
| 04 | Topology-aware v1 optimizer | `code/step04_build_topology_aware_v1.py` | `results/step04_topology_aware_v1/` |

`results/step04_initial_partition/` keeps the legacy Step 04 initial partition
snapshot. Step 04 v1 topology-aware code is preserved, but the full v1
topology-aware result snapshot is not committed. Use the README in
`results/step04_topology_aware_v1/` to regenerate it.
