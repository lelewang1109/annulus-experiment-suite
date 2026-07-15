# v2 center-point workflow snapshot

v2 is the center-point workflow snapshot from the original experiment suite.
Step 1 maps each cell center into the annulus, then Step 2-4 build topology
constraints and a topology-aware partition.

## Code

- `scripts/workflow/step01_build_center_polar_annulus.py`
- `scripts/workflow/step02_build_center_topology.py`
- `scripts/workflow/step03_build_line_md_scaffold.py`
- `scripts/workflow/step04_build_line_md_layer_column_partition.py`
- `scripts/workflow/step04_build_topology_aware_partition.py`
- `scripts/workflow/step04_build_topology_aware_partition_v2.py`

## Results

Primary result directory:

- `results/workflow/step04_topology_aware_v2/`

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

The metrics are stored in
`results/workflow/step04_topology_aware_v2/step04_topology_aware_metrics_v2.json`.
