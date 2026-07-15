# v1 center-point topology-aware workflow

v1 preserves the first topology-aware optimizer for the annular partition
workflow.

## Code

- `scripts/workflow/step04_build_topology_aware_partition.py`
- `scripts/workflow/annulus_topology_utils.py`

The script reads the canonical center-polar cells from
`data/input/beijing_grid_cells.csv` and the topology constraints expected under
`results/workflow/step03/step03_topology_constraints.csv` when run from a full
workflow checkout.

## Result status

This folder records the v1 code and canonical input files, but does not include a
committed full v1 result snapshot. The public result snapshot in this repository
is v2, and v3 is added as a separate four-vertex experiment.

To regenerate v1 results from this package, provide or copy the Step 3 topology
constraints into `results/workflow/step03/`, then run:

```bash
python3 scripts/workflow/step04_build_topology_aware_partition.py
```

The default output directory is `results/workflow/step04_topology_aware_v1/`.
