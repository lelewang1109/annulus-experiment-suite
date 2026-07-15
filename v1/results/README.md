# v1 results

This folder contains the Step 3 topology constraints needed to rerun the v1 topology-aware optimizer.

A full v1 output snapshot is not committed in the current repository. The standardized committed result snapshot is v2, and v3 is the 0715 four-vertex supplement.

To generate v1 output, run from `v1/`:

```bash
python3 scripts/workflow/step04_build_topology_aware_partition.py
```

The default generated output directory is:

- `results/workflow/step04_topology_aware_v1/`
