# Step 04 topology-aware v1 results

The full v1 output snapshot is not committed in this project.

To regenerate it from `workflow-v1/`, run:

```bash
python3 code/step04_build_topology_aware_v1.py \
  --input-csv data/input/beijing_grid_cells.csv \
  --topology-csv results/step03_line_md_scaffold/step03_topology_constraints.csv \
  --output-dir results/step04_topology_aware_v1
```
