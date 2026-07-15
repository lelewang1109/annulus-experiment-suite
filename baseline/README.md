# Baseline experiments

This folder contains the four baseline experiments. They test continuous or
geometry-driven annulus mappings before the discrete topology-aware workflow.

## Structure

```text
baseline/
├── data/
├── code/
└── results/
```

## Code and results

Baseline methods are independent experiments, so they are named by method rather
than by step.

| Method | Code | Results |
|---|---|---|
| Direct four-corner polar mapping | `code/direct_polar_mapping.py` | `results/direct_polar_mapping/` |
| Uniform harmonic annulus | `code/uniform_harmonic_annulus.py` | `results/uniform_harmonic_annulus/` |
| Cotangent/conformal annulus | `code/cotangent_conformal_annulus.py` | `results/cotangent_conformal_annulus/` |
| Area-preserving optimization | `code/area_preserving_annulus.py` | `results/area_preserving_annulus/` |

`code/baseline_common.py` is shared by the harmonic, conformal, and
area-preserving experiments.
