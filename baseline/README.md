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

| Step | Method | Code | Results |
|---|---|---|---|
| 01 | Direct four-corner polar mapping | `code/step01_direct_polar_mapping.py` | `results/step01_direct_polar_mapping/` |
| 02 | Uniform harmonic annulus | `code/step02_uniform_harmonic_annulus.py` | `results/step02_uniform_harmonic_annulus/` |
| 03 | Cotangent/conformal annulus | `code/step03_cotangent_conformal_annulus.py` | `results/step03_cotangent_conformal_annulus/` |
| 04 | Area-preserving optimization | `code/step04_area_preserving_annulus.py` | `results/step04_area_preserving_annulus/` |

`code/baseline_common.py` is shared by Steps 02-04.
