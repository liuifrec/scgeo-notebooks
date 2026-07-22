# Synthetic revision benchmarks

## Purpose and inference status

These notebooks review checksum-pinned synthetic manuscript-profile outputs.
They test controlled method behavior with held-out evaluation and make no
biological claim. One simulation job/seed is the independent evaluation unit.

## Notebook order

1. `00_manuscript_benchmark_overview.ipynb`
2. `01_robust_estimator_comparison.ipynb`
3. `02_representation_and_dynamics_validation.ipynb`
4. `03_framework_ablation.ipynb`
5. `04_synthetic_geometry_stress_test.ipynb`
6. `05_synthetic_dynamics_stress_test.ipynb`

## Implementation, inputs, and outputs

Run through `scripts/execute_revision_notebooks.py`. Inputs are the frozen
benchmark directory selected by `SCGEO_BENCHMARK_DIR` and
`configs/manuscript_benchmark_v1.json`. Generated review copies, figures,
source CSVs, alt text, metadata, and execution logs are written under ignored
`results/revision_synthetic_benchmark/` paths.

## Frozen findings and limitations

The evidence includes robust-estimator, representation, local-geometry,
dynamics, and framework-ablation results, including negative findings for seed
dependence, instability recall, local-distortion recall, and uncertainty
coverage. Controlled corruption sensitivity is not general out-of-distribution
detection, and synthetic performance does not establish biological validity.
