# Dataset map

| Dataset/workflow | Entry point | Main scripts | Analysis unit | Status |
|---|---|---|---|---|
| Synthetic benchmark | `notebooks/benchmarks/00_manuscript_benchmark_overview.ipynb` | `scripts/execute_revision_notebooks.py` | One frozen simulation job/seed | Controlled synthetic evidence only |
| Pancreas | `notebooks/public_validation/pancreas/00_acquire_and_validate.ipynb` | `scripts/execute_pancreas_validation.py`, `scripts/pancreas_validation_common.py` | Cells/states in a public developmental dataset | Descriptive dynamics validation |
| Dataset B, GSE249479 | `notebooks/public_validation/gse249479/00_metadata_replication_and_memory_audit.ipynb` | `scripts/gse249479_*`, staged executors | Cells and marker-inferred states; no valid biological sample key | `descriptive_only` |
| Dataset C, GSE211713 | `notebooks/public_validation/gse211713/00_study_design_replication_and_file_audit.ipynb` | `scripts/gse211713_*`, staged executors | Independent mouse/GSM libraries | Replicate-aware primary associations; cross-sectional |
| Original recovery workflow, GSE280305 | `notebooks/data_prep/01_gse280305_paths.ipynb` | notebook-local and exploration steps | See original design metadata | Manuscript workflow with semi-manual steps |

Dataset B uses official R Augur as its primary external comparator. Dataset C
does not use dose, time, irradiation status, condition, or mouse as a
batch-removal covariate. Neither Dataset B nor Dataset C uses Harmony or
Scanorama because no eligible independent technical batch was identified.
