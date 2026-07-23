# Dataset map

| Dataset/workflow | Entry point | Main scripts | Analysis unit | Status |
|---|---|---|---|---|
| Synthetic benchmark | `notebooks/benchmarks/00_manuscript_benchmark_overview.ipynb` | `scripts/execute_revision_notebooks.py` | One frozen simulation job/seed | Controlled synthetic evidence only |
| Public pancreas | `notebooks/public_validation/pancreas/00_acquire_and_validate.ipynb` | `scripts/execute_pancreas_validation.py`, `scripts/pancreas_validation_common.py` | Cells/states in the public pancreas developmental-dynamics dataset | Descriptive dynamics validation |
| GSE249479 — human HSPC inflammatory xenograft dataset | `notebooks/public_validation/gse249479/00_metadata_replication_and_memory_audit.ipynb` | `scripts/gse249479_*`, staged executors | Cells and marker-inferred states; no valid biological sample key | `descriptive_only` |
| GSE211713 — mouse-lung radiation dataset | `notebooks/public_validation/gse211713/00_study_design_replication_and_file_audit.ipynb` | `scripts/gse211713_*`, staged executors | Independent mouse/GSM libraries | Replicate-aware primary associations; cross-sectional |
| Original GSE280305 radiation-recovery case study | `notebooks/data_prep/01_gse280305_paths.ipynb` | notebook-local and exploration steps | See original design metadata | Manuscript workflow with semi-manual steps |

GSE249479 uses official R Augur as its principal external comparator. GSE211713
does not use dose, time, irradiation status, condition, or mouse as a
batch-removal covariate. Neither GSE249479 nor GSE211713 uses Harmony or
Scanorama because no eligible independent technical batch was identified.

Internal development labels such as “Dataset B” and “Dataset C” are not used in
reviewer-facing documentation; public datasets are identified by accession and
biological system.
