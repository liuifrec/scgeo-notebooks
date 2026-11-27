# Workflow catalog and dataset map

Find the workflow that matches the biological question or methodological
component you want to inspect. Each linked folder README gives the notebook
order, inputs, implementation, and inference scope.

[Start here](START_HERE.md) · [Results guide](RESULTS_OVERVIEW.md) ·
[Reproducibility](REPRODUCIBILITY_GUIDE.md) · [Cite and reuse](PUBLICATION_GUIDE.md)

## Workflow entry points

| Workflow | Start here | What to inspect | Analysis unit and scope |
|---|---|---|---|
| ScGeo 0.2.0 learning track | [Tutorial index](../notebooks/tutorials/v0.2.0/README.md) | Released analysis and plotting APIs, including sparse coverage and molecular explanations | Independent synthetic samples; optional public descriptive rendering; not manuscript evidence |
| Small synthetic report | [Tutorial notebook](../notebooks/tutorials/01_quickstart_perturbation_report.ipynb) | Descriptive shift ranking and a figure/source-table bundle | Generated samples; conceptual example, not manuscript evidence |
| GSE280305 hematopoietic recovery | [Original manuscript guide](../notebooks/manuscript/README.md); [first preparation notebook](../notebooks/data_prep/01_gse280305_paths.ipynb) | Preparation, reference annotation, displacement, and geometry–velocity alignment | Follow original design metadata; comparisons are among post-irradiation time points |
| Synthetic benchmarks | [Benchmark guide](../notebooks/benchmarks/README.md); [overview notebook](../notebooks/benchmarks/00_manuscript_benchmark_overview.ipynb) | Estimator robustness, representation/dynamics diagnostics, and framework ablation | Frozen simulation job/seed; controlled method evidence |
| GSE132188-derived pancreatic development | [Pancreas guide](../notebooks/public_validation/pancreas/README.md); [acquisition notebook](../notebooks/public_validation/pancreas/00_acquire_and_validate.ipynb) | Geometry–dynamics relationships in scVelo/CellRank context | Cells/states; descriptive dynamics validation |
| GSE249479 human HSPC inflammatory xenograft | [GSE249479 guide](../notebooks/public_validation/gse249479/README.md); [design audit](../notebooks/public_validation/gse249479/00_metadata_replication_and_memory_audit.ipynb) | TNF/LPS geometry, abundance, representation sensitivity, and official R Augur comparison | No valid biological sample key; `descriptive_only` |
| GSE211713 mouse-lung radiation | [GSE211713 guide](../notebooks/public_validation/gse211713/README.md); [study-design audit](../notebooks/public_validation/gse211713/00_study_design_replication_and_file_audit.ipynb) | Mouse-level shifts, early/late contrasts, representation consensus, and fibroblast sensitivity | Independent mouse/GSM libraries; replicate-aware, cross-sectional associations |
| Manuscript evidence assembly | [Finalization guide](../notebooks/revision_finalization/README.md) | Definitions, evidence ledgers, figure inventory, and reproducibility audit | Presentation and evidence assembly from frozen numerical inputs |

## Implementation and environments

| Workflow | Main implementation or executor | Environment/configuration entry point |
|---|---|---|
| ScGeo 0.2.0 tutorials | [`execute_scgeo_v020_tutorials.py`](../scripts/execute_scgeo_v020_tutorials.py), exact release example snapshots | [`scgeo_v020_tutorials.txt`](../environment/scgeo_v020_tutorials.txt) (separate new environment) |
| Synthetic benchmark review | [`execute_revision_notebooks.py`](../scripts/execute_revision_notebooks.py) | [`manuscript_benchmark_v1.json`](../configs/manuscript_benchmark_v1.json) and the frozen package environment |
| Pancreas | [`execute_pancreas_validation.py`](../scripts/execute_pancreas_validation.py), [`pancreas_validation_common.py`](../scripts/pancreas_validation_common.py) | [`pancreas_environment.yml`](../environment/pancreas_environment.yml) |
| GSE249479 | [`execute_gse249479_validation.py`](../scripts/execute_gse249479_validation.py) and staged `gse249479_*` modules | [`gse249479_environment.yml`](../environment/gse249479_environment.yml) and dataset-specific configurations |
| GSE211713 | [`execute_gse211713_validation.py`](../scripts/execute_gse211713_validation.py) and staged `gse211713_*` modules | [`gse211713_environment.yml`](../environment/gse211713_environment.yml) and dataset-specific configurations |
| Original GSE280305 | Data-preparation, exploration, and manuscript notebooks | Recorded notebook metadata and semi-manual annotation steps |

Environment specifications describe workflow components; the full frozen run
also depends on the recorded package checkout, interpreter, inputs, and
versions. Read the [reproducibility guide](REPRODUCIBILITY_GUIDE.md) before
launching a full-data executor.

## Data provenance and comparison design

Public accessions: [GSE280305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280305),
[GSE132188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188),
[GSE249479](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249479), and
[GSE211713](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211713).
The pancreas workflow uses the processed public CellRank input recorded by its
acquisition stage rather than assuming that a GEO download alone recreates
that object.

Official R Augur is the principal GSE249479 comparator. Its Python approximation
is retained as implementation sensitivity. Neither GSE249479 nor GSE211713 uses
Harmony or Scanorama because no eligible independent technical batch was
identified. In GSE211713, dose, time, irradiation status, condition, and mouse
are biological design variables and are not used as batch-removal covariates.

Historical filenames retain labels such as `dataset_b`, `dataset_c`, and
`dataset_d` for compatibility. Public-facing descriptions identify each dataset
by accession and biological system; these labels do not define extra datasets.
