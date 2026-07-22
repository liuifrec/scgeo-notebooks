# scgeo-notebooks

Reproducibility companion for the
[`scgeo`](https://github.com/liuifrec/scgeo) Python package and the ScGeo
manuscript revision.

## 1. Relationship to the ScGeo package

The package repository contains the installable library, tests, and API
contracts. This repository contains source notebooks, execution wrappers,
public-data provenance, figure assembly, evidence ledgers, and reviewer-facing
audits. Revision workflows pin the package at
[`9a0ed16`](https://github.com/liuifrec/scgeo/tree/9a0ed16cbaa57f935f9c9bc87d1643a25b51012c).

Start with [Start here](docs/START_HERE.md). Package users who do not need to
reproduce the manuscript can instead follow the package
[quick start](https://github.com/liuifrec/scgeo/blob/main/docs/QUICKSTART.md).

## 2. Repository map

| Path | Role |
|---|---|
| `notebooks/benchmarks/` | Frozen synthetic revision benchmark review |
| `notebooks/public_validation/` | Public pancreatic-development, GSE249479 human HSPC inflammatory xenograft, and GSE211713 mouse-lung radiation validation |
| `notebooks/data_prep/` | Original GSE280305 preprocessing |
| `notebooks/exploration/` | Semi-manual reference preparation and annotation |
| `notebooks/manuscript/` | Original manuscript-oriented analyses and figures |
| `notebooks/revision_finalization/` | Reviewer metric, aggregation, and terminology audits |
| `scripts/` | Testable implementations and clean-kernel runners |
| `configs/` | Frozen workflow settings and checksums |
| `environment/` | Workflow-specific environment specifications |
| `results/` | Mostly ignored generated artifacts; a small set of tracked manifests may remain |

See [Dataset map](docs/DATASET_MAP.md) and [Results overview](docs/RESULTS_OVERVIEW.md).

## 3. Quick-start decision tree

1. **Learning the package API?** Use the package quick start and tutorials.
2. **Reviewing revision evidence?** Read the frozen evidence summaries before
   executing anything.
3. **Reproducing a public dataset?** Open that dataset’s folder README and run
   its numbered source notebooks through the corresponding executor.
4. **Rebuilding manuscript presentation only?** Use the assembly notebook after
   verifying its checksum-pinned numerical inputs.
5. **Reproducing the original GSE280305 workflow?** Follow the data-preparation,
   semi-manual annotation, and manuscript order in
   [`notebooks/manuscript/README.md`](notebooks/manuscript/README.md).

Do not substitute conditions, libraries, cells, or artificial partitions for a
biological sample identity.

## 4. Dataset table

| Workflow | Biological system | Inference status |
|---|---|---|
| Synthetic benchmark | Controlled simulations with stored truth and held-out evaluation | No biological claim |
| GSE280305 | Irradiated LSK hematopoietic recovery | Original biological case study |
| Public pancreas | Pancreatic developmental dynamics with scVelo/CellRank context | Descriptive dynamics validation |
| GSE249479 | Human HSPC inflammatory xenograft perturbation; 34,432 retained cells | `descriptive_only`; no recoverable biological-replicate identity |
| GSE211713 | Mouse-lung radiation response; 20 independent mouse libraries and 131,157 retained cells | Replicate-aware association; cross-sectional, not longitudinal |

Official R Augur is the principal GSE249479 comparator. The Python implementation
is retained only as an **Augur-inspired Python approximation** sensitivity
analysis.

## 5. Execution environments

Environment specifications are under `environment/`. Execute from the
repository root so relative paths and imports resolve consistently. The frozen
workflows used named kernels or explicit Python executables recorded in their
execution metadata; byte-level reproduction requires the corresponding package
versions as well as the same numerical inputs.

Common commands include:

```bash
python scripts/execute_revision_notebooks.py
python scripts/execute_pancreas_validation.py
python scripts/execute_gse249479_validation.py
python scripts/execute_gse211713_validation.py
```

Each dataset README identifies the appropriate staged executors and inputs.

## 6. Thin-notebook architecture

Most revision notebooks are deliberately short, output-free orchestration
entry points. Their corresponding modules under `scripts/` hold the reusable,
testable implementation. This keeps notebook diffs readable and permits
clean-kernel execution without embedding generated outputs in Git.

See [Thin notebook architecture](docs/THIN_NOTEBOOK_ARCHITECTURE.md).

## 7. Generated-versus-tracked artifact policy

Tracked source notebooks must have zero outputs and null execution counts.
Executed review copies, figures, source CSVs, alt text, logs, H5AD objects,
models, downloaded data, and caches are written under ignored result or local
data directories. Checksum manifests and source configuration may be tracked;
large generated artifacts are not.

See [Reproducibility guide](docs/REPRODUCIBILITY_GUIDE.md).

## 8. Frozen revision checkpoints

- ScGeo package: `9a0ed16cbaa57f935f9c9bc87d1643a25b51012c`
- Companion merge checkpoint: `958a0ada568c53e78864ede022c050830cb55a36`
- Dataset-specific numerical inputs and outputs: pinned in execution metadata,
  validation reports, or evidence ledgers rather than inferred from filenames.

Frozen thresholds and numerical results must not be changed during a
presentation or documentation pass.

## 9. Manuscript-oriented workflow

The original manuscript pipeline and the major-revision evidence package are
separate but related:

1. original GSE280305 preparation and semi-manual annotation;
2. original manuscript analyses;
3. frozen synthetic and public-data revision validation;
4. GSE249479 and GSE211713 manuscript assembly;
5. reviewer evidence ledgers and reproducibility audit.

The [Results overview](docs/RESULTS_OVERVIEW.md) summarizes accepted and negative
findings without recomputing them.

## 10. Limitations and semi-manual steps

- Original reference preparation and annotation include semi-manual decisions.
- GSE249479 lacks a recoverable biological-replicate identity and remains
  descriptive.
- GSE211713 has independent mouse libraries but is cross-sectional; it does not
  establish causality, persistence, reversal, or within-mouse change.
- Nested PCA dimensions share a basis and are not independent confirmations.
- UMAP is display-only in quantitative validation workflows.
- Some frozen validation scripts retain absolute workstation paths. Supported
  environment-variable overrides and remaining fixed paths are documented in
  the [Reproducibility guide](docs/REPRODUCIBILITY_GUIDE.md).
