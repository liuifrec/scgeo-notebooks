<p align="center">
  <img src="https://raw.githubusercontent.com/liuifrec/scgeo/main/docs/assets/ScGeo_logo.png" width="360" alt="ScGeo logo">
</p>

<h1 align="center">ScGeo workflows</h1>

<p align="center">Explore the analyses behind geometry-aware single-cell interpretation.</p>

[![Paper DOI](https://img.shields.io/badge/Paper-10.1093%2Fbioadv%2Fvbag266-blue)](https://doi.org/10.1093/bioadv/vbag266)
[![Python package](https://img.shields.io/pypi/v/scgeo)](https://pypi.org/project/scgeo/)

**Read the paper, follow the evidence, and adapt the workflows to your own study.**
This repository connects ScGeo's published analyses to their source notebooks,
configurations, execution scripts, and figure-generation steps. It is the
reproducibility companion to the [`scgeo` Python package](https://github.com/liuifrec/scgeo).

[ScGeo 0.2.0 tutorials](notebooks/tutorials/v0.2.0/README.md) ·
[Start here](docs/START_HERE.md) ·
[Workflow catalog](docs/DATASET_MAP.md) ·
[Results guide](docs/RESULTS_OVERVIEW.md) ·
[Plotting workflows](docs/PLOTTING_WORKFLOWS.md) ·
[Reproduce the paper](docs/REPRODUCIBILITY_GUIDE.md) ·
[Cite and reuse](docs/PUBLICATION_GUIDE.md)

## Find your starting point

| Your goal | Start with | What you get |
|---|---|---|
| Use ScGeo on your own AnnData | [ScGeo 0.2.0 tutorials](notebooks/tutorials/v0.2.0/README.md) and [package quick start](https://github.com/liuifrec/scgeo/blob/v0.2.0/docs/QUICKSTART.md) | Twelve independent API notebooks with actual figure and source-data exports |
| Understand the report format | [Small synthetic tutorial](notebooks/tutorials/01_quickstart_perturbation_report.ipynb) | A no-download conceptual example with a figure, source table, and metadata |
| Redraw stored evidence tables | [Released plotting workflows](docs/PLOTTING_WORKFLOWS.md) | Source-aware state, contrast, sample, and representation figures |
| Understand the published evidence | [Results guide](docs/RESULTS_OVERVIEW.md) | Main findings, their interpretation, and the workflow behind each claim |
| Reproduce a dataset or figure | [Workflow catalog](docs/DATASET_MAP.md) and [reproducibility guide](docs/REPRODUCIBILITY_GUIDE.md) | Inputs, environments, execution order, and frozen checkpoints |
| Cite or adapt an analysis | [Citation and reuse guide](docs/PUBLICATION_GUIDE.md) | Paper citation, reference-manager files, and a provenance checklist |

The new [0.2.0 learning track](notebooks/tutorials/v0.2.0/README.md) calls the published library directly and covers robust sample-aware analysis, representation comparisons, PAGA pies, matched velocity fields, annotation associations, figure composition, the atlas, gene/program response curves, sparse states, and the stored-evidence compass. Eleven lessons use synthetic inputs; the final lesson optionally uses checksum-matched public PBMC/pancreas files. No private data or manuscript downloads are needed.

The older conceptual tutorial remains available unchanged; it constructs its report with NumPy, pandas and Matplotlib rather than calling the library.

## Analyses included

ScGeo separates state displacement, sample uncertainty, representation
stability, local geometry, distributional change, and dynamics agreement.
The workflows below show how those evidence layers can be interpreted across
different biological designs.

| Workflow | Biological question or purpose | Scope of the evidence |
|---|---|---|
| [GSE280305: hematopoietic recovery](notebooks/manuscript/README.md) | How does hematopoietic state space change across post-irradiation time points? | Original biological case study; D8 is a post-irradiation reference |
| [Synthetic benchmarks](notebooks/benchmarks/README.md) | How do estimators and diagnostics behave under controlled perturbations? | Stored ground truth with separate calibration and held-out evaluation |
| [Pancreatic development](notebooks/public_validation/pancreas/README.md) | How does representation geometry agree with developmental dynamics? | Descriptive validation using related scVelo/CellRank context |
| [GSE249479: human HSPC inflammation](notebooks/public_validation/gse249479/README.md) | How do TNF- and LPS-associated shifts relate to abundance and separability? | 34,432 retained cells; `descriptive_only` because biological-replicate identity is unavailable |
| [GSE211713: mouse-lung radiation](notebooks/public_validation/gse211713/README.md) | Which early and late radiation associations are stable across representations? | 20 independent mouse libraries; 131,157 retained cells; replicate-aware, cross-sectional contrasts |

The [results guide](docs/RESULTS_OVERVIEW.md) preserves both supported findings
and limitations. Official R Augur is the principal GSE249479 comparator; the
retained Python approximation is an implementation-sensitivity analysis.

## Current use and paper reproduction

For learning the package and starting a new analysis, install the published
library in a fresh environment (Python 3.12 recommended):

```bash
python -m pip install "scgeo[scanpy,notebook]==0.2.0"
```

See the [package installation guide](https://github.com/liuifrec/scgeo/blob/main/docs/INSTALLATION.md)
for optional workflows and conda-managed environments. This repository is a
collection of workflows, not an additional Python package to install.

For the new tutorial track, use its [separate environment and runner](notebooks/tutorials/v0.2.0/README.md). The runner verifies source snapshots and figure/table hashes, and preserves output-free source notebooks.

For the paper's frozen evidence, use the recorded environments, numerical
inputs, and checkpoints:

| Component | Frozen checkpoint |
|---|---|
| ScGeo package | [`9a0ed16cbaa57f935f9c9bc87d1643a25b51012c`](https://github.com/liuifrec/scgeo/tree/9a0ed16cbaa57f935f9c9bc87d1643a25b51012c) |
| Companion merge checkpoint | [`958a0ada568c53e78864ede022c050830cb55a36`](https://github.com/liuifrec/scgeo-notebooks/tree/958a0ada568c53e78864ede022c050830cb55a36) |

A current package installation does not replace those recorded environments.
Some frozen executors retain workstation paths and named interpreters; the
[reproducibility guide](docs/REPRODUCIBILITY_GUIDE.md) identifies supported
path overrides and remaining portability limits before you run a workflow.

## What the workflows produce

Where supported by the workflow, generated bundles include:

- **Executed notebooks** showing the analysis in order.
- **PNG/SVG figures and source tables** for inspection and reuse.
- **Alt text, execution metadata, and checksums** linking presentation to its inputs.
- **Evidence ledgers** connecting a conclusion with its analysis unit and limitations.

Source notebooks are output-free. Generated copies, large data, models, and
figures are kept in ignored result or local-data directories. The repository
contains the instructions and source needed to generate these outputs; it is
not a complete download of all manuscript data and executed artifacts.

## Repository map

| Location | Contents |
|---|---|
| [`notebooks/tutorials/v0.2.0/`](notebooks/tutorials/v0.2.0/README.md) | Published API tutorials, with independent execution and exports |
| [`notebooks/tutorials/`](notebooks/tutorials/) | Earlier conceptual and application tutorials (unchanged) |
| [`notebooks/data_prep/`](notebooks/data_prep/), [`notebooks/exploration/`](notebooks/exploration/) | Original GSE280305 preprocessing and semi-manual reference annotation |
| [`notebooks/manuscript/`](notebooks/manuscript/README.md) | Original manuscript analyses and figure utilities |
| [`notebooks/benchmarks/`](notebooks/benchmarks/README.md) | Frozen synthetic benchmark review |
| [`notebooks/public_validation/`](notebooks/public_validation/) | Pancreas, GSE249479, and GSE211713 workflows |
| [`notebooks/revision_finalization/`](notebooks/revision_finalization/README.md) | Evidence ledgers, metric definitions, and manuscript audit trail |
| [`scripts/`](scripts/), [`configs/`](configs/), [`environment/`](environment/) | Reusable implementations, frozen settings, and environment specifications |

[Thin notebook architecture](docs/THIN_NOTEBOOK_ARCHITECTURE.md) explains how
source notebooks, reusable code, and generated outputs fit together.

## Cite ScGeo

Liu Y-C, Yoshida K. **ScGeo reveals non-canonical trajectories beyond RNA
velocity in radiation-induced hematopoietic recovery.** *Bioinformatics Advances*.
2026;6(1):vbag266. [doi:10.1093/bioadv/vbag266](https://doi.org/10.1093/bioadv/vbag266).

Use [BibTeX](CITATION.bib), [RIS](CITATION.ris), or GitHub's **Cite this repository**
menu, powered by [CITATION.cff](CITATION.cff). The [citation and reuse guide](docs/PUBLICATION_GUIDE.md)
explains how to report software versions, workflow commits, datasets, and
adaptations alongside the paper citation.

## Interpretation and support

Independent biological sample IDs determine the scope of inference. Related
PCA dimensions share a basis, UMAP is used for display in the quantitative
validation workflows, and agreement with a dynamics model inherits that
model's assumptions. GSE249479 remains descriptive; GSE211713 comparisons are
between different mice and do not establish longitudinal change or causality.

For a workflow question or reproducibility issue, open an
[issue](https://github.com/liuifrec/scgeo-notebooks/issues) with the notebook,
repository commits, environment, and relevant error or output. Package API
questions belong in the [ScGeo package issue tracker](https://github.com/liuifrec/scgeo/issues).
