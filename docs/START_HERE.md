# Start here

This companion helps you move from the ScGeo paper to the workflow behind a
result. Choose a path below according to what you want to do.

[Repository home](../README.md) · [Workflow catalog](DATASET_MAP.md) ·
[Results guide](RESULTS_OVERVIEW.md) · [Cite and reuse](PUBLICATION_GUIDE.md)

## Use ScGeo on your own data

Start with the package [installation guide](https://github.com/liuifrec/scgeo/blob/main/docs/INSTALLATION.md)
and [quick start](https://github.com/liuifrec/scgeo/blob/main/docs/QUICKSTART.md).
They introduce the AnnData inputs, sample identifiers, analysis functions, and
result tables. You can use the package without downloading the manuscript
analysis data or cloning this companion.

For a new study, choose the biological contrast, sample unit, state labels,
and quantitative representations before interpreting the output. The
[analysis-layer guide](https://github.com/liuifrec/scgeo/blob/main/docs/ANALYSIS_LAYERS.md)
explains what each evidence layer measures.

## Try a small report example

The [synthetic report tutorial](../notebooks/tutorials/01_quickstart_perturbation_report.ipynb)
creates its input at runtime. It demonstrates descriptive state-shift ranking,
sample-level resampling, and export of a figure, source CSV, alt text, and
metadata. Its calculations use NumPy and pandas; it is a conceptual report
example rather than an API tutorial or a reproduction of the paper.

In a fresh Python 3.12 environment:

```bash
git clone https://github.com/liuifrec/scgeo-notebooks.git
cd scgeo-notebooks
python -m pip install "scgeo[scanpy,notebook]==0.1.0"
python -m ipykernel install --user --name scgeo-tutorial --display-name "ScGeo tutorial"
python scripts/execute_revision_notebooks.py --notebook notebooks/tutorials/01_quickstart_perturbation_report.ipynb --kernel scgeo-tutorial --output-dir results/tutorial_quickstart --keep-executed
```

The selected tutorial needs no external dataset or frozen benchmark directory.
The runner writes outputs beneath `results/tutorial_quickstart/`, including an
executed copy under `executed_notebooks/`. The tracked notebook stays
output-free. See [notebook architecture](THIN_NOTEBOOK_ARCHITECTURE.md) for the
source/output convention.

## Redraw existing result tables

For the development plotting API, see [plotting workflows](PLOTTING_WORKFLOWS.md).
It includes an executable renderer for existing CSVs and a link to the package's
synthetic visual gallery. These new interfaces require the matching development
checkout; they are not in the published `scgeo==0.1.0` distribution.

## Read the published evidence

Begin with the [results guide](RESULTS_OVERVIEW.md), then use the linked dataset
README to find its ordered notebooks and underlying implementation. Read each
result with its analysis unit and inference scope:

| Evidence | Interpretation |
|---|---|
| Synthetic benchmark | Controlled method behavior under specified simulations |
| GSE280305 | Original post-irradiation case study, with D8 as a post-exposure reference |
| Public pancreatic development | Descriptive comparison with developmental dynamics |
| GSE249479 human HSPC inflammation | Descriptive treatment geometry; valid biological-replicate identity is unavailable |
| GSE211713 mouse-lung radiation | Replicate-aware, cross-sectional associations across independent mice |

The evidence-ledger and figure-source locations are described in the results
guide. Large generated artifacts are usually ignored by Git, so a source
notebook may refer to files that must first be generated or recovered from the
recorded analysis environment.

## Reproduce a dataset or rebuild a figure

1. Choose a workflow from the [dataset map](DATASET_MAP.md).
2. Read its folder README and identify the required raw/processed inputs.
3. Follow the [reproducibility guide](REPRODUCIBILITY_GUIDE.md) for checkpoints,
   environments, checksums, and path configuration.
4. Run the numbered stages through the matching executor; inspect each
   stage's metadata before moving to the next.
5. For figure assembly, use the recorded numerical inputs and verify their
   checksums before regenerating presentation outputs.

Full manuscript workflows require the recorded data, package environment, and
sometimes workstation paths. A current pip installation is a learning and
new-analysis route; it does not by itself reproduce the frozen manuscript.

## Cite or adapt a workflow

Use the [citation and reuse guide](PUBLICATION_GUIDE.md) for the paper citation,
reference-manager downloads, and the information to record in your Methods.
When adapting a workflow, identify changed parameters or preprocessing steps
and describe the resulting analysis as your adaptation.
