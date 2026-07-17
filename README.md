# scgeo-notebooks

Reproducibility companion notebooks for the **ScGeo manuscript** (radiation-induced hematopoietic recovery; **GSE280305**).

This repository is organized so reviewers can run a minimal, manuscript-focused pipeline end-to-end, while still keeping side analyses and demos available.

## Quick start

1. **Create and activate a Python environment** (example with `conda`):
   ```bash
   conda create -n scgeo python=3.10 -y
   conda activate scgeo
   ```
2. **Install core notebook dependencies** used across these analyses:
   ```bash
   pip install jupyterlab scanpy scvelo cellrank anndata scanorama pandas numpy matplotlib seaborn
   ```
3. **Launch Jupyter** from the repo root:
   ```bash
   jupyter lab
   ```
4. **Run notebooks in the execution order below** (data prep, then reference prep/annotation, then manuscript notebooks).

> Notes:
> - Package versions may need to match the manuscript environment exactly for byte-level reproducibility.
> - If you already use a lab/cluster environment, keep your existing setup and use the notebook order in this README as the source of truth.

## Repository structure

- `notebooks/data_prep/` — **preprocessing and dataset preparation** required before downstream analyses.
- `notebooks/exploration/` — **intermediate, analysis-relevant notebooks**, including manual reference preparation and manual annotation steps used upstream of final manuscript outputs.
- `notebooks/manuscript/` — **final figure/result generation** notebooks for manuscript outputs.
- `notebooks/benchmarks/` — **synthetic revision benchmark notebooks** that read frozen manuscript-profile simulation outputs.
- `notebooks/public_validation/` — **public validation notebooks** using external public datasets only.
- `notebooks/tutorials/` — **user-facing demos** and learning-oriented walkthroughs.

## Synthetic revision benchmark workflow

The revision benchmark notebooks use only the frozen synthetic manuscript-profile outputs. They do not rerun the manuscript simulation suite and do not tune or propose new thresholds.

Configure the frozen benchmark source with:

```bash
export SCGEO_BENCHMARK_DIR=/path/to/manuscript_v1
```

If the variable is unset, the notebooks use the relative default recorded in `configs/manuscript_benchmark_v1.json`:

```text
../scgeo/results/simulation/manuscript_v1
```

The pinned revision config records the `manuscript_v1` protocol, source commit, expected job counts, calibration seeds, held-out evaluation seeds, required audit files, and checksum manifests. The committed manifests are:

- `results_manifest/benchmark_files.csv`
- `results_manifest/checksums.sha256`

Run the clean-kernel execution test from the repository root:

```bash
python scripts/execute_revision_notebooks.py
```

The runner executes these notebooks in fresh kernels:

1. `notebooks/benchmarks/00_manuscript_benchmark_overview.ipynb`
2. `notebooks/benchmarks/01_robust_estimator_comparison.ipynb`
3. `notebooks/benchmarks/02_representation_and_dynamics_validation.ipynb`
4. `notebooks/benchmarks/03_framework_ablation.ipynb`
5. `notebooks/benchmarks/04_synthetic_geometry_stress_test.ipynb`
6. `notebooks/benchmarks/05_synthetic_dynamics_stress_test.ipynb`
7. `notebooks/tutorials/01_quickstart_perturbation_report.ipynb`

Generated revision artifacts are written under `results/revision_synthetic_benchmark/`:

- `figures/` contains SVG and PNG figures.
- `figure_sources/` contains figure-source CSV tables.
- `alt_text/` contains deterministic alt text.
- `metadata/` records Python, scgeo, package, repository commit, source commit, protocol, and timestamp metadata.
- `execution/revision_notebook_execution_report.json` records clean-kernel runtimes and artifact inventory.

The benchmark notebooks keep calibration seeds (`0-4`) and held-out evaluation seeds (`5-19`) visibly separated. One simulation job/seed is treated as the independent unit. The notebooks explicitly report the negative revision findings captured by the frozen audit outputs: balanced-replicate seed dependence, incomplete state-level instability recall, zero local-distortion recall in representation-corruption jobs, and imperfect bootstrap uncertainty coverage.

## Public pancreas Dataset D validation

Dataset D is a public endocrine-pancreas developmental-dynamics validation using the official CellRank pancreas dataset. It does not modify the ScGeo package, does not change frozen synthetic thresholds or protocol settings, and does not create artificial treatment/control labels.

Create the optional environment from the repository root:

```bash
conda env create -f environment/pancreas_environment.yml
conda activate scgeo-pancreas-dataset-d
```

The public data and output locations are configurable:

```bash
export SCGEO_PANCREAS_DATA_DIR=data/public/pancreas_dataset_d
export SCGEO_PANCREAS_OUTPUT_DIR=results/public_validation/pancreas_dataset_d
```

Run the clean-kernel public validation workflow:

```bash
python scripts/execute_pancreas_validation.py
```

The runner executes these source-output-free notebooks and saves executed review copies under `results/public_validation/pancreas_dataset_d/executed_notebooks/`:

1. `notebooks/public_validation/pancreas/00_acquire_and_validate.ipynb`
2. `notebooks/public_validation/pancreas/01_scvelo_dynamical_velocity.ipynb`
3. `notebooks/public_validation/pancreas/02_cellrank_fates.ipynb`
4. `notebooks/public_validation/pancreas/03_scgeo_representation_dynamics.ipynb`
5. `notebooks/public_validation/pancreas/04_manuscript_figures.ipynb`

The workflow records official public dataset checksums, scVelo dynamical velocity outputs, a CellRank VelocityKernel/GPCCA comparator, ScGeo-style representation-dynamics evidence across PCA20, PCA30, PCA50, diffusion map, and UMAP diagnostics, negative controls, PNG/SVG figures, figure-source CSVs, deterministic alt text, metadata, and version records.

CellRank output is treated as a complementary probabilistic comparator derived from scVelo RNA velocity. It is not reported as evidence independent of RNA velocity.

## Workflow order (manuscript-focused)

For manuscript-oriented reruns, use the notebook groups in this order:

1. **Data preparation**: run `notebooks/data_prep/` notebooks to preprocess and prepare inputs.
2. **Reference preparation and annotation (semi-manual)**:
   - `notebooks/exploration/06_Ref_prep.ipynb` — manual reference preparation.
   - `notebooks/exploration/07_Reference_based_annotation.ipynb` — manual reference-based annotation.
3. **Final manuscript outputs**: run `notebooks/manuscript/` notebooks for final figures/graphs and reported result outputs.

> Reproducibility note: core preprocessing and manuscript notebooks are scripted, while the two exploration notebooks above include manual/semi-manual decisions that are part of the analysis flow.

## External references for annotation

Some notebooks (particularly in `notebooks/exploration/`) use external reference resources to assist with cell type annotation and interpretation:

- Hematopoietic reference atlas (HemAtlas)
- Azimuth reference mapping (HuBMAP): https://azimuth.hubmapconsortium.org/

These resources were used for manual or semi-guided annotation steps and are not part of the ScGeo framework itself. They serve as biological references to support interpretation of embedding structure and inferred trajectories.

## Execution order (reviewer-focused)

Run notebooks in this exact order.

### 1) Data preparation

1. `notebooks/data_prep/01_gse280305_paths.ipynb`
2. `notebooks/data_prep/02_gse280305_pathC_velocity.ipynb`  
   - corresponds to the manuscript velocity-prep step sometimes referenced as `02_gse280305_path_velocity.ipynb`
3. `notebooks/data_prep/03_gse280305_cellrank_sparse.ipynb`
4. `notebooks/data_prep/04_scgeo_gse280305_phase1_qc.ipynb`

### 2) Reference preparation and annotation (semi-manual, upstream)

1. `notebooks/exploration/06_Ref_prep.ipynb`
2. `notebooks/exploration/07_Reference_based_annotation.ipynb`

### 3) Manuscript analyses (IMPORTANT ORDER)

1. `notebooks/manuscript/05_OOD.ipynb` (**Figure 2**)
2. `notebooks/manuscript/Velocity_shift_alignment.ipynb` (**Figure 3**)
3. `notebooks/manuscript/Final_summary.ipynb` (**intermediate integrated summary**)
4. `notebooks/manuscript/test_driver_genes.ipynb` (**Figure 4**)

> Note: `notebooks/manuscript/Final_summary.ipynb` is intentionally run **before** `notebooks/manuscript/test_driver_genes.ipynb` despite its filename.

## Figure-to-notebook mapping

- **Figure 2** → `notebooks/manuscript/05_OOD.ipynb`
- **Figure 3** → `notebooks/manuscript/Velocity_shift_alignment.ipynb`
- **Figure 4** → `notebooks/manuscript/test_driver_genes.ipynb`
- **Intermediate integrated summary** → `notebooks/manuscript/Final_summary.ipynb`

## Notebook naming notes

To avoid breaking existing links/references, notebook filenames are preserved in this repository.

- Manuscript notebooks are documented above with explicit figure mapping for readability.
- Exploration notebooks include intermediate analysis steps; specifically, `06_Ref_prep.ipynb` and `07_Reference_based_annotation.ipynb` are upstream of final manuscript outputs.
- Tutorial notebooks are demo-oriented and not required for manuscript reproduction.
- Existing numeric prefixes (for ordering) are kept where already meaningful.

## Data provenance (GSE280305)

- Dataset: **NCBI GEO accession GSE280305** (radiation-induced hematopoietic recovery).
- Biological context in these notebooks: mouse LSK cells across post-irradiation timepoints (D8, D11, D14, D21).
- Raw/intermediate data files are not stored in this repository.

Recommended provenance workflow:
1. Download source data from GEO (`GSE280305`) and record download date + GEO file checksums in your local run log.
2. Place data in your local/project data location expected by the data prep notebooks.
3. Execute `notebooks/data_prep/`, then `notebooks/exploration/06_Ref_prep.ipynb` and `notebooks/exploration/07_Reference_based_annotation.ipynb`, before running manuscript notebooks.

## Scope summary

- **manuscript** = reproducible paper figures/results
- **exploration** = intermediate analysis notebooks (including manual reference prep/annotation)
- **tutorials** = user demos
- **data_prep** = preprocessing pipeline
