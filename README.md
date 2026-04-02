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
4. **Run notebooks in the execution order below** (data prep first, then manuscript notebooks).

> Notes:
> - Package versions may need to match the manuscript environment exactly for byte-level reproducibility.
> - If you already use a lab/cluster environment, keep your existing setup and use the notebook order in this README as the source of truth.

## Repository structure

- `notebooks/data_prep/` — **preprocessing pipeline** required before manuscript analyses.
- `notebooks/manuscript/` — **reproducible paper figures/results**.
- `notebooks/exploration/` — **side analyses** (helpful context, not required for core reproduction).
- `notebooks/tutorials/` — **user-facing demos** and learning-oriented walkthroughs.

## Execution order (reviewer-focused)

Run notebooks in this exact order.

### 1) Data preparation

1. `notebooks/data_prep/01_gse280305_paths.ipynb`
2. `notebooks/data_prep/02_gse280305_pathC_velocity.ipynb`  
   - corresponds to the manuscript velocity-prep step sometimes referenced as `02_gse280305_path_velocity.ipynb`
3. `notebooks/data_prep/03_gse280305_cellrank_sparse.ipynb`
4. `notebooks/data_prep/04_scgeo_gse280305_phase1_qc.ipynb`

### 2) Manuscript analyses (IMPORTANT ORDER)

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
- Exploration notebooks are treated as optional/non-required.
- Tutorial notebooks are demo-oriented and not required for manuscript reproduction.
- Existing numeric prefixes (for ordering) are kept where already meaningful.

## Data provenance (GSE280305)

- Dataset: **NCBI GEO accession GSE280305** (radiation-induced hematopoietic recovery).
- Biological context in these notebooks: mouse LSK cells across post-irradiation timepoints (D8, D11, D14, D21).
- Raw/intermediate data files are not stored in this repository.

Recommended provenance workflow:
1. Download source data from GEO (`GSE280305`) and record download date + GEO file checksums in your local run log.
2. Place data in your local/project data location expected by the data prep notebooks.
3. Execute `notebooks/data_prep/` notebooks in order before running manuscript notebooks.

## Scope summary

- **manuscript** = reproducible paper figures/results
- **exploration** = side analyses (not required)
- **tutorials** = user demos
- **data_prep** = preprocessing pipeline
