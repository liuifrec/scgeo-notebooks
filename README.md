# scgeo-notebooks

This repository contains analysis workflows and figure-generation notebooks for the ScGeo framework.

## Study

Radiation-induced hematopoietic recovery (GSE280305)

- Mouse LSK cells
- Timepoints: D8, D11, D14, D21
- Integration: Scanorama
- Dynamics: scVelo + CellRank
- Geometry: ScGeo

## Key findings

- Perturbation induces structured displacement in embedding space
- Early states occupy out-of-distribution regions
- Recovery is only partially aligned with RNA velocity
- Discordant states reveal stress-response and adaptive programs

## Notebook map

### Data preparation
- `data_prep/01–04`: preprocessing, velocity, integration

### Manuscript analyses
- `manuscript/05_OOD.ipynb` → Figure 2 (OOD + redistribution)
- `manuscript/Velocity_shift_alignment.ipynb` → Figure 3
- `manuscript/test_driver_genes.ipynb` → Figure 4
- `manuscript/Final_summary.ipynb` → integrated results

### Tutorials
- PBMC demo notebooks

## Data

Raw and intermediate data are not included.

To reproduce:
1. Download GSE280305
2. Run notebooks in `notebooks/data_prep/`