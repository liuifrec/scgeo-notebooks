# GSE249479 — human HSPC inflammatory xenograft dataset

## Purpose and inference status

GSE249479, a public human HSPC inflammatory xenograft dataset, is used to
validate treatment-associated geometry and comparator complementarity. The
compact analysis retained 34,432 cells, but no recoverable biological-replicate
identity could be linked to cells. Every GSE249479 result is therefore
`descriptive_only`.

## Notebook order

1. `00_metadata_replication_and_memory_audit.ipynb`
2. `00a_download_and_reconstruct.ipynb` (provenance/reconstruction entry point;
   not required when the frozen workstation H5AD is present)
3. `01_sparse_qc_and_compact_object.ipynb`
4. `02_annotation_and_signatures.ipynb`
5. `03a`–`03d`: representation generation and quality
6. `04a`–`04c`: descriptive ScGeo treatment geometry
7. `05_augur_comparator.ipynb`
8. `06_scgeo_augur_comparison.ipynb`
9. `07_manuscript_figure_assembly.ipynb`

## Implementation, inputs, and outputs

Implementations are in `scripts/gse249479_*`, with staged executors for audit,
representations, ScGeo, and comparator work. The frozen source and compact H5AD
locations can be overridden with `SCGEO_GSE249479_H5AD` and
`SCGEO_GSE249479_COMPACT_H5AD`; generated artifacts use
`SCGEO_GSE249479_OUTPUT_DIR` where supported. Large H5ADs, scVI models,
executed notebooks, figures, tables, and logs remain ignored.

## Principal accepted findings

Frozen descriptive evidence includes TNF shifts in Activated HSC and Lymphoid,
broader LPS displacement, an abundance-dominant LPS HSC/quiescent pattern, high
rank agreement among the primary representations, and diffusion-map
sensitivity with limited influence on the primary conclusions. Official R
Augur 1.0.3 is the principal comparator. The Python method is an Augur-inspired
Python approximation used only for implementation sensitivity.

## Limitations

Cells, conditions, samples, condition batches, libraries, and SouporCell clades
are not biological replicates here. Signatures and marker-inferred labels do not
provide replication. PCA20/PCA30/PCA50 are nested views, not independent
confirmations, and UMAP is display-only.
