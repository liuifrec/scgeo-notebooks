# GSE211713 — mouse-lung radiation dataset

## Purpose and inference status

GSE211713, a public mouse-lung radiation dataset with 20 independent mouse
libraries, contains 131,157 retained cells. The three frozen primary contrasts
are replicate-aware associations based on mouse-level state centers. The design
is cross-sectional, not longitudinal, and does not establish causality,
persistence, reversal, or within-mouse change.

## Notebook order

1. `00_study_design_replication_and_file_audit.ipynb`
2. `01_download_and_sparse_reconstruction.ipynb`
3. `02_sparse_qc_and_compact_object.ipynb`
4. `03_coarse_annotation_and_fibroblast_states.ipynb`
5. `04a`–`04d`: representation ensemble and quality
6. `05a`–`05e`: replicate-aware ScGeo and cross-contrast summaries
7. `06_manuscript_figure_assembly.ipynb`

## Implementation, inputs, and outputs

Implementations are in `scripts/gse211713_*`, with staged validation,
representation, ScGeo, and assembly executors. Supported data and canonical-v2
path overrides are listed in `docs/REPRODUCIBILITY_GUIDE.md`. Generated MEX
downloads, H5ADs, scVI models, executed notebooks, numerical outputs, figures,
and logs remain under ignored local data or result paths.

## Principal accepted findings

Under the frozen primary consensus, no early-versus-control major state met the
stable-effect rule. Late-versus-control associations were stable for
Epithelial, Endothelial, Proliferating, Myeloid, and Lymphoid states, while
Fibroblast/stromal was representation-unstable. All six major compartments
differed stably between early and late irradiated mice. Fibroblast-subtype
conclusions were more representation-sensitive.

## Limitations

Independent mice support replicate-aware association, not causality. Early and
late groups contain different mice, so the comparison is cross-sectional and
does not demonstrate longitudinal reversal or persistence. Exact permutation
resolution is discrete, related contrasts reuse controls, annotations are
marker-inferred, and UMAP is display-only.
