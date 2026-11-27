# GSE280305: native plotting application (ScGeo 0.2.1)

The goal of this example is to demonstrate the **native ScGeo plotting API** on the recovered public data, not to expand the biological discovery analysis or replace the original paper. The checkpoint and existing descriptive geometry, molecular and local-geometry results are read-only and checksum-gated.

[Open the application notebook](08_gse280305_native_visualization.ipynb).

## What runs

| Native function | Demonstration |
| --- | --- |
| `embedding_facets` | Four matched time-point panels in each of the two original embeddings |
| `representation_stability_heatmap` | Stored D8–D21 normalized displacements and consensus, with local diagnostic flags |
| `consensus_state_map` | Stored consensus classifications on Scanorama UMAP |
| `local_distortion_map` | State-median neighbor preservation and local shape distortion; no per-cell values manufactured |
| `state_evidence_panel` | Stored raw-PCA effects, consensus and full-roster local diagnostics |
| `geometry_dynamics_atlas` | Adaptive state/evidence layout, explicitly without a dynamics panel |
| `shift_curves` | Native expression/contribution maps and pooled curves for selected genes in clusters 4 and 9 |
| `shift_embedding` | Schematic center connectors using the stored analysis weights |

Each figure uses `return_result=True` and `scgeo.pl.export_figure` to produce PNG, SVG, PDF, source CSVs and metadata. Full output stays local. The molecular plotting API requires a compatible in-memory explanation and cell weights, so the notebook recreates only the D8–D21 explanation, compares all 42,000 gene rows and 14 state estimates with the reviewed exports, and stops on a mismatch. It does not refit PCA or Scanorama, recompute velocity, or conduct new marker screens.

## Run

Use the existing isolated ScGeo 0.2.1 environment (tested candidate source: `683228b1c52808d3cc5df2de5f6e104932c29328`). Set the five paths explicitly:

```bash
export SCGEO_A_H5AD=/path/to/GSE280305_scgeo_phase1_qc.h5ad
export SCGEO_A_GEOMETRY=/path/to/gse280305_geometry_v1
export SCGEO_A_MOLECULAR=/path/to/gse280305_molecular_v1
export SCGEO_A_LOCAL=/path/to/gse280305_local_review
export SCGEO_A_OUTPUT=/path/to/separate/native-plot-results
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
```

Open the notebook from the companion checkout with the matching kernel. Alternatively, execute through `nbclient` using a kernel registered with the exact environment interpreter. The notebook creates a fresh result directory and verifies the reviewed input hashes before and after execution.

## Evidence boundaries

D8 is a post-irradiation reference. There is one recorded library per time point, with no verified biological replicate identities. All analyses are descriptive and no biological-sample bootstrap intervals are included. Cluster IDs are retained rather than inventing biological labels. The original scaled-HVG expression matrix defines the PCA feature units; molecular contributions are not log fold changes. Local geometry was assessed on all four time points at k=15, and the state-level maps display those medians. Curves are not temporal trajectories. Original velocity vectors are not available, so velocity-dependent panels are explicitly deferred.

## Publishing

After the real-data run succeeds, review the native PlotResult PNG/SVG images and their source tables and metadata. Publish a compact selection of those **actual native figures** in the repository's results-first gallery, with links to this notebook. Do not publish the large H5AD or unreviewed cell-level exports. The workstation run completed all six code cells, all eleven native plotting calls and the 42,000-row molecular equivalence check. See the [results-first gallery](../results/gse280305/README.md) for unchanged native SVGs, original plot metadata, the aggregate state report and review scope. The tracked source notebook remains output-free.
