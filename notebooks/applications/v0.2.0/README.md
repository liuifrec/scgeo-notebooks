# ScGeo public biological applications (B and C)

This separate, opt-in track redraws previously computed public-data evidence with
the released plotting API. It does not replace the frozen paper reproduction
workflows. It requires the original, checksum-matched local H5ADs and revision
output tables; it does not download or refit them.

| Notebook | Design | Views |
| --- | --- | --- |
| [04 GSE249479 HSPC inflammation](04_gse249479_inflammation.ipynb) | 34,432 cells, PBS/TNF/LPS, descriptive only | Matched condition/score facets, relative density, annotation associations, frozen effects, representation matrices, neighbor-adjacency pies, stored compass, geometry/evidence atlas |
| [05 GSE211713 mouse lung](05_gse211713_radiation_response.ipynb) | 131,157 retained cells from 20 libraries; primary contrasts use 13 mice (5/4/4) | Primary and fibroblast facets, stored scores, relative density, individual qualified mouse centers, frozen effects and representation matrices, equal-mouse pies, stored compass and geometry/evidence atlas |

These are runnable sources, not pre-executed public figures. The two datasets
are the first applications; the original hematopoietic recovery, held-out
simulation benchmark and pancreas workflows remain separate until their
verified local objects have been recovered. The current tutorial series
remains the download-free starting point.

## Set up

Use an independent Python 3.12 environment. From this companion checkout:

```bash
python -m pip install -r environment/scgeo_v020_tutorials.txt "scgeo[scanpy,notebook]==0.2.0"
python -m pip check
python -m ipykernel install --user --name scgeo-v020-applications --display-name "ScGeo 0.2.0 Applications"
```

Set the following local paths explicitly (replace placeholders; do not commit
local paths or the input files):

```bash
export SCGEO_APPLICATION_B_ROOT=/path/to/gse249479/results
export SCGEO_APPLICATION_C_ROOT=/path/to/gse211713
export SCGEO_APPLICATION_FROZEN_RESULTS=/path/to/original/scgeo-notebooks/results/public_validation
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
unset PYTHONPATH PYTHONHOME
python scripts/execute_scgeo_v020_applications.py --validate-only
python scripts/execute_scgeo_v020_applications.py --dataset B --output-dir /path/to/new/results
python scripts/execute_scgeo_v020_applications.py --dataset C --output-dir /path/to/new/results
```

The output directory must be outside all input directories. Use the original
revision outputs rather than a new empty clone's results directory. Every run
creates a new directory. To run interactively, additionally set
`SCGEO_APPLICATION_OUTPUT_ROOT` to an existing separate output directory
before starting JupyterLab, select the registered kernel and run one notebook.
The runner retains partial executed notebooks on failure. Open the executed
notebooks in its reported run folder to see the figures without rerunning.

## Opt-in 0.2.1 visualization profile

The default runner remains pinned to the published 0.2.0 package. The same
notebook sources also accept `--version 0.2.1` in a separate environment.
The tested pre-release source commit is
[683228b](https://github.com/liuifrec/scgeo/commit/683228b1c52808d3cc5df2de5f6e104932c29328).
The maintainer's GSE211713 preview of that commit completed 33 plot comparisons,
99/99 selected unchanged-source checks and zero collapsed-layout warnings.
This preview used the old application source and a separate preview driver;
the integrated notebook profile still requires full real-data execution.

For candidate testing only, install into a fresh environment (not the
published 0.2.0 environment), using:

```bash
python -m pip install "scgeo[scanpy,notebook] @ git+https://github.com/liuifrec/scgeo.git@683228b1c52808d3cc5df2de5f6e104932c29328" nbformat nbclient ipykernel
python -m pip check
python -m ipykernel install --user --name scgeo-v021-applications --display-name "ScGeo 0.2.1 Applications"
export SCGEO_APPLICATION_EXPECTED_SHA=683228b1c52808d3cc5df2de5f6e104932c29328
python scripts/execute_scgeo_v020_applications.py --version 0.2.1 --kernel scgeo-v021-applications --dataset C --output-dir /path/to/new/results
```

For interactive notebooks set `SCGEO_APPLICATION_VERSION=0.2.1` before
launching Jupyter. When using the eventual published wheel, unset
`SCGEO_APPLICATION_EXPECTED_SHA`; do not claim a VCS commit for a wheel
without separate release provenance. The profile uses adaptive compasses and
atlases, larger pies, stable category colors, shared heatmap/effect scales and
PC1 axes within each declared hierarchy/contrast family. The abundance axis
uses the maximum observed state fraction among all retained samples in that
hierarchy with headroom (minimum 5 percentage points, maximum 100%). This is
a display-only improvement over the preview's full 0–100% range, not a new
denominator or inference. Numerical source exports remain unchanged in meaning.
No new estimator or private data is introduced.

## Molecular shift explanations (0.2.1)

A separate opt-in track applies the actual `robust_shift`, `explain_shift`,
`shift_curves` and `shift_embedding` APIs to the original public count data.
See [molecular shift guide](MOLECULAR_SHIFTS.md) for strict feature-to-coordinate
verification, coverage rules and execution instructions. The existing geometry
notebooks remain unchanged.

## Inputs and safeguards

- B: `gse249479_revision_representations.h5ad` and
  `gse249479_revision_qc_hvg.h5ad`. The source evidence lives under
  `gse249479_dataset_b/scgeo/04_{tnf,lps}_detailed_representation_state_evidence.csv`.
- C: `gse211713_revision_representations_v2.h5ad` (not the invalidated v1),
  `gse211713_revision_qc_hvg.h5ad` and `gse211713_revision_annotated.h5ad`.
  Use `gse211713_dataset_c/scgeo_v2/per_contrast/<contrast>__<hierarchy>/`,
  with `state_evidence.csv`, `individual_mouse_centers.csv`, and
  `report.json`. These are the original contrast-specific files, not
  invalidated or aggregated substitute results.

Five exact H5AD SHA256 values and shapes are pinned in the helper. Cell
identities and order must match. Expression matrices are opened backed and
read-only; only annotation metadata are used. The B evidence CSVs are hashed
and cross-checked against recorded cell coverage and labels, but no independent
expected CSV hash is asserted. The C CSVs must match their own original report
hashes and exact per-state eligible and excluded mouse membership/counts.
These checks verify consistency; the original report is local source
provenance, not an independently signed artifact. All source hashes are
rechecked after rendering.

The B graph is **new**: Scanpy neighbors on frozen PCA50 (30 neighbors,
Euclidean, seed 1729). B builds this graph in a temporary coordinate-only
AnnData, preserving the original object (which stores diffusion coordinates but
not diffusion eigenvalues). No missing diffusion metadata is invented. C uses
its existing PCA50 graph induced on the primary
roster. Both aggregate symmetric edge mass using state membership matrix H
and cell graph W: M = H.T W H, then divide each off-diagonal M[i,j] by
sqrt(sum(M[i,:]) * sum(M[j,:])). The diagonal is zero. The PAGA-style renderers
display this descriptive adjacency; **this is not Scanpy PAGA inference**.
No velocities, OOD scores, or causal drivers are fabricated. Graph context
includes all retained cells in the supplied groups. In the fibroblast view,
non-target cells are pooled into an explicit context category, not an inferred
biological state. The sample-abundance pie uses the complete retained-cell
denominator, not a fibroblast-only or qualified-state subset.

B has no valid biological replicate IDs. No sample bootstrap, biological
intervals or pseudo-replicates are introduced. C retains the frozen eligibility
criteria (30 cells per mouse for major states, 20 for fibroblasts, at least
three eligible mice per group). Sample-center plots show only PC1, with no
inferred pairing. The compass separates stored PCA50 magnitudes from display
connectors. The normalized effect intervals use the original normalized
interval columns. Representations are not independent confirmations.

Neither dataset has a saved expression-to-coordinate map in the H5ADs.
Signature/marker score facets are allowed, but exact signed gene contributions
remain deferred until original feature order, preprocessing and PCA/SVD
weights can be recovered and verified. These applications do not establish
a lineage, recovery mechanism or causal regulation.

## Outputs and review

Each plot is exported to PNG, SVG and PDF with its source tables and
PlotResult manifest. The notebook verifies each exported hash and saves an
application manifest containing source digests, package versions, graph origin
and the companion commit. Keep executed notebooks and large outputs in the
separate local results directory. Review figures, manifests and any public
cell/sample identifiers before publishing selected gallery previews. Do not
commit raw H5ADs, local paths, or unreviewed executed outputs. No real-data CI
execution is implied by the source validation.
