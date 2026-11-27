# Molecular shift explanations on public datasets

This is a separate **new analysis** using ScGeo 0.2.1, not a replacement for
the original manuscript results or the existing successful geometry gallery.

| Notebook | Coverage-qualified outputs | Predeclared focus for figures |
| --- | --- | --- |
| [06 — GSE249479](06_gse249479_shift_explanations.ipynb) | Seven detailed states × two contrasts (coverage permitting) | HSC/quiescent, activated HSC, MPP/progenitor |
| [07 — GSE211713](07_gse211713_shift_explanations.ipynb) | Six major and three fibroblast states × three primary contrasts | Fibroblast/stromal, myeloid, the three refined fibroblast states |

For each supported state, the workflow calculates a fresh robust point
estimate on the verified linear representation, then calls the released
explain_shift API. If the basis is frozen, the new point estimate must match
the stored frozen displacement (relative tolerance 1e-4, absolute 1e-5).
The signed gene sum must reconstruct the new observed displacement.
No frozen source or original estimator implementation is changed.

## Exact expression and mapping

The checksum-gated input paths are the same as in the [application guide](README.md).
Both original expression matrices store raw counts. They are read backed in
sparse row blocks and subset to their original ordered HVGs. **Normalization
differs between datasets**.

- B: 3,000 HVGs, divide by the original per-cell total_counts, multiply by
  10,000, log1p. The original *uncentered* randomized SVD checkpoint
  gse249479_dataset_b/representations/03a_pca_representations.npz includes
  loadings, feature order, cell order and PCA50 scores. The corresponding
  audit/03a_pca_summary.json records its digest and the source count digest.
  Both are verified; checkpoint arrays must match the frozen H5AD.
- C: first 3,500 flagged HVGs in original feature order, divide by the
  retained HVG count sum, multiply by 10,000, log1p. The corrected v2 centered
  ARPACK PCA stores scores but not loadings. If available,
  c6_v2/metadata/pca_stage.json must match the feature-order hash. For
  centered PCA, with normalized matrix Xc, score matrix Z and loading L,
  Xc.T @ Z = L @ (Z.T @ Z). Recover L from the cross-product in bounded
  sparse blocks. This is an algebraic recovery attempt, not a claim that
  the saved scores alone identify loadings for arbitrary nonlinear maps.

For either basis, check loading orthogonality (maximum absolute Gram error
5e-4) and all-cell affine projection residual against 1e-6 + 1e-5 times
maximum centered coordinate magnitude. The B intercept must also be near
zero. ScGeo independently checks the selected contrast's affine mapping
again. Original feature IDs and readable feature labels are exported,
together with the recovered loadings and intercept.

### Verified global error budget for recovered C loadings

The original state-relative default may reject a narrow state even after the
full-data map has passed, because both the coordinate scale and the fitted
intercept change when the cells are subset. With global residual
r_i = z_i - x_i L - b, and E = max_i ||r_i||_infinity, the residual after
recentering on any qualified subset S is r_i - mean_S(r), bounded by 2 E.
Only for the globally verified recovered centered PCA map, use a fixed
absolute state tolerance of 2 E + 1e-6 and relative tolerance zero. The
additional 1e-6 is the original absolute numerical allowance. E is measured
over every frozen cell before any state explanation; it must satisfy the
original global threshold and orthogonality check. No failed state can
inflate this budget, and no fallback or representation refit is triggered
by a state-level failure. B's original SVD and explicitly refitted bases
retain the original per-state default.

The mapping audit records the policy, global scale, residual and threshold.
The state summary records actual local residual, applied threshold, original
state-relative threshold and policy. The original coordinates, qualified
mouse roster, frozen point-estimate comparison and signed contribution
conservation check are unchanged. This is finite-precision molecular
accounting on a verified approximate affine map, not an assertion of
bitwise-exact recovery of the historical PCA model.

If the original B checkpoint is absent or a verified projection cannot be
recovered, the default is to **stop**. A checksum or identity mismatch always
stops. To deliberately permit a fresh linear fit instead, set
SCGEO_MOLECULAR_ALLOW_REFIT=1. This fallback uses the original count
normalization but produces new scores and loadings, verifies their map, and
labels all outputs new_representation_not_frozen. It must not be presented
as an explanation of the original frozen geometry. Prefer the strict run
first. The all-cell acceptance tolerance is never relaxed.

## Scientific interpretation

A positive gene contribution supports the observed treated-minus-reference
direction, a negative one opposes it. Values are in the linear representation
units, not log fold changes or causal-driver scores. Their sum reconstructs
the shift magnitude. Genes defining an expression-derived axis are not
independent evidence about that axis.

The original signature and annotation marker definitions are reused as
predefined programs, restricted to unambiguously matched HVGs with at least
two members. The membership audit includes absent and ambiguous members.
Program sums are not enrichment tests, and overlapping programs do not
partition the displacement. For focus-state plots, select the largest
positive and negative gene contributions and the program with the largest
absolute sum. These plots are explicitly selected for illustration, not
independent confirmation of a response. The shift curve x-axis is a geometric
coordinate, not time or inferred pseudotime. Schematic connectors are not
velocity or lineage trajectories.

B is pooled descriptive analysis with no invented replicate IDs or sample
bootstrap intervals. C uses exactly the qualified mouse IDs and cell counts
from the frozen state-specific reports (30/20 cells per mouse for major/
fibroblast states; at least three eligible mice per group). Each state is
analyzed separately, with all qualified cells, and 200 conditional sample
bootstrap draws for molecular contributions. Curve ribbons show observed
sample ranges, not confidence intervals. Unqualified states are recorded as
skipped; no unsupported shift direction or missing expression is invented.
There is no local-key subdivision, regulatory network inference, or
causal validation. Refined fibroblast states are separate contrasts, not
an exact decomposition of a parent displacement.

## Run on the existing 0.2.1 application environment

From the companion repository root, set the input paths and register the
kernel as in README.md. The exact tested source environment from the 0.2.1
preview may be reused; with a published distribution, unset the VCS-source
guard unless independent provenance is available.

~~~bash
unset SCGEO_MOLECULAR_ALLOW_REFIT
python scripts/execute_scgeo_v020_applications.py --validate-only
python scripts/execute_scgeo_v020_applications.py --molecular --version 0.2.1 \
  --kernel scgeo-v021-applications --dataset B --timeout 7200 \
  --output-dir /path/to/new/molecular-results
python scripts/execute_scgeo_v020_applications.py --molecular --version 0.2.1 \
  --kernel scgeo-v021-applications --dataset C --timeout 7200 \
  --output-dir /path/to/new/molecular-results
~~~

The timeout is per notebook cell. Separate dataset runs avoid keeping both
large expression matrices in memory. The original data and geometry
outputs are read-only. The runner checks the plot-export manifests and
molecular table hashes, and retains partial executed notebooks on failure.
Inspect mapping_audit.json and molecular_state_summary.csv before reviewing
the gene/program curves. Full gene and program tables are gzip-compressed.
Do not commit raw data or unreviewed cell-level outputs.

Real-data execution and mapping recovery are not established by the
synthetic contract tests. These workflows require the original local
files and explicit workstation verification.
