# Reuse ScGeo's plotting workflows

The package now has development interfaces for evidence panels, contrast
comparisons, biological sample views, and reproducible figure export. This
companion entry point redraws stored tables so presentation changes can be
reviewed separately from numerical analyses.

[Start here](START_HERE.md) ·
[Package plotting guide](https://github.com/liuifrec/scgeo/blob/main/docs/PLOTTING.md) ·
[Visual gallery](https://github.com/liuifrec/scgeo/blob/main/docs/PLOTTING_GALLERY.md) ·
[Reproducibility](REPRODUCIBILITY_GUIDE.md)

## Use the matching development checkout

These features are **unreleased** and are not included in `scgeo==0.1.0` on
PyPI. Install the ScGeo source checkout containing the plotting guide and record
its commit. In an environment with that checkout installed:

```bash
python -m pip install -e /path/to/scgeo
python /path/to/scgeo/examples/plotting_gallery.py --output-dir results/plotting_gallery
```

The gallery is entirely synthetic. It produces actual figures, CSVs, and JSON
manifests through the package API without downloading manuscript datasets.
Related rotated representations illustrate sensitivity; they are not
independent validation experiments.

## Redraw an archived state report

Use the source CSV from the analysis you intend to present:

```bash
python scripts/render_plot_tables.py state path/to/state_report.csv \
  --output-stem results/presentation/state_evidence \
  --metric delta_norm \
  --title "State evidence"
```

The renderer reads the table once and records its filename and SHA-256 digest.
It exports PDF, SVG, PNG, the plotted source tables, and a JSON manifest beneath
the chosen output stem. Optional `--metadata analysis_metadata.json` retains
study identifiers, comparison direction, representation, bootstrap unit,
analysis commit, and original input checksums. A CSV alone does not recover
metadata that its original exporter omitted.
Choose an output directory different from the input table directory; the
renderer enforces this to protect archived source tables.

For example, a state-report metadata file can supply the comparison labels and
the exact quantitative representation without guessing them from column names:

```json
{
  "package_commit": "record-the-exact-commit-used-for-your-analysis",
  "provenance": {
    "group0": "control",
    "group1": "treated",
    "condition_key": "condition",
    "effect_representation": "X_pca",
    "sample_key": "sample",
    "bootstrap_unit": "sample"
  }
}
```

Replace these values with the recorded analysis settings. Omitted provenance
remains missing. Literal identity labels such as a state named `NA` are retained;
missing scale or interval declarations do not become valid provenance labels.

The evidence panel preserves stored interval endpoints and distinguishes
availability from a neutral or unstable result. Raw magnitude intervals are
never assigned to normalized effects. The display cannot turn magnitude
uncertainty into a hypothesis test.

## Compare multiple contrasts

Prepare one tidy table with these fields:

| Column | Meaning |
|---|---|
| `state` | Cell-state label |
| `contrast` | Explicit comparison and direction |
| `metric` | Measured quantity, such as `delta_norm` |
| `estimate` | Stored effect estimate |
| `interval_low`, `interval_high` | Optional paired endpoints in the estimate's units |
| `interval_type` | Required for intervals; e.g. `sample_bootstrap_ci95` or `representation_iqr` |
| `representation`, `unit` | Declared comparable coordinate system and units |
| `scale_id` | Optional explicit identifier of a common quantitative scale |
| `evidence_status` | Optional assessed, insufficient, unavailable, or other recorded status |

```bash
python scripts/render_plot_tables.py effects path/to/contrast_effects.csv \
  --metric delta_norm --output-stem results/presentation/contrast_comparison
```

Shared axes require a declared compatible scale. A PCA name alone does not
make different datasets' coordinate distances comparable. Duplicate
state/contrast rows must be resolved explicitly before rendering.

Existing GSE211713 assembly tables include across-representation summaries.
Where a source table labels `normalized_magnitude_q25/q75`, those endpoints
describe representation quartiles; do not rename them as biological-sample
bootstrap confidence intervals. GSE249479 remains descriptive when valid
biological-replicate identity is unavailable.

## Inspect genuine biological samples

New sample-aware `robust_shift` runs retain exact sample-center matrices:

```python
samples = sg.get.sample_centers(adata)
samples.to_csv("results/sample_centers.csv", index=False)
```

```bash
python scripts/render_plot_tables.py replicates results/sample_centers.csv \
  --coordinate 0 --output-stem results/presentation/sample_coordinate_0
```

Coordinates are zero-based. One point represents one biological sample; cell
counts are not replicate weights. The displayed mean or median is descriptive.
The CLI draws independent points. If an analysis has genuine pairing, use the
package API with an explicitly joined `pair_key`; never infer pairing from
similar sample names or connect independent cross-sectional mice.

Old results lacking sample-center storage cannot supply these plots from
aggregate estimates. Retain the original result and rerun only in a separate,
documented analysis output when the source data and parameters are available.

## Representation and dynamics views

Render a tidy per-representation table with `node`, `rep`, and the chosen metric:

```bash
python scripts/render_plot_tables.py representation path/to/per_rep_state.csv \
  --metric alignment_cosine --output-stem results/presentation/representations
```

The package's `velocity_shift_alignment(..., mode='stored', analysis_rep=...)`
keeps quantitative alignment separate from its display overlay. It requires
the stored alignment table and, for arrows, compatible display velocity
coordinates. A presentation change to UMAP does not recompute the reported
alignment.

## Relate the new views to the manuscript

Use these figures for new analyses, supplementary inspection, or clearly
identified presentation adaptations. The original assembly scripts and frozen
sources remain the route to reproduce the published figures. Record both
package and companion commits, retain source hashes, and preserve the design's
inference scope. GSE280305's D8 reference is post irradiation; GSE211713 is
cross-sectional; nested representations share information.

For credits and Methods reporting, use the [publication guide](PUBLICATION_GUIDE.md).
