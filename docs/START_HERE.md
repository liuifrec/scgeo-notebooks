# Start here

## Choose the workflow

- **Package user:** install
  [`scgeo`](https://github.com/liuifrec/scgeo) and follow its quick start.
- **Reviewer:** begin with [Results overview](RESULTS_OVERVIEW.md), then follow
  links to the relevant evidence ledger and output-free source notebook.
- **Reproducer:** read the dataset folder README, create the named environment,
  set supported path overrides, and use the clean-kernel executor.
- **Manuscript editor:** use frozen evidence ledgers and assembly outputs; do not
  rerun numerical analyses to adjust presentation.

## Repository conventions

Numbered notebooks show workflow order. Thin notebooks call implementations in
`scripts/`; tracked notebooks have no embedded outputs. Generated review copies
and large artifacts belong under ignored result or workstation data paths.

## Scientific interpretation

Always identify the analysis unit and inference scope. Synthetic jobs support
controlled method evaluation. The pancreas workflow is descriptive. Dataset B
is `descriptive_only`. Dataset C primary contrasts use independent mouse
libraries but are cross-sectional. Representation agreement does not turn
nested PCA views into independent evidence.

Next: [Dataset map](DATASET_MAP.md), [Reproducibility guide](REPRODUCIBILITY_GUIDE.md),
and [Thin notebook architecture](THIN_NOTEBOOK_ARCHITECTURE.md).
