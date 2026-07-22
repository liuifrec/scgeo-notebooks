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
controlled method evaluation. The public pancreatic-development workflow is
descriptive. GSE249479, a public human HSPC inflammatory xenograft dataset, is
`descriptive_only` because no recoverable biological-replicate identity is
available. GSE211713, a public mouse-lung radiation dataset with 20 independent
mouse libraries, supports replicate-aware primary contrasts that are
cross-sectional rather than longitudinal. Representation agreement does not
turn nested PCA views into independent evidence.

Next: [Dataset map](DATASET_MAP.md), [Reproducibility guide](REPRODUCIBILITY_GUIDE.md),
and [Thin notebook architecture](THIN_NOTEBOOK_ARCHITECTURE.md).
