# Reproducibility guide

## Pinned repositories

The revision evidence package uses ScGeo commit
`9a0ed16cbaa57f935f9c9bc87d1643a25b51012c`. Verify the package commit, the
companion branch, numerical checksums, and a clean worktree before executing a
frozen workflow.

## Supported path overrides

The public-data overrides cover GSE249479, a public human HSPC inflammatory
xenograft dataset, and GSE211713, a public mouse-lung radiation dataset with 20
independent mouse libraries.

| Variable | Purpose |
|---|---|
| `SCGEO_SOURCE_REPO` | ScGeo package checkout |
| `SCGEO_BENCHMARK_DIR` | Frozen synthetic benchmark directory |
| `SCGEO_GSE249479_H5AD` | GSE249479 source H5AD |
| `SCGEO_GSE249479_COMPACT_H5AD` | GSE249479 compact QC/HVG H5AD |
| `SCGEO_GSE249479_OUTPUT_DIR` | GSE249479 generated output root |
| `SCGEO_GSE211713_DATA_DIR` | GSE211713 local data root |
| `SCGEO_GSE211713_COMPACT_H5AD` | GSE211713 compact QC/HVG H5AD |
| `SCGEO_GSE211713_ANNOTATED_H5AD` | GSE211713 annotated H5AD |
| `SCGEO_GSE211713_REPRESENTATION_V2_H5AD` | GSE211713 canonical representation object |
| `SCGEO_GSE211713_C6_V2_OUTPUT_DIR` | GSE211713 representation-v2 outputs |
| `SCGEO_GSE211713_C7_V2_OUTPUT_DIR` | GSE211713 ScGeo-v2 outputs |

Set only variables supported by the invoked script. Dataset-specific executors
also record their resolved inputs and checksums.

## Frozen workstation paths

Some validated executors and gates still contain absolute paths under
`/home/liuyuchen`. Examples include GSE211713 v2 validation/execution,
`scripts/finalize_revision_evidence.py`, the recorded scVI interpreter, and
GSE249479 execution/comparator wrappers. These paths reproduce the frozen
workstation audit. This documentation pass does not rewrite them; portability
changes require a separate implementation and numerical-equivalence review.

## Source and generated notebooks

Tracked source notebooks must contain no outputs and null execution counts.
Executors write clean-kernel review copies beneath ignored `results/`
directories. H5ADs, models, downloads, caches, and generated figures are not
committed unless repository policy explicitly tracks a small manifest or
documentation artifact.

## Validation checklist

1. verify Git commits and worktree cleanliness;
2. parse JSON/YAML configs;
3. verify source notebook outputs and code-cell immutability;
4. validate numerical input/output checksums;
5. reject duplicate canonical keys and invalidated-run inputs;
6. run `git diff --check` and the relevant test suite.
