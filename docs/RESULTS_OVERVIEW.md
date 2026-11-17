# Guide to the published results

Use this page to connect the paper's evidence with its source workflow. The
summaries describe frozen analyses; exact numerical claims should be checked
against the associated evidence ledger and figure-source table.

[Read the paper](https://doi.org/10.1093/bioadv/vbag266) ·
[Workflow catalog](DATASET_MAP.md) · [Reproduce a workflow](REPRODUCIBILITY_GUIDE.md) ·
[Cite and reuse](PUBLICATION_GUIDE.md)

## What each analysis contributes

| Analysis | Contribution | How to read it |
|---|---|---|
| GSE280305 | Biological motivation: changes in hematopoietic state geometry during post-irradiation recovery | Time-point differences relative to a post-irradiation reference |
| Synthetic benchmark | Controlled tests of estimator robustness, uncertainty, representation corruption, local geometry, and dynamics agreement | Behavior under the specified simulation and held-out evaluation protocol |
| Pancreatic development | A distinct developmental context for geometry–dynamics comparison | Descriptive validation with related velocity-derived comparators |
| GSE249479 | Human inflammatory HSPC perturbations and complementarity with separability/abundance measures | Descriptive evidence without recoverable biological replication |
| GSE211713 | Mouse-level uncertainty and representation-sensitive early/late radiation associations | Replicate-aware association in a cross-sectional design |

## Original hematopoietic case study: GSE280305

The original workflow follows hematopoietic state-space changes across
post-irradiation time points and compares geometric displacement with RNA
velocity. D8 is the earliest sampled post-irradiation reference, so reference
deviation describes departure from that sampled state rather than distance
from an untreated healthy baseline. The original Scanorama representation is a
case-study coordinate system; the public validation workflows assess
representation dependence separately.

Start with the [manuscript guide](../notebooks/manuscript/README.md). Original
reference preparation and annotation include semi-manual decisions documented
in the preparation/exploration notebooks.

## Synthetic benchmark

The [benchmark workflows](../notebooks/benchmarks/README.md) keep calibration
and held-out evaluation separate. They examine robust estimators, state
shifts, representation and local-geometry diagnostics, dynamics agreement, and
framework ablation against stored truth.

Retain the documented negative findings when discussing performance:
balanced-replicate seed dependence, incomplete state-level instability recall,
limited local-distortion recall (including zero recall in the explicit
corruption jobs recorded by the protocol), and imperfect bootstrap coverage.
These results define the tested behavior and failure modes;
they do not establish biological truth.

## Public pancreatic development

The [pancreas workflow](../notebooks/public_validation/pancreas/README.md)
compares representation geometry with public developmental dynamics and
negative controls. scVelo and CellRank provide related dynamical context;
CellRank is not independent evidence from RNA velocity. Conclusions are
descriptive and do not establish causal transitions.

## Human HSPC inflammation: GSE249479

The [GSE249479 workflow](../notebooks/public_validation/gse249479/README.md)
retained **34,432 cells**. No recoverable biological-replicate identity was
available, so every result remains **`descriptive_only`**.

Frozen evidence includes TNF-associated shifts in Activated HSC and Lymphoid,
broader LPS displacement, an abundance-dominant LPS HSC/quiescent pattern, high
rank agreement among the primary representations, and diffusion-map
sensitivity with limited influence on the primary conclusions. These evidence
layers answer complementary questions rather than forming a single composite
perturbation score.

Official R Augur is the principal comparator. The retained Python method is an
Augur-inspired approximation used only for supplementary implementation
sensitivity.

## Mouse-lung radiation: GSE211713

The [GSE211713 workflow](../notebooks/public_validation/gse211713/README.md)
retained **131,157 cells from 20 independent mouse libraries**. Primary
comparisons use biological mouse-level centers and uncertainty estimates.

Under the frozen consensus rule, early 17 Gy versus control was
representation-sensitive. Late 17 Gy versus control was stable in **five of six
major compartments**, and all six major compartments differed stably between
early and late irradiated mice. Fibroblast/stromal and fibroblast-subtype
conclusions showed greater representation sensitivity.

These are **cross-sectional associations**: different mice were sampled in the
early and late groups. They do not establish causality, persistence, reversal,
longitudinal change, or within-mouse change.

## Find the evidence behind a figure

The workflow's numbered notebooks lead to its figure sources and metadata.
After the corresponding analyses and assembly have run, the principal
locations are:

| Evidence | Generated location relative to the repository |
|---|---|
| Synthetic figures, source tables, and review copies | `results/revision_synthetic_benchmark/` |
| Pancreas outputs | `results/public_validation/pancreas_dataset_d/` |
| GSE249479 manuscript ledger | `results/public_validation/gse249479_dataset_b/manuscript/dataset_b_evidence_ledger.csv` |
| GSE211713 manuscript ledger | `results/public_validation/gse211713_dataset_c/manuscript/dataset_c_evidence_ledger.csv` |
| Cross-workflow evidence ledger and final inventory | `results/revision_finalization/` |

These paths describe generated artifacts, most of which are ignored by Git and
are not included in a fresh clone. Source implementations and checksum policies
are linked from the [workflow catalog](DATASET_MAP.md). Figure assembly reads
frozen numerical inputs; changing labels or layout should not alter those
inputs or the scientific conclusion they support.
