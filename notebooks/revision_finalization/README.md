# Revision-finalization notebooks

## Purpose and inference status

These output-free notebooks assemble reviewer-facing definitions, inventories,
terminology corrections, and reproducibility evidence from frozen artifacts.
They are presentation and audit workflows, not new numerical analyses.

## Notebook order

1. `01_reviewer2_metric_repairs.ipynb`
2. `02_distribution_and_aggregation_audit.ipynb`
3. `03_terminology_and_claim_audit.ipynb`

## Implementation, inputs, and outputs

The testable implementation is `scripts/assemble_reviewer2_repairs.py`, with
final evidence assembly in `scripts/finalize_revision_evidence.py`. Inputs are
frozen synthetic, pancreas, Dataset B, official R Augur, and Dataset C evidence
artifacts. Generated tables, figures, captions, alt text, ledgers, and audit
reports are written under ignored `results/revision_finalization/` paths.

## Accepted findings and limitations

The outputs clarify mixscore aggregation, distribution comparison units,
representation metrics, ScGeo–scPASI distinctions, and corrected claims. They
must retain negative evidence, distinguish official R Augur from the Python
approximation, keep Dataset B descriptive, and keep Dataset C cross-sectional.
They do not change any numerical result or replace scientific manuscript and
response-letter writing.
