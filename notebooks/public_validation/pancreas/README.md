# Public pancreatic-development workflow

## Purpose and inference status

The public pancreatic-development workflow uses the official CellRank pancreas dataset to compare ScGeo
representation–dynamics diagnostics with scVelo and CellRank context. It is a
descriptive dynamics validation; CellRank is derived from RNA-velocity context
and is not an independent biological confirmation.

## Notebook order

1. `00_acquire_and_validate.ipynb`
2. `01_scvelo_dynamical_velocity.ipynb`
3. `02_cellrank_fates.ipynb`
4. `03_scgeo_representation_dynamics.ipynb`
5. `04_manuscript_figures.ipynb`

## Implementation, inputs, and outputs

Use `scripts/execute_pancreas_validation.py` and
`scripts/pancreas_validation_common.py`. The official public input and its
checksum are recorded by the acquisition stage. Generated H5ADs, executed
notebooks, figures, figure-source tables, alt text, and metadata are written
under ignored `results/public_validation/pancreas_dataset_d/` paths.

## Frozen findings and limitations

The accepted evidence concerns representation-dependent agreement with public
developmental dynamics and its negative controls. It does not establish causal
transitions, independent confirmation beyond the velocity-derived comparator,
or universal representation preservation.
