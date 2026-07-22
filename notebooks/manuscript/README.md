# Original manuscript notebooks

## Purpose and inference status

These notebooks assemble the original GSE280305 manuscript analyses after data
preparation and semi-manual reference annotation. Their biological scope is
defined by the original recovery design and should not be generalized to the
public revision datasets.

## Notebook order

1. `05_OOD.ipynb` — retained filename; interpret the analysis as reference
   support, representation instability, or extrapolation sensitivity rather
   than general OOD detection.
2. `Velocity_shift_alignment.ipynb`
3. `Final_summary.ipynb`
4. `test_driver_genes.ipynb`

`lite_generate.ipynb` and `Tif_convert.ipynb` are presentation utilities rather
than independent biological analyses.

## Implementation, inputs, and outputs

Inputs are prepared under `notebooks/data_prep/`, followed by the semi-manual
`notebooks/exploration/06_Ref_prep.ipynb` and
`07_Reference_based_annotation.ipynb` steps. Notebook-local code and companion
scripts generate manuscript figures and intermediate files under local result
paths.

## Accepted findings and limitations

The notebooks preserve the frozen manuscript analyses; this README does not
restate numerical conclusions. Reference preparation and annotation contain
semi-manual decisions, filenames preserve historical terminology, and weak
reference support is not evidence of general out-of-distribution detection.
