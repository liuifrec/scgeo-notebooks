# How the notebooks and scripts fit together

ScGeo's companion keeps scientific explanation close to each workflow while
placing reusable execution logic in scripts. This makes an analysis easier to
read, review, and adapt without storing large outputs in notebook source.

[Start here](START_HERE.md) · [Workflow catalog](DATASET_MAP.md) ·
[Reproducibility](REPRODUCIBILITY_GUIDE.md)

## Three parts of a workflow

| Part | Role | Where it lives |
|---|---|---|
| Source notebook | Explains the biological question, inference scope, inputs, outputs, and ordered analysis steps | Tracked under `notebooks/` |
| Implementation | Performs reusable calculations, validation, execution, or figure assembly | Tracked under `scripts/`, with settings under `configs/` |
| Generated bundle | Holds executed copies, figures, source tables, alt text, metadata, and logs | Ignored result or local-data directories |

Most revision notebooks are short orchestration documents that call the
matching script. Some original exploratory/manuscript notebooks contain more
notebook-local implementation and semi-manual decisions; follow their folder
guide rather than assuming every notebook has the same structure.

## Released API tutorials

The separate [0.2.0 learning track](../notebooks/tutorials/v0.2.0/README.md) uses explicit API cells alongside byte-identical, checksum-recorded examples from the release source. Its execution runner starts fresh kernels, writes executed copies only to ignored results, and checks the hash of every figure/source table in each PlotResult manifest. The optional public-data lesson records non-execution if no files are supplied. These tutorials do not modify the frozen workflow implementations or configs.

## Reading a source notebook

The ordered cells describe what the workflow expects and why a result has its
stated interpretation. Tracked notebooks retain stable code-cell source and
order, with null execution counts and no embedded outputs. Read them alongside
the linked scripts and the [results guide](RESULTS_OVERVIEW.md).

Source notebooks do not serve as a complete rendered archive of all manuscript
results. Their generated inputs and outputs may need to be recovered from the
recorded environment or produced by earlier stages.

## Running and inspecting outputs

The corresponding executor starts a clean kernel and records execution status.
Where supported, request a separate executed notebook copy; for example,
`execute_revision_notebooks.py` uses `--keep-executed`. Generated review copies
are written under ignored result paths so source notebooks remain suitable for
version control.

Figure-producing stages can write PNG/SVG images, source CSVs, alt text, and
metadata together. The source table and configuration explain what the image
summarizes; the version and checksum records identify the inputs used to create
it. Consult the dataset README for the specific output locations.

## Adapting a workflow

Change reusable numerical logic in its script and document any altered design,
parameters, or inputs. Validate the affected results before using them in a new
study. Keep presentation-only edits distinct from numerical changes when
rebuilding a manuscript figure from frozen evidence.

For an editorial pass, the readability validator compares source notebooks with
a chosen base revision and checks cell types, code source/order, outputs, and
execution counts. A numerical implementation change requires its own validation;
passing a readability check does not establish scientific equivalence.
