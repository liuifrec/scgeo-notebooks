# Reproducing and adapting ScGeo workflows

A reproducible result needs its source code, numerical inputs, configuration,
and execution environment. This guide explains how those pieces are recorded
in the companion and where a fresh workstation needs additional setup.

[Start here](START_HERE.md) · [Workflow catalog](DATASET_MAP.md) ·
[Results guide](RESULTS_OVERVIEW.md) · [Citation and reporting](PUBLICATION_GUIDE.md)

## Choose the appropriate starting point

| Goal | Package and workflow choice |
|---|---|
| Learn the API or analyze a new dataset | Use the published package and [current quick start](https://github.com/liuifrec/scgeo/blob/main/docs/QUICKSTART.md); record the version and your analysis choices |
| Explore the report format | Run the [small synthetic tutorial](START_HERE.md#try-a-small-report-example), which creates its own data |
| Reproduce frozen manuscript evidence | Use the checkpoints, recorded environment, and checksum-pinned inputs below |
| Adapt a manuscript workflow | Start from a recorded revision, document your changes, and validate the new inputs/design |

The environment YAMLs under [`environment/`](../environment/) are
workflow-specific specifications, not a single lockfile for all analyses.
Some stages use separately recorded Python/R environments or scVI
interpreters. Consult the dataset README and execution metadata before
assuming that one YAML is sufficient for the complete pipeline.

## Frozen repositories and numerical inputs

| Component | Checkpoint |
|---|---|
| ScGeo package | [`9a0ed16cbaa57f935f9c9bc87d1643a25b51012c`](https://github.com/liuifrec/scgeo/tree/9a0ed16cbaa57f935f9c9bc87d1643a25b51012c) |
| Companion merge checkpoint | [`958a0ada568c53e78864ede022c050830cb55a36`](https://github.com/liuifrec/scgeo-notebooks/tree/958a0ada568c53e78864ede022c050830cb55a36) |

The package pin is also recorded in the workflow configs. Dataset-specific
checksums and execution metadata identify the numerical objects used by later
stages. A matching filename is not a substitute for a matching checksum.
Historical assembly and finalization scripts may impose additional exact-path
or clean-worktree checks; inspect those requirements before running them.

For a frozen reproduction, record both Git commits, confirm a clean worktree,
and preserve the original environment specification or package-version record.
If a documented artifact is unavailable, report the missing input rather than
silently substituting another run.

## Execution sequence

1. Open the relevant [dataset README](DATASET_MAP.md#workflow-entry-points) and
   identify its numbered stages and expected input objects.
2. Create the recorded environment(s), select the correct Jupyter kernel, and
   confirm that the intended package checkout is imported.
3. Acquire or locate the required data, then validate its metadata and checksums.
4. Configure supported paths and run the matching executor from the repository
   root. Use the executor's `--help` for its actual options.
5. Inspect the execution report and stage outputs before proceeding to
   representation construction, ScGeo analysis, or figure assembly.
6. Preserve output checksums and version records with the resulting figure or
   table bundle.

Main executors include
[`execute_revision_notebooks.py`](../scripts/execute_revision_notebooks.py),
[`execute_pancreas_validation.py`](../scripts/execute_pancreas_validation.py),
[`execute_gse249479_validation.py`](../scripts/execute_gse249479_validation.py),
and [`execute_gse211713_validation.py`](../scripts/execute_gse211713_validation.py).
The dataset guides identify the additional staged executors. These commands
are not interchangeable, and their kernel option names differ.

## Supported path overrides

Set only variables read by the script you are invoking. Dataset-specific
executors record resolved inputs and checksums; an override supported by one
stage may not replace a fixed path in another stage.

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

## Portability of the frozen workflows

Some validated executors and gates retain absolute paths under
`/home/liuyuchen`. Examples include GSE211713 v2 validation/execution,
[`finalize_revision_evidence.py`](../scripts/finalize_revision_evidence.py),
the recorded scVI interpreter, and GSE249479 execution/comparator wrappers.
These paths reflect the validated workstation environment. A fresh clone and
package installation alone are therefore insufficient for every frozen stage.

For another machine, inspect the paths in the selected workflow and distinguish
supported overrides from implementation changes. Any change to a frozen path
assumption, preprocessing stage, numerical configuration, or package version
should be recorded and checked for numerical equivalence before being described
as an exact reproduction.

## Source notebooks and generated artifacts

Tracked source notebooks contain no outputs and have null execution counts.
Executors can write separate clean-kernel review copies under ignored results
paths. Figures, source CSVs, alt text, metadata, and logs accompany the relevant
workflow. H5ADs, models, downloads, and caches remain outside the tracked source
unless policy explicitly retains a small manifest or documentation artifact.

A figure bundle is most useful when its image travels with the source table,
caption/alt text, configuration, and version record. See
[Thin notebook architecture](THIN_NOTEBOOK_ARCHITECTURE.md) and
[Citation and reporting](PUBLICATION_GUIDE.md#report-the-analysis-you-actually-ran).

## Verification checklist

- Verify both repository commits and worktree cleanliness.
- Parse the JSON/YAML configurations and inspect the selected kernel/interpreter.
- Check that source notebooks remain output-free and their code is unchanged.
- Verify numerical input/output checksums and reject invalidated-run inputs.
- Confirm unique canonical result keys before assembling tables or figures.
- Run `git diff --check` and the validation appropriate to the changed component.

The existing `scripts/validate_readability_only.py --base <revision>` checks a
historical documentation-only contract. Its allowlist is intentionally narrow;
it does not replace validation of a new feature, citation format, or numerical
implementation change.
