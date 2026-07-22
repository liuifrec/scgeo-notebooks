# Thin notebook architecture

## Why notebooks are thin

Revision notebooks are reviewable orchestration documents. They state purpose,
inference scope, frozen inputs, outputs, findings, and limitations, then call a
testable implementation under `scripts/`. This separates scientific narration
from reusable execution logic and keeps code review focused.

## Tracked source notebook

A source notebook is committed with:

- stable cell order and code-cell source;
- null execution counts;
- no embedded outputs;
- Markdown that explains the frozen workflow;
- calls into the corresponding script.

## Generated review copy

The executor starts a clean kernel, runs the source notebook, and writes the
executed copy beneath an ignored results directory. Figures, source CSVs, alt
text, logs, and metadata are generated alongside it. The executed copy is for
review and validation, not version control.

## Where to change logic

Numerical or scientific implementation changes belong in the corresponding
script and require their own validation task. Documentation passes may edit
Markdown only. A notebook-readability validator compares source notebooks with
the base branch and rejects code-source, code-order, cell-type, execution-count,
or output changes.
