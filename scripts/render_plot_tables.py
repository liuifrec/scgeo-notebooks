"""Render existing result CSVs with the development ScGeo plotting API.

This presentation-only entry point reads a source table, records its hash, and
exports figures and source data to a separate destination. It does not run an
analysis or change the manuscript's frozen numerical outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from io import BytesIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scgeo as sg


def render_table(kind, input_path, output_stem, *, metric=None, coordinate=None,
                 state_order=None, formats=("pdf", "svg", "png"), title=None,
                 metadata=None, dpi=300):
    """Render one table and return the figure/source/provenance export mapping."""
    if not hasattr(sg.pl, "export_figure"):
        raise RuntimeError(
            "This renderer needs the ScGeo development checkout containing the plotting upgrades. "
            "Install that checkout with python -m pip install -e /path/to/scgeo; "
            "the released scgeo==0.1.0 does not contain these APIs."
        )
    input_path = Path(input_path)
    if Path(output_stem).resolve().parent == input_path.resolve().parent:
        raise ValueError("Choose an output directory separate from the source table directory.")
    source_bytes = input_path.read_bytes()
    columns = pd.read_csv(BytesIO(source_bytes), nrows=0).columns
    label_columns = {name: str for name in ("state", "node", "contrast", "metric", "representation",
                     "rep", "condition", "sample", "interval_type", "unit", "scale_id") if name in columns}
    missing_tokens = ["", "NA", "N/A", "NaN", "nan", "NULL", "null", "None", "<NA>"]
    identity_columns = {"state", "node", "contrast", "representation", "rep", "condition", "sample"}
    na_values = {name: ([""] if name in identity_columns else missing_tokens) for name in columns}
    table = pd.read_csv(BytesIO(source_bytes), dtype=label_columns,
                        keep_default_na=False, na_values=na_values)
    if metadata is not None and not isinstance(metadata, dict):
        raise ValueError("metadata must be a JSON object.")
    source = dict(metadata or {})
    source.update({
        "source_table": input_path.name,
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "presentation_scope": "Render stored table values; no analysis or inference recomputed.",
    })
    table.attrs.update(source)
    common = dict(state_order=state_order, title=title, return_result=True, show=False)
    with sg.pl.style_context("paper"):
        if kind == "state":
            result = sg.pl.state_evidence_panel(
                table, effect_metric=metric or "auto", show_sample_counts=True,
                show_inference=True, evidence_mode="full", max_states=None, **common,
            )
        elif kind == "effects":
            result = sg.pl.evidence_comparison(table, metric=metric, **common)
        elif kind == "replicates":
            if coordinate is None:
                raise ValueError("A replicate plot requires --coordinate (zero-based); it is never inferred.")
            result = sg.pl.replicate_effect(table, coordinate=coordinate, **common)
        elif kind == "representation":
            result = sg.pl.representation_stability_heatmap(
                table, metric=metric or "normalized_delta_norm", **common,
            )
        else:
            raise ValueError("kind must be state, effects, replicates, or representation.")
        try:
            result.metadata.update(source)
            return result.save(output_stem, formats=formats, dpi=dpi, bbox_inches="tight")
        finally:
            plt.close(result.fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["state", "effects", "replicates", "representation"])
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--metric")
    parser.add_argument("--coordinate", type=int)
    parser.add_argument("--state-order", nargs="+")
    parser.add_argument("--formats", nargs="+", choices=["pdf", "svg", "png"], default=["pdf", "svg", "png"])
    parser.add_argument("--title")
    parser.add_argument("--metadata", type=Path, help="JSON object with study/analysis provenance to retain")
    parser.add_argument("--dpi", type=float, default=300)
    args = parser.parse_args()
    metadata = json.loads(args.metadata.read_text(encoding="utf-8")) if args.metadata else None
    output = render_table(
        args.kind, args.input, args.output_stem, metric=args.metric, coordinate=args.coordinate,
        state_order=args.state_order, formats=tuple(args.formats), title=args.title,
        metadata=metadata, dpi=args.dpi,
    )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
