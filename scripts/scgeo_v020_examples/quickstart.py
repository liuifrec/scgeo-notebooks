"""Run a small ScGeo analysis with synthetic coordinates and no downloads.

From a checkout, after installing ScGeo::

    python examples/quickstart.py --output-dir results/quickstart

The fixture illustrates the API, not a biological finding or benchmark. Each
condition contains four independent simulated samples and three cell states.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
import scgeo as sg


def make_toy_data(seed: int = 7) -> AnnData:
    """Create 720 cells in a single synthetic quantitative representation."""
    rng = np.random.default_rng(seed)
    centers = {"State A": [0, 0, 0], "State B": [3, 0, 0], "State C": [0, 3, 0]}
    shifts = {"State A": [1.1, 0.3, 0], "State B": [0.4, 0, 0], "State C": [0, 0, 0]}
    coordinates, rows = [], []
    for condition in ("control", "treated"):
        for sample in range(4):
            sample_offset = rng.normal(0, 0.12, 3)
            for state, center in centers.items():
                shift = np.asarray(shifts[state]) if condition == "treated" else np.zeros(3)
                coordinates.append(
                    np.asarray(center) + sample_offset + shift + rng.normal(0, 0.22, (30, 3))
                )
                rows.extend(
                    {"condition": condition, "sample": f"{condition}_{sample}", "cell_type": state}
                    for _ in range(30)
                )
    obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))])
    # X is a placeholder: this example analyzes obsm coordinates, not expression.
    adata = AnnData(X=np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
    adata.obsm["X_toy"] = np.vstack(coordinates)
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/quickstart"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    adata = make_toy_data()
    sg.tl.robust_shift(
        adata, rep="X_toy", condition_key="condition", group0="control",
        group1="treated", by="cell_type", sample_key="sample",
        center="geometric_median", n_boot=100, seed=7,
    )
    report = sg.get.state_report(adata, node_key="cell_type")
    report.to_csv(args.output_dir / "state_report.csv", index=False)
    with plt.rc_context({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False}):
        figure = sg.pl.state_evidence_panel(
            report, show_coverage=False, figsize=(7.5, 3.5),
            title="Synthetic example · four samples per condition", show=False,
        )
        for suffix in ("png", "svg"):
            figure.savefig(args.output_dir / f"state_evidence_panel.{suffix}", dpi=180, bbox_inches="tight")
        plt.close(figure)
    metadata = {
        "scgeo_version": sg.__version__, "seed": 7, "cells": adata.n_obs,
        "states": 3, "samples_per_condition": 4, "representation": "X_toy",
        "sample_key": "sample", "n_boot": 100, "bootstrap_unit": "sample",
        "scope": "Synthetic API demonstration; one representation; no biological conclusions.",
        "not_computed": ["representation_stability", "local_geometry", "dynamics"],
    }
    (args.output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    columns = ["state", "delta_norm", "n_samples0", "n_samples1", "bootstrap_unit"]
    print(report[columns].to_string(index=False))
    print(f"Wrote the report, figure and metadata to {args.output_dir.resolve()}")
    print("Representation stability, local geometry and dynamics were not computed.")


if __name__ == "__main__":
    main()
