"""Generate ScGeo's plotting gallery from small, fully synthetic inputs.

Run from an installed development checkout::

    python examples/plotting_gallery.py --output-dir results/plotting_gallery

The three conditions and their sample identities are simulated. Related
rotated coordinates illustrate representation sensitivity, not independent
validation. No external data are downloaded and no manuscript outputs are read.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

import scgeo as sg


def make_gallery_data(seed=17):
    rng = np.random.default_rng(seed)
    centers = {"State A": [0, 0, 0], "State B": [3, 0, 0], "State C": [0, 3, 0]}
    effects = {"State A": [1.1, .3, .2], "State B": [.4, -.1, .15], "State C": [0, 0, 0]}
    coordinates, velocities, records = [], [], []
    for condition, strength in (("control", 0), ("mild", .45), ("treated", 1)):
        for sample in range(4):
            offset = rng.normal(0, .12, 3)
            for state, center in centers.items():
                n = 20
                delta = np.asarray(effects[state]) * strength
                coordinates.append(center + offset + delta + rng.normal(0, .2, (n, 3)))
                direction = np.array([1, .2, .1]) if state != "State B" else np.array([-1, 0, 0])
                velocities.append(np.tile(direction, (n, 1)) * .02)
                records.extend({"state": state, "condition": condition, "sample": f"{condition}_{sample}"}
                               for _ in range(n))
    obs = pd.DataFrame(records, index=[f"cell_{i}" for i in range(len(records))])
    adata = AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
    coords, velocity = np.vstack(coordinates), np.vstack(velocities)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    adata.obsm["X_toy"] = coords
    adata.obsm["X_rotated"] = coords @ rotation
    adata.obsm["V_toy"] = velocity
    adata.obsm["V_rotated"] = velocity @ rotation
    # A declared linear display projection, not an inferred UMAP trajectory.
    adata.obsm["X_display"] = coords[:, :2].copy()
    adata.obsm["velocity_display"] = velocity[:, :2].copy()
    return adata


def run_gallery(output_dir, formats=("png", "svg")):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = make_gallery_data()
    states = ["State A", "State B", "State C"]
    source_metadata = {
        "scope": "Synthetic API demonstration; no biological findings.",
        "seed": 17, "cells": adata.n_obs, "samples_per_condition": 4,
        "quantitative_representation": "X_toy",
        "representation_roles": {"X_toy": "primary", "X_rotated": "related rotated sensitivity view"},
        "display": "First two coordinates of X_toy; linear projection only.",
    }

    def save(name, result, *, provenance=None):
        result.metadata.update(source_metadata if provenance is None else provenance)
        exported = result.save(output_dir / name, formats=formats, dpi=180, bbox_inches="tight")
        plt.close(result.fig)
        return exported

    shift_kwargs = dict(rep="X_toy", condition_key="condition", group0="control",
                        by="state", sample_key="sample", center="geometric_median", n_boot=40, seed=17)
    sg.tl.robust_shift(adata, group1="treated", **shift_kwargs)
    sg.tl.representation_stability(
        adata, reps=["X_toy", "X_rotated"], node_key="state", condition_key="condition",
        group0="control", group1="treated", sample_key="sample", center="geometric_median",
        n_boot=12, min_cells=10, seed=17,
    )
    report = sg.get.state_report(adata, node_key="state")
    # This first report intentionally records local geometry and dynamics as
    # not computed. Subsequent gallery examples add those analyses explicitly.
    with sg.pl.style_context("paper"):
        save("state_evidence", sg.pl.state_evidence_panel(
            report, effect_metric="delta_norm", representation="X_toy", state_order=states,
            evidence_mode="full", show_sample_counts=True, show_inference=True,
            figsize=(13, 3.8), title="Synthetic evidence overview · treated versus control",
            return_result=True, show=False,
        ))
        sample_table = sg.get.sample_centers(adata)
        save("biological_replicates", sg.pl.replicate_effect(
            sample_table, coordinate=0, state_order=states, condition_order=["control", "treated"],
            title="Synthetic samples · one point per biological unit", return_result=True, show=False,
        ))

        comparisons = []
        for treatment in ("mild", "treated"):
            copy = adata.copy()
            sg.tl.robust_shift(copy, group1=treatment, **shift_kwargs)
            current = sg.get.state_report(copy, node_key="state")
            comparison = current[["state", "delta_norm", "magnitude_ci95_low", "magnitude_ci95_high"]].rename(
                columns={"delta_norm": "estimate", "magnitude_ci95_low": "interval_low",
                         "magnitude_ci95_high": "interval_high"})
            comparison["contrast"] = f"{treatment} vs control"
            comparison["metric"] = "raw displacement magnitude"
            comparison["interval_type"] = "sample_bootstrap_ci95"
            comparison["representation"] = "X_toy"
            comparison["unit"] = "toy coordinate units"
            comparison["scale_id"] = "same_simulated_coordinate_system"
            comparisons.append(comparison)
        comparison_table = pd.concat(comparisons, ignore_index=True)
        comparison_table.attrs = dict(source_metadata)
        save("contrast_comparison", sg.pl.evidence_comparison(
            comparison_table, state_order=states, title="Synthetic contrasts · common coordinate system",
            return_result=True, show=False,
        ))

        sg.tl.representation_stability(
            adata, reps=["X_toy", "X_rotated"], node_key="state", condition_key="condition",
            group0="control", group1="treated", sample_key="sample", center="geometric_median",
            velocity_keys={"X_toy": "V_toy", "X_rotated": "V_rotated"},
            n_boot=12, min_cells=10, seed=17,
        )
        save("representation_alignment", sg.pl.representation_stability_heatmap(
            adata, metric="alignment_cosine", state_order=states, rep_order=["X_toy", "X_rotated"],
            title="Stored alignment · related representation views", return_result=True, show=False,
        ))
        save("consensus_map", sg.pl.consensus_state_map(
            adata, node_key="state", basis="display", max_points=300, random_state=17,
            title="Stored consensus on a linear display projection", return_result=True, show=False,
        ))
        save("geometry_dynamics", sg.pl.velocity_shift_alignment(
            adata, node_key="state", condition_key="condition", group0="control", group1="treated",
            basis="display", mode="stored", analysis_rep="X_toy", min_cells=10,
            figsize=(10, 4.5), bg_alpha=.35, velocity_color="#16868b", shift_color="#34495e",
            title="Stored quantitative alignment and display overlay", return_result=True, show=False,
        ))
        save("composition", sg.pl.composition_drift(
            adata, node_key="state", condition_key="condition", group0="control", group1="treated",
            basis="display", state_order=states,
            title="Synthetic composition · balanced cell counts by construction", return_result=True, show=False,
        ))
        sg.tl.local_geometry_stability(
            adata, reps=["X_toy", "X_rotated"], node_key="state", sample_key="sample",
            k_values=(5,), n_boot=3, max_exact_cells=800, store_per_cell=True, seed=17,
        )
        save("local_geometry", sg.pl.local_distortion_map(
            adata, basis="display", node_key="state", max_points=300, random_state=17,
            title="Local distortion · related rotated coordinates", return_result=True, show=False,
        ))
        # Missingness fixture is deliberately separate from measured geometry.
        matrix = pd.DataFrame({"state": ["State A", "State A", "State A", "State B", "State B", "State B"],
                               "setting": ["one", "two", "three"] * 2,
                               "value": [.82, np.nan, .91, .61, .82, .73]})
        save("robustness_coverage", sg.pl.robustness_matrix(
            matrix, row_key="state", col_key="setting", value_key="value", summary="pass_rate", pass_threshold=.7,
            col_order=["one", "two", "three"], figsize=(10, 3.5),
            title="Diagnostic fixture · robustness and coverage", return_result=True, show=False,
        ), provenance={
            "scope": "Constructed diagnostic missingness fixture; not derived from cell measurements.",
            "settings": 3, "states": 2,
        })
    print(f"Wrote synthetic gallery figures, source tables, and metadata to {output_dir.resolve()}")
    return output_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/plotting_gallery"))
    parser.add_argument("--formats", nargs="+", choices=["png", "svg", "pdf"], default=["png", "svg"])
    args = parser.parse_args()
    run_gallery(args.output_dir, formats=tuple(args.formats))


if __name__ == "__main__":
    main()
