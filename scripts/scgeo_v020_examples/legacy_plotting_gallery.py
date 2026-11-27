"""Exercise upgraded legacy plotting APIs on the existing synthetic generator.

No public/private biological data or downloads are used. Stored summaries are
calculated by unchanged ScGeo methods on explicitly synthetic coordinates.
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
import scgeo as sg
from plotting_gallery import make_gallery_data


def run(output_dir, formats=("png", "svg", "pdf")):
    output_dir = Path(output_dir)
    adata = make_gallery_data()
    adata.obs["state"] = pd.Categorical(adata.obs["state"], categories=["State A", "State B", "State C"])
    adata.uns["state_colors"] = ["#427CB3", "#D88843", "#4E9B85"]
    adata.obs["illustrative_score"] = np.linalg.norm(adata.obsm["X_toy"], axis=1)
    rng = np.random.default_rng(17)
    adata.obs["map_confidence"] = rng.uniform(.1, 1., adata.n_obs)
    adata.obs["map_ood_score"] = adata.obs["illustrative_score"]
    predictions = adata.obs["state"].astype(object).to_numpy().copy()
    query_positions = np.flatnonzero(adata.obs["condition"].eq("treated"))
    predictions[query_positions[:20]] = "Novel label"
    predictions[query_positions[20:30]] = None
    adata.obs["map_pred"] = pd.Categorical(predictions)
    common = dict(condition_key="condition", group0="control", group1="treated", by="state")
    sg.tl.shift(adata, rep="X_toy", sample_key="sample", **common)
    sg.tl.density_overlap(adata, rep="X_display", **common)
    sg.tl.distribution_test(adata, rep="X_toy", sample_key="sample", n_perm=19, seed=17, **common)
    outputs = {}

    def save(name, result):
        result.metadata["gallery"] = {"dataset": "fully_synthetic", "seed": 17,
            "display": "linear first-two-coordinate projection; not fitted UMAP",
            "inference": "unchanged ScGeo methods on synthetic data; 19 permutations for renderer demonstration"}
        outputs[name] = result.save(output_dir / name, formats=formats, dpi=150)
        plt.close(result.fig)

    with sg.pl.style_context("paper"):
        save("distribution_report", sg.pl.distribution_report(
            adata, basis="display", group0="control", group1="treated",
            score_key="illustrative_score", return_result=True, show=False))
        fig, axs = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
        fig.suptitle("Synthetic stored summaries · exact rows and ordering")
        overlap = sg.pl.density_overlap(adata, axes={"bc": axs[0], "hellinger": axs[1]},
                                        return_result=True, show=False)
        ranking = sg.pl.delta_rank(adata, ax=axs[2], return_result=True)
        save("overlap_and_rank", sg.pl.PlotResult(
            fig, {**overlap.axes, **ranking.axes}, data={"overlap": overlap.data, "ranking": ranking.data},
            metadata={"overlap": overlap.metadata, "ranking": ranking.metadata}))
        save("separate_distribution_metrics", sg.pl.distribution_test(
            adata, level="global", separate_metrics=True, return_result=True, show=False))
        qc = sg.pl.mapping_qc_panel(adata, basis="display", palette_from="state",
            condition_key="condition", query_value="treated", return_result=True, show=False)
        qc.fig.suptitle("Synthetic mapping display · illustrative scores; no mapping model fitted")
        save("mapping_legend", sg.pl.legend_from_data(qc.data["legend"].to_dict(orient="records"),
            title="Query predictions · synthetic", return_result=True, show=False))
        save("mapping_qc", qc)
        save("ood_landscape", sg.pl.ood_landscape(adata, basis="display", ood_key="map_ood_score",
            groupby="state", contour_gridsize=12, contour_quantile=.75,
            title="Synthetic OOD display · descriptive score threshold",
            return_result=True, show=False))
        fig, axs = plt.subplots(1, 3, figsize=(10, 4.5), gridspec_kw={"width_ratios": [4, .2, 2]}, constrained_layout=True)
        structure = sg.pl.consensus_structure(adata, basis="display", groupby="state",
            score_key="illustrative_score", topk=150, ax=axs[0], cax=axs[1],
            title="Synthetic top-score structure", print_legend=False, return_result=True, show=False)
        legend = sg.pl.legend_from_data(structure.data["legend"].to_dict(orient="records"), ax=axs[2],
            title="Selected cells", return_result=True, show=False)
        save("consensus_structure", sg.pl.PlotResult(fig, {**structure.axes, **legend.axes},
            data={"structure": structure.data, "legend": legend.data},
            metadata={"structure": structure.metadata, "legend": legend.metadata}))
        # Constructed undirected topology for renderer demonstration, not inferred PAGA.
        adata.uns["paga"] = {"groups": "state", "node_order": ["State A", "State B", "State C"],
            "connectivities": np.array([[0., .8, .5], [.8, 0., 0.], [.5, 0., 0.]])}
        compass = sg.pl.recovery_compass(adata, node_key="state", condition_key="condition",
            group0="control", group1="treated", basis="display", ood_key="map_ood_score",
            bg_alpha=.25,
            title="Synthetic display-space compass · constructed graph",
            return_result=True, show=False)
        compass.metadata["graph_source"] = "constructed synthetic topology; not inferred PAGA"
        save("recovery_compass", compass)
    return outputs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/legacy_plotting"))
    parser.add_argument("--formats", nargs="+", choices=("png", "svg", "pdf"), default=["png", "svg", "pdf"])
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir, tuple(args.formats)), indent=2))
