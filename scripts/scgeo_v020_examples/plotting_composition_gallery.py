"""Compose existing ScGeo renderers on one caller-owned canvas; no downloads.

All coordinates, annotations and scores come from the existing fully synthetic
faceting gallery. They are not biological evidence or a fitted UMAP.
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
from embedding_facets_gallery import synthetic


def run(output_dir, formats=("png", "svg", "pdf")):
    output_dir = Path(output_dir)
    adata = synthetic()
    with sg.pl.style_context("paper"):
        fig = plt.figure(figsize=(15, 4.5), constrained_layout=True)
        grid = fig.add_gridspec(2, 4, height_ratios=[1, 0.10])
        atlas_axes = {"states": fig.add_subplot(grid[0, 0])}
        colorbar_ax = fig.add_subplot(grid[1, 0])
        facet_axes = {f"panel_{i:02d}": fig.add_subplot(grid[0, i + 1]) for i in range(3)}
        facet_axes["legend"] = fig.add_subplot(grid[1, 1:])
        fig.suptitle("Fully synthetic example · composed atlas and matched velocity facets")
        atlas = sg.pl.geometry_dynamics_atlas(
            adata, color="illustrative_score", panels=["states"],
            axes=atlas_axes, show_colorbar=False, show=False)
        atlas.axes["states"].set_title("All cells · illustrative score")
        # A caller can place a differently oriented colorbar using the returned
        # mappable, without redrawing points or changing its normalization.
        atlas.artists["states_colorbar"] = fig.colorbar(
            atlas.artists["states_points"], cax=colorbar_ax,
            orientation="horizontal", label="Illustrative score")
        facets = sg.pl.embedding_facets(
            adata, facet="condition", color="state", velocity="grid",
            velocity_options={"grid_size": 22, "n_neighbors": 24, "min_neighbors": 5},
            axes=facet_axes, context="all", return_result=True, show=False)
        # Each child result describes only its named panels. Bundle both sources
        # explicitly so the exported parent figure has complete provenance.
        result = sg.pl.PlotResult(
            fig=fig,
            axes={"atlas_states": atlas.axes["states"], **facets.axes},
            data={"atlas": atlas.data, "facets": facets.data},
            metadata={"plot": "plotting_composition_gallery", "dataset": "fully_synthetic",
                      "seed": 31, "atlas": atlas.metadata, "facets": facets.metadata},
            warnings=atlas.warnings + facets.warnings,
            artists={"atlas": atlas.artists, "facets": facets.artists})
        output = result.save(output_dir / "composed_atlas_facets", formats=formats, dpi=180)
        plt.close(fig)
        fig = plt.figure(figsize=(10, 10), constrained_layout=True)
        grid = fig.add_gridspec(3, 3, width_ratios=[1.2, 1, .04])
        fig.suptitle("Fully synthetic score panels · one scale, all-input histograms")
        results = {}
        values = adata.obs["illustrative_score"]
        norm = plt.Normalize(values.min(), values.max())
        for row, name in enumerate(("alignment_panel", "ambiguity_panel", "consensus_subspace_panel")):
            axes = dict(zip(("embedding", "histogram", "colorbar"),
                            [fig.add_subplot(grid[row, col]) for col in range(3)]))
            results[name] = getattr(sg.pl, name)(
                adata, score_key="illustrative_score", axes=axes, topk=20, bins=30,
                max_points=240, norm=norm, return_result=True, show=False)
            axes["embedding"].set_title(name.replace("_", " "))
        bundle = sg.pl.PlotResult(
            fig=fig,
            axes={f"{name}_{key}": ax for name, result in results.items() for key, ax in result.axes.items()},
            data={name: result.data for name, result in results.items()},
            metadata={"plot": "score_panel_gallery", "dataset": "fully_synthetic",
                      "scores": "illustrative random scores, not fitted alignment or consensus evidence",
                      "panels": {name: result.metadata for name, result in results.items()}})
        score_output = bundle.save(output_dir / "composed_score_panels", formats=formats, dpi=180)
        plt.close(fig)
        adata.uns["scgeo"] = {"paga_composition_stats": {
            "table": pd.DataFrame({"node": ["State A", "Missing effect", "Missing probability", "State D"],
                                   "effect": [1.4, np.nan, -.8, .1], "q": [.003, .02, np.nan, .8]}),
            "params": {"source": "explicitly synthetic stored values; no test fitted"}}}
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
        fig.suptitle("Synthetic composition table · missing evidence remains unavailable")
        table_result = sg.pl.paga_composition_panel(
            adata, axes=dict(zip(("volcano", "bar"), axes)), return_result=True, show=False)
        table_result.metadata["dataset"] = "fully_synthetic"
        table_output = table_result.save(output_dir / "composition_missing_evidence", formats=formats, dpi=180)
        plt.close(fig)
    return {"atlas_facets": output, "score_panels": score_output, "composition_evidence": table_output}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/plotting_composition"))
    parser.add_argument("--formats", nargs="+", choices=("png", "svg", "pdf"), default=["png", "svg", "pdf"])
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir, tuple(args.formats)), indent=2))
