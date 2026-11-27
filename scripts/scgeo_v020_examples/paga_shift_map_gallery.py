"""Render a branched, fully synthetic state-graph and displacement gallery.

Run from the development checkout containing the upgraded plotting API::

    python -m pip install -e .
    python examples/paga_shift_map_gallery.py --output-dir results/paga_shift_map

All coordinates, cell states, sample identities, velocities, and graph edges
are generated here. ``X_display`` is a declared linear projection, not UMAP.
The symmetric graph is a constructed fixture, not an inferred PAGA result.
ScGeo computes the stored quantitative evidence from the synthetic 3D data.
No downloads, Scanpy, scVelo, or external data are required.
"""

from __future__ import annotations

import argparse
import inspect
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

import scgeo as sg

STATES = [
    "Reservoir",
    "Branch point",
    "Upper early",
    "Upper late",
    "Lower early",
    "Lower late",
    "Side state",
    "Treatment-only",
]
STATE_PALETTE = {
    "Reservoir": "#356890",
    "Branch point": "#789DB5",
    "Upper early": "#1C8A85",
    "Upper late": "#69AA83",
    "Lower early": "#B06C8C",
    "Lower late": "#7B6C9D",
    "Side state": "#BF9252",
    "Treatment-only": "#7D8990",
}
CONDITION_PALETTE = {"control": "#387EA3", "treated": "#E19954"}
GRAPH_EDGES = [
    ("Reservoir", "Branch point", 0.92),
    ("Branch point", "Upper early", 0.80),
    ("Upper early", "Upper late", 0.87),
    ("Branch point", "Lower early", 0.83),
    ("Lower early", "Lower late", 0.90),
    ("Reservoir", "Side state", 0.42),
    ("Lower late", "Treatment-only", 0.50),
]


def make_paga_demo(seed=41):
    """Build eight anisotropic state clouds with four samples per condition.

    A third, unrelated condition has displaced coordinates and deliberately
    different velocities. It tests the distinction between a display background
    and the control/treated data used for displacement, pies, and velocity.
    The final state has no control cells and no estimable two-condition shift.
    """
    rng = np.random.default_rng(seed)
    centers = np.array(
        [
            [-4.5, 0.0, 0.05],
            [-2.2, 0.05, 0.15],
            [0.15, 1.55, 0.55],
            [2.55, 2.45, 0.95],
            [0.05, -1.50, -0.35],
            [2.60, -2.35, -0.80],
            [-2.9, 2.65, 0.35],
            [4.8, -1.15, -0.15],
        ]
    )
    effects = np.array(
        [
            [0.20, 0.08, 0.03],
            [0.42, 0.12, 0.18],
            [0.56, 0.32, 0.10],
            [0.62, 0.17, 0.05],
            [-0.30, 0.45, 0.25],
            [0.22, 0.55, 0.12],
            [0.0, 0.0, 0.0],
            [-0.18, 0.32, 0.10],
        ]
    )
    directions = np.array(
        [
            [1.0, 0.1, 0.08],
            [0.8, 0.25, 0.3],
            [0.7, 0.7, 0.1],
            [1.0, 0.15, 0.1],
            [0.6, -0.7, -0.1],
            [1.0, -0.38, 0.0],
            [-0.4, 0.8, 0.1],
            [0.6, 0.8, 0.1],
        ]
    )
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    counts = {
        "control": [45, 50, 40, 35, 42, 38, 20, 0],
        "treated": [36, 44, 48, 52, 35, 30, 18, 12],
        "context": [6] * 8,
    }
    angles = np.array([0.0, 0.15, 0.65, 0.3, -0.55, -0.2, 1.0, 0.45])
    records, coordinates, velocity = [], [], []
    for condition in ("control", "treated", "context"):
        strength = 1.0 if condition == "treated" else 0.0
        for sample_number in range(4):
            # Sample identities are independent across conditions, not paired.
            sample_offset = rng.normal(0, 0.075, 3)
            for state_number, state in enumerate(STATES):
                n_cells = counts[condition][state_number]
                if not n_cells:
                    continue
                angle = angles[state_number]
                rotation_2d = np.array(
                    [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
                )
                noise_xy = rng.normal(size=(n_cells, 2)) * [0.34, 0.20]
                noise_xy = noise_xy @ rotation_2d.T
                noise = np.column_stack([noise_xy, rng.normal(0, 0.22, n_cells)])
                state_sample_offset = rng.normal(0, 0.035, 3)
                center = centers[state_number] + effects[state_number] * strength
                if condition == "context":
                    center = center + [0.0, 8.5, 1.4]
                points = center + sample_offset + state_sample_offset + noise
                if condition == "context":
                    vectors = np.tile([0.8, -0.6, 0.3], (n_cells, 1))
                else:
                    speed = 0.025 * (1.08 if condition == "treated" else 1.0)
                    vectors = speed * directions[state_number] + rng.normal(
                        0, 0.0015, (n_cells, 3)
                    )
                coordinates.append(points)
                velocity.append(vectors)
                records.extend(
                    {
                        "state": state,
                        "condition": condition,
                        "sample": f"{condition}_sample_{sample_number + 1}",
                    }
                    for _ in range(n_cells)
                )
    obs = pd.DataFrame(
        records, index=[f"synthetic_cell_{i:04d}" for i in range(len(records))]
    )
    obs["state"] = pd.Categorical(obs["state"], categories=STATES, ordered=True)
    obs["condition"] = pd.Categorical(
        obs["condition"], categories=["control", "treated", "context"]
    )
    adata = AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
    adata.var_names = ["unused_placeholder"]
    primary = np.vstack(coordinates)
    vectors = np.vstack(velocity)
    rotation_3d = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    adata.obsm["X_synthetic"] = primary
    adata.obsm["X_rotated"] = primary @ rotation_3d
    adata.obsm["V_synthetic"] = vectors
    adata.obsm["V_rotated"] = vectors @ rotation_3d
    adata.obsm["X_display"] = primary[:, :2].copy()
    adata.obsm["velocity_display"] = vectors[:, :2].copy()
    graph = np.zeros((len(STATES), len(STATES)))
    for left, right, weight in GRAPH_EDGES:
        i, j = STATES.index(left), STATES.index(right)
        graph[i, j] = graph[j, i] = weight
    adata.uns["paga_fixture"] = {
        "connectivities": graph,
        "groups": "state",
        "node_order": list(STATES),
        "source": "Constructed symmetric topology fixture; no PAGA inference was run.",
    }
    adata.uns["synthetic_demo"] = {
        "seed": int(seed),
        "samples_per_condition": 4,
        "sample_design": "Independent sample identities across conditions; no pairing.",
        "display": "First two coordinates of X_synthetic; not a UMAP embedding.",
        "topology": "Constructed symmetric graph fixture; not inferred PAGA connectivity.",
        "injected_shifts": {
            state: effects[i].tolist() for i, state in enumerate(STATES)
        },
        "cells_per_sample": counts,
    }
    return adata


def compute_demo_evidence(adata, seed=41):
    """Compute actual sample-aware ScGeo results on the synthetic coordinates."""
    common = dict(
        condition_key="condition",
        group0="control",
        group1="treated",
        sample_key="sample",
        center="geometric_median",
        seed=seed,
    )
    sg.tl.robust_shift(adata, rep="X_synthetic", by="state", n_boot=80, **common)
    sg.tl.representation_stability(
        adata,
        reps=["X_synthetic", "X_rotated"],
        node_key="state",
        velocity_keys={"X_synthetic": "V_synthetic", "X_rotated": "V_rotated"},
        n_boot=20,
        min_cells=15,
        **common,
    )
    return adata


def _shared_display_limits(adata):
    mask = adata.obs["condition"].isin(["control", "treated"]).to_numpy()
    coordinates = adata.obsm["X_display"][mask]
    low, high = coordinates.min(axis=0) - 0.55, coordinates.max(axis=0) + 0.55
    return tuple(zip(low, high))


def run_gallery(output_dir, formats=("png", "svg", "pdf"), seed=41):
    """Save four views with exact plotted tables and interpretation metadata."""
    parameters = inspect.signature(sg.pl.paga_shift_map).parameters
    if not {"graph_node_order", "return_result"} <= set(parameters):
        raise RuntimeError(
            "This gallery requires the development checkout with the upgraded "
            "paga_shift_map API. Install that checkout with `python -m pip install -e .` "
            "using the same Python interpreter that runs this example."
        )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = compute_demo_evidence(make_paga_demo(seed), seed)
    xlim, ylim = _shared_display_limits(adata)
    provenance = {
        "scope": "Fully synthetic API demonstration; no biological findings.",
        "simulation": adata.uns["synthetic_demo"],
        "quantitative_representation": "X_synthetic",
        "representation_roles": {
            "X_synthetic": "primary synthetic 3D representation",
            "X_rotated": "related orthogonal sensitivity view, not independent validation",
        },
        "graph_source": adata.uns["paga_fixture"]["source"],
        "graph_matrix_row_order": list(STATES),
        "display_centers": "Cell-level mean centers in X_display, not multivariate sample-aware robust centers.",
        "node_size": "Observed cell counts; not numbers of biological samples.",
        "velocity_display": "When shown: mean projected velocity across control and treated cells, multiplied by 20 for display.",
        "condition_scope": ["control", "treated"],
        "excluded_context": "context cells are excluded from contrast, pies, and velocity summaries.",
    }
    common = dict(
        node_key="state",
        condition_key="condition",
        group0="control",
        group1="treated",
        basis="display",
        paga_key="paga_fixture",
        graph_node_order=STATES,
        min_cells=15,
        agg="mean",
        palette=STATE_PALETTE,
        background="contrast",
        max_points=1500,
        random_state=seed,
        xlim=xlim,
        ylim=ylim,
        aspect="equal",
        connectivity_threshold=0.10,
        background_size=5,
        background_alpha=0.23,
        node_size=175,
        endpoint_size=48,
        edge_lw=1.8,
        edge_alpha=0.32,
        arrow_linewidth=1.7,
        arrow_head_size=10,
        arrow_scale=1.0,
        show_endpoints=True,
        label_fontsize=9,
        label_mode="auto",
        legend_loc="outside",
        figsize=(10.5, 7.0),
        return_result=True,
        show=False,
    )

    def save(stem, result, caption):
        result.metadata["synthetic_fixture"] = provenance
        result.metadata["gallery_caption"] = caption
        result.fig.text(
            0.055,
            0.032,
            "SYNTHETIC DISPLAY + CONSTRUCTED GRAPH",
            fontsize=8,
            weight="bold",
            color="#627584",
        )
        result.fig.text(0.055, 0.012, caption, fontsize=7, color="#627584")
        exported = result.save(
            output_dir / stem, formats=formats, dpi=220, bbox_inches="tight"
        )
        plt.close(result.fig)
        return exported

    exports = {}
    with sg.pl.style_context("paper"):
        overview = sg.pl.paga_shift_map(
            adata,
            **common,
            show_velocity=True,
            velocity_basis="display",
            velocity_scale=20,
            velocity_color="#178C86",
            title="State graph and condition displacement",
        )
        exports["overview"] = save(
            "paga_shift_map_overview",
            overview,
            "Control → treated centers; velocity ×20. One state has no control cells, so no displacement is drawn.",
        )

        aligned = sg.pl.paga_shift_map(
            adata,
            **common,
            node_color_mode="alignment",
            alignment_source="stored",
            analysis_rep="X_synthetic",
            cmap="coolwarm",
            vmin=-1,
            vmax=1,
            title="Stored alignment in the primary representation",
        )
        exports["alignment"] = save(
            "paga_shift_map_alignment",
            aligned,
            "Colors read stored 3D alignment; arrows show display-space displacement. The triangle marks a state without a reference condition.",
        )

        composition = sg.pl.paga_shift_map(
            adata,
            **common,
            node_color_mode="alignment",
            alignment_source="stored",
            analysis_rep="X_synthetic",
            cmap="coolwarm",
            vmin=-1,
            vmax=1,
            pie_key="condition",
            pie_categories=["control", "treated"],
            pie_palette=CONDITION_PALETTE,
            pie_size_scale=1.35,
            title="Condition composition and stored alignment",
        )
        exports["composition"] = save(
            "paga_shift_map_composition",
            composition,
            "Sectors: observed control/treated cell fractions. Rim: stored 3D alignment. No sample-level uncertainty is estimated.",
        )

        fig, ax = plt.subplots(figsize=(10.5, 7.0))
        # Allocate the final axes geometry before automatic label placement.
        fig.subplots_adjust(left=0.08, right=0.76, bottom=0.16, top=0.88)
        focused = sg.pl.paga_shift_map(
            adata,
            **common,
            ax=ax,
            label_top_n=3,
            highlight_nodes=["Upper late", "Lower early", "Treatment-only"],
            highlight_edgecolor="#182F40",
            highlight_lw=2.4,
            title="Focus selected states on the same coordinate limits",
        )
        exports["focus"] = save(
            "paga_shift_map_focus",
            focused,
            "Caller-supplied axis and shared limits; top labels follow display displacement. Highlights are chosen explicitly.",
        )

    print(
        f"Wrote four synthetic state-graph views, source tables, and metadata to {output_dir.resolve()}"
    )
    return exports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/paga_shift_map")
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=["png", "svg", "pdf"],
        default=["png", "svg", "pdf"],
    )
    parser.add_argument("--seed", type=int, default=41)
    args = parser.parse_args()
    run_gallery(args.output_dir, formats=tuple(args.formats), seed=args.seed)


if __name__ == "__main__":
    main()
