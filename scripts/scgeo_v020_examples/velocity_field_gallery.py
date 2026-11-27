"""Render supported local fields on a fully synthetic branched display.

Run from the development checkout containing ``scgeo.pl.velocity_field``::

    python -m pip install -e .
    python examples/velocity_field_gallery.py --output-dir results/velocity_fields

Coordinates, graph edges, sample identities, and vectors are generated here.
They are not fitted UMAP coordinates, inferred PAGA topology, RNA velocity,
biological trajectories, or evidence of a treatment effect. No downloads,
Scanpy, or scVelo are required.
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
from matplotlib.colors import Normalize

import scgeo as sg

STATES = [
    "Source",
    "Trunk",
    "Junction",
    "Upper arm",
    "Upper island",
    "Lower bend",
    "Lower tip",
    "Rare pocket",
]
STATE_PALETTE = {
    "Source": "#496C8D",
    "Trunk": "#6D94B0",
    "Junction": "#6C9B9A",
    "Upper arm": "#359887",
    "Upper island": "#78A970",
    "Lower bend": "#A17496",
    "Lower tip": "#79699D",
    "Rare pocket": "#B89161",
}
GRID_EXTENT = (-5.65, 5.10, -3.45, 3.60)
FIELD_OPTIONS = {
    "grid_size": 39,
    "grid_extent": GRID_EXTENT,
    "n_neighbors": 45,
    "radius": 0.45,
    "min_neighbors": 12,
    "min_coherence": 0.30,
}


def _corridor(t, arm):
    """Return a constructed curve and its analytic tangent."""
    if arm == "trunk":
        x = -5.0 + 4.4 * t
        y = 0.23 * np.sin(np.pi * t) - 0.10
        dx = np.full_like(t, 4.4)
        dy = 0.23 * np.pi * np.cos(np.pi * t)
    else:
        sign = 1.0 if arm == "upper" else -1.0
        x = -0.6 + 5.0 * t
        y = sign * 2.45 * np.sin(np.pi * t / 2) - 0.05
        dx = np.full_like(t, 5.0)
        dy = sign * 2.45 * np.pi / 2 * np.cos(np.pi * t / 2)
    tangent = np.column_stack([dx, dy])
    tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
    return np.column_stack([x, y]), tangent


def make_velocity_demo(seed=71):
    """Generate curved corridors with a true gap and local counterflow.

    Four independent sample labels per condition are recorded for provenance.
    The plotting estimator pools eligible cells; it does not use those labels
    for biological-sample inference. The unrelated context condition has
    overlapping coordinates and large reversed vectors to expose scope errors.
    """
    rng = np.random.default_rng(seed)
    points, vectors, records = [], [], []
    for condition in ("control", "treated", "context"):
        treated = condition == "treated"
        for sample_number in range(4):
            sample_offset = rng.normal(0, 0.035, 2)
            for arm in ("trunk", "upper", "lower"):
                n = 145 if arm == "trunk" else 185
                if condition == "context":
                    n = 65
                t = rng.uniform(0, 1, n)
                curve, tangent = _corridor(t, arm)
                normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
                xy = curve + sample_offset
                xy += rng.normal(0, 0.08, (n, 1)) * tangent
                xy += rng.normal(0, 0.14, (n, 1)) * normal
                if treated:
                    xy += np.column_stack(
                        [np.full(n, 0.065), 0.10 + 0.055 * np.sin(np.pi * t)]
                    )
                # A coordinate-level exclusion leaves a real empty region even
                # after the sample offsets and condition displacements above.
                in_gap = (xy[:, 0] > 1.10) & (xy[:, 0] < 2.55) & (xy[:, 1] > 0.90)
                keep = ~in_gap
                xy, tangent, t = xy[keep], tangent[keep], t[keep]
                angle = (0.18 if treated else 0.0) * np.sin(np.pi * t) ** 2
                rotated = np.column_stack(
                    [
                        tangent[:, 0] * np.cos(angle) - tangent[:, 1] * np.sin(angle),
                        tangent[:, 0] * np.sin(angle) + tangent[:, 1] * np.cos(angle),
                    ]
                )
                speed = 0.047 + 0.014 * (1 + np.sin(2 * np.pi * t)) / 2
                speed *= 1.18 if treated else 1.0
                velocity = rotated * speed[:, None]
                opposing_patch = (arm == "lower") & (t > 0.38) & (t < 0.64)
                patch_positions = np.flatnonzero(opposing_patch)
                # Mixed forward/reverse vectors coexist in the same patch.
                reverse_positions = patch_positions[::2]
                velocity[reverse_positions] *= -1
                velocity += rng.normal(0, 0.0025, velocity.shape)
                if condition == "context":
                    velocity = -3.5 * rotated * speed[:, None]
                if arm == "trunk":
                    state = np.where(
                        t < 0.32, "Source", np.where(t < 0.72, "Trunk", "Junction")
                    )
                elif arm == "upper":
                    state = np.where(xy[:, 0] < 1.10, "Upper arm", "Upper island")
                else:
                    state = np.where(t < 0.55, "Lower bend", "Lower tip")
                points.append(xy)
                vectors.append(velocity)
                records.extend(
                    {
                        "state": str(state[i]),
                        "condition": condition,
                        "sample": f"{condition}_sample_{sample_number + 1}",
                        "corridor": arm,
                        "curve_parameter": float(t[i]),
                        "opposing_patch": bool(opposing_patch[i]),
                    }
                    for i in range(len(xy))
                )
            # Only four cells per selected condition: too few for the field's
            # minimum local count even when the two conditions are pooled.
            rare_xy = np.array([[-3.10, 2.55]]) + rng.normal(0, 0.05, (1, 2))
            points.append(rare_xy)
            vectors.append(np.array([[0.025, 0.045]]))
            records.append(
                {
                    "state": "Rare pocket",
                    "condition": condition,
                    "sample": f"{condition}_sample_{sample_number + 1}",
                    "corridor": "rare",
                    "curve_parameter": np.nan,
                    "opposing_patch": False,
                }
            )
    obs = pd.DataFrame(
        records, index=[f"field_cell_{i:05d}" for i in range(len(records))]
    )
    obs["state"] = pd.Categorical(obs.state, categories=STATES, ordered=True)
    obs["condition"] = pd.Categorical(
        obs.condition, categories=["control", "treated", "context"]
    )
    adata = AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
    adata.var_names = ["unused_placeholder"]
    adata.obsm["X_display"] = np.vstack(points)
    adata.obsm["velocity_display"] = np.vstack(vectors)
    graph = np.zeros((len(STATES), len(STATES)))
    edges = [
        ("Source", "Trunk", 0.90),
        ("Trunk", "Junction", 0.88),
        ("Junction", "Upper arm", 0.78),
        ("Junction", "Lower bend", 0.81),
        ("Lower bend", "Lower tip", 0.84),
    ]
    for left, right, weight in edges:
        i, j = STATES.index(left), STATES.index(right)
        graph[i, j] = graph[j, i] = weight
    adata.uns["paga_fixture"] = {
        "connectivities": graph,
        "groups": "state",
        "node_order": list(STATES),
        "source": "Constructed symmetric graph; upper island and rare pocket are intentionally disconnected.",
    }
    adata.uns["synthetic_field_demo"] = {
        "seed": int(seed),
        "scope": "Synthetic display-field API demonstration; no fitted embedding, velocity model, or biological findings.",
        "display_coordinates": "Constructed noisy curved corridors in two dimensions; not UMAP.",
        "vectors": "Analytic curve tangents, synthetic speed changes, noise, and a deliberate opposing-vector patch.",
        "topology": "Constructed symmetric fixture; not inferred PAGA connectivity.",
        "sample_design": "Four independent sample identities per condition; the descriptive field pools eligible cells.",
        "gap": {"x_open_interval": [1.10, 2.55], "y_above": 0.90},
        "opposing_patch": "Lower corridor, curve parameter in (0.38, 0.64); mixed forward/reverse vectors.",
        "rare_pocket": "Four cells per condition; below the local neighbor-support threshold.",
        "context": "Overlapping coordinates with large reversed vectors; excluded from the requested fields.",
        "condition_counts": {
            str(k): int(v) for k, v in obs.condition.value_counts().items()
        },
    }
    return adata


def run_gallery(output_dir, formats=("png", "svg", "pdf"), seed=71):
    """Save an overview, a support diagnostic, and a matched-condition figure."""
    if (
        not hasattr(sg.pl, "velocity_field")
        or "velocity_mode" not in inspect.signature(sg.pl.paga_shift_map).parameters
    ):
        raise RuntimeError(
            "This example needs the unreleased local-velocity-field API. "
            "Install the source checkout containing this example with python -m pip install -e ."
        )
    adata = make_velocity_demo(seed)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    xmin, xmax, ymin, ymax = GRID_EXTENT
    provenance = adata.uns["synthetic_field_demo"]

    def save(stem, result, caption):
        result.metadata["synthetic_fixture"] = provenance
        result.metadata["gallery_caption"] = caption
        # Preserve sample identities and injected-region flags beside the
        # plotted coordinate/vector sources, including the excluded context.
        result.data["fixture_observations"] = adata.obs.copy()
        result.fig.text(
            0.055,
            0.034,
            "SYNTHETIC DISPLAY + SUPPLIED VECTORS",
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
    common = dict(
        basis="display",
        velocity_basis="display",
        condition_key="condition",
        background="selected",
        max_points=2200,
        random_state=seed,
        xlim=(xmin, xmax),
        ylim=(ymin, ymax),
        aspect="equal",
        background_size=4.0,
        background_alpha=0.20,
        return_result=True,
        show=False,
    )
    with sg.pl.style_context("paper"):
        hero = sg.pl.paga_shift_map(
            adata,
            node_key="state",
            condition_key="condition",
            group0="control",
            group1="treated",
            basis="display",
            paga_key="paga_fixture",
            graph_node_order=STATES,
            min_cells=15,
            palette=STATE_PALETTE,
            background="contrast",
            max_points=2200,
            random_state=seed,
            xlim=(xmin, xmax),
            ylim=(ymin, ymax),
            aspect="equal",
            node_size=125,
            endpoint_size=36,
            edge_lw=1.5,
            edge_alpha=0.30,
            arrow_linewidth=1.4,
            arrow_head_size=8,
            background_size=4,
            background_alpha=0.18,
            show_velocity=True,
            velocity_basis="display",
            velocity_mode="stream",
            velocity_color="#14857F",
            velocity_options={
                **FIELD_OPTIONS,
                "stream_density": 1.35,
                "stream_linewidth": 1.15,
            },
            label_mode="auto",
            label_fontsize=8.5,
            legend_loc="outside",
            figsize=(11, 7.8),
            title="Local flow on a supplied state graph",
            return_result=True,
            show=False,
        )
        exports["overview"] = save(
            "velocity_paga_stream",
            hero,
            "Paths interpolate a supported display field; they are not cell trajectories or elapsed time. Graph edges are a constructed fixture.",
        )

        # Reserve the colorbar and caption gutters before drawing the field.
        support_fig = plt.figure(figsize=(10.7, 7.7))
        support_ax = support_fig.add_axes([0.07, 0.18, 0.78, 0.69])
        support_cax = support_fig.add_axes([0.90, 0.22, 0.022, 0.58])
        support = sg.pl.velocity_field(
            adata,
            **common,
            **FIELD_OPTIONS,
            groups=["control", "treated"],
            mode="grid",
            scale=4.0,
            color_by="coherence",
            cmap="YlGnBu",
            norm=Normalize(0, 1),
            arrow_linewidth=1.25,
            arrow_head_size=5.5,
            legend=False,
            ax=support_ax,
            cax=support_cax,
            title="Local support and direction agreement",
        )
        ax = support.axes["field"]
        grid = support.data["grid"]
        low_agreement = (
            grid.count_supported
            & ~grid.coherence_supported
            & np.isfinite(grid.coherence)
        )
        if low_agreement.any():
            support.artists["low_agreement_nodes"] = ax.scatter(
                grid.loc[low_agreement, "x"],
                grid.loc[low_agreement, "y"],
                s=28,
                facecolors="none",
                edgecolors="#C58A35",
                linewidths=0.8,
                zorder=5,
                label="Enough neighbors; low direction agreement",
            )
            ax.legend(loc="lower left", fontsize=7.5, frameon=False)
        note_style = dict(
            fontsize=8,
            color="#344C5C",
            ha="center",
            va="center",
            arrowprops={"arrowstyle": "-", "color": "#81919C", "linewidth": 0.8},
        )
        ax.annotate(
            "Empty coordinate gap\nfinite radius stops interpolation",
            xy=(1.82, 1.70),
            xytext=(1.4, 3.13),
            **note_style,
        )
        ax.annotate(
            "Mixed forward / reverse vectors\nlow local agreement",
            xy=(1.93, -2.12),
            xytext=(0.8, -3.07),
            **note_style,
        )
        ax.annotate(
            "Rare pocket\ninsufficient neighbors",
            xy=(-3.10, 2.55),
            xytext=(-4.45, 2.9),
            **note_style,
        )
        support.metadata["diagnostic_annotations"] = [
            "constructed empty gap",
            "injected opposing-vector patch",
            "rare pocket below minimum neighbors",
        ]
        support.metadata["low_agreement_marker_rule"] = (
            "count_supported & ~coherence_supported & finite(coherence)"
        )
        exports["support"] = save(
            "velocity_field_grid",
            support,
            "Arrows ×4. Color measures local direction agreement, not confidence. Only supported nodes passing the threshold have arrows.",
        )

        # Final geometry and a dedicated shared colorbar are set before plotting.
        fig = plt.figure(figsize=(12.2, 6.2))
        axes = [
            fig.add_axes([0.055, 0.20, 0.425, 0.66]),
            fig.add_axes([0.535, 0.20, 0.425, 0.66]),
        ]
        cax = fig.add_axes([0.38, 0.113, 0.24, 0.025])
        fig.suptitle(
            "Condition-specific fields on the same display", x=0.51, y=0.97, fontsize=15
        )
        speed_norm = Normalize(0, 0.08)
        facets = {}
        for i, condition in enumerate(("control", "treated")):
            facets[condition] = sg.pl.velocity_field(
                adata,
                **common,
                **FIELD_OPTIONS,
                groups=[condition],
                mode="grid",
                scale=4.0,
                color_by="speed",
                cmap="viridis_r",
                norm=speed_norm,
                arrow_linewidth=1.1,
                arrow_head_size=4.5,
                legend=False,
                ax=axes[i],
                cax=cax if i else None,
                show_colorbar=bool(i),
                title=condition.capitalize(),
            )
        comparison = sg.pl.PlotResult(
            fig=fig,
            axes={"control": axes[0], "treated": axes[1], "colorbar": cax},
            data={
                f"{condition}_{key}": table
                for condition, result in facets.items()
                for key, table in result.data.items()
            },
            metadata={
                "plot_type": "velocity_field_condition_comparison",
                "panels": {
                    condition: result.metadata for condition, result in facets.items()
                },
                "shared_grid_extent": list(GRID_EXTENT),
                "shared_display_scale": 4.0,
                "shared_speed_limits": [0, 0.08],
            },
            warnings=[
                warning for result in facets.values() for warning in result.warnings
            ],
            artists={condition: result.artists for condition, result in facets.items()},
        )
        exports["conditions"] = save(
            "velocity_condition_comparison",
            comparison,
            "Matched grid, finite radius, support thresholds, and speed scale; arrows ×4 in both panels. Descriptive pooled-cell fields.",
        )
    print(
        f"Wrote three synthetic velocity-field figures and full sources to {output_dir.resolve()}"
    )
    return exports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/velocity_fields")
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=["png", "svg", "pdf"],
        default=["png", "svg", "pdf"],
    )
    parser.add_argument("--seed", type=int, default=71)
    args = parser.parse_args()
    run_gallery(args.output_dir, formats=tuple(args.formats), seed=args.seed)


if __name__ == "__main__":
    main()
