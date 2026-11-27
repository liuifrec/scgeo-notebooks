"""Generate stored-evidence compass and readability examples on synthetic data.

Run from an installed development checkout. No external data are downloaded.
The sparse-state input, display coordinates, graph and velocities are synthetic;
actual ScGeo analysis APIs generate the stored numerical evidence. Optional
timings are observational measurements, never CI performance thresholds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from copy import deepcopy
from io import BytesIO
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from condition_sparse_states_gallery import STATES, make_example_data
from matplotlib.colors import Normalize

import scgeo as sg


def prepare_example(seed=83):
    """Use real analysis outputs, with explicit constructed-vector provenance."""
    adata = make_example_data(seed)
    rng = np.random.default_rng(seed + 1)
    rotation = np.array(
        [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]],
        dtype=float,
    )
    adata.obsm["X_rotated"] = adata.obsm["X_pca"] @ rotation
    velocity = np.tile([0.05, 0.01, -0.01, 0.02], (adata.n_obs, 1))
    velocity[adata.obs.state.eq("Other")] *= -1
    adata.obsm["V_pca"] = velocity
    adata.obsm["V_rotated"] = velocity @ rotation
    # Independently declared synthetic display vectors, not fitted RNA kinetics.
    adata.obsm["velocity_display"] = np.tile([0.03, 0.015], (adata.n_obs, 1))
    adata.obs["ood_score"] = rng.uniform(0.0, 1.0, adata.n_obs)
    adata.obs.loc[adata.obs.state.eq("Other"), "ood_score"] = np.nan
    adata.obs["marker_score"] = np.asarray(adata.X[:, 0]).reshape(-1)

    sg.tl.robust_shift(
        adata,
        rep="X_pca",
        condition_key="condition",
        group0="control",
        group1="treated",
        by="state",
        sample_key="sample",
        center="geometric_median",
        n_boot=40,
        seed=seed,
    )
    sg.tl.representation_stability(
        adata,
        reps=["X_pca", "X_rotated"],
        node_key="state",
        condition_key="condition",
        group0="control",
        group1="treated",
        sample_key="sample",
        center="geometric_median",
        velocity_keys={"X_pca": "V_pca", "X_rotated": "V_rotated"},
        n_boot=16,
        min_cells=10,
        seed=seed,
    )
    report = sg.get.state_report(adata, node_key="state")
    report["state_name"] = report.state.map(
        {
            "Shared": "Shared population",
            "Enriched": "Enriched population",
            "Other": "Other population",
        }
    )
    source = adata.uns["scgeo"]["representation_stability"]
    alignment = (
        source["per_rep_state"].loc[source["per_rep_state"].rep.eq("X_pca")].copy()
    )
    alignment.attrs["params"] = deepcopy(source["params"])
    # These describe the actual implementation used above, not a new estimator.
    alignment.attrs["provenance"] = {
        "velocity_scope": "contrast",
        "velocity_aggregation": "cell_mean",
        "velocity_weighting": "cells",
    }
    return adata, report, alignment


def compass(adata, report=None, alignment=None, **kwargs):
    options = {
        "node_key": "state",
        "condition_key": "condition",
        "group0": "control",
        "group1": "treated",
        "basis": "display",
        "paga_key": "paga_fixture",
        "graph_node_order": STATES,
        "node_position": "pooled",
        "sample_key": "sample",
        "min_cells": 10,
        "node_size_mode": "constant",
        "node_size_scale": 520,
        "arrow_color_mode": "constant",
        "arrow_color": "#455E70",
        "label_fontsize": 9,
        "ood_key": "ood_score",
        "ood_threshold": 0.8,
        "return_result": True,
        "show": False,
    }
    if report is not None:
        options.update(
            mode="stored",
            shift_df=report,
            alignment_df=alignment,
            analysis_rep="X_pca",
            effect_metric="delta_norm",
        )
    options.update(kwargs)
    return sg.pl.recovery_compass(adata, **options)


def save_result(result, out, name, formats, *, dpi=110):
    result.metadata["synthetic_example"] = {
        "data": "generic simulated samples and continuous expression",
        "topology": "constructed symmetric fixture",
        "velocity": "constructed vectors, not fitted RNA kinetics",
        "representations": "PCA and a related orthogonal rotation; no independent biological replication",
        "validation": "software semantics only; no biological mechanism inferred",
    }
    result.fig.text(
        0.01,
        0.005,
        "Synthetic data · supplied topology and vectors · stored analysis evidence",
        fontsize=6.5,
        color="#607080",
    )
    exported = result.save(out / name, formats=formats, dpi=dpi)
    plt.close(result.fig)
    return exported


def add_unrelated_condition(adata):
    """Append opposite vectors without changing either selected condition."""
    positions = np.flatnonzero(adata.obs.state.eq("Shared"))[:30]
    extra = adata.obs.iloc[positions].copy()
    extra["condition"] = "unrelated"
    extra["sample"] = "U1"
    obs = pd.concat([adata.obs, extra], ignore_index=True)
    obs.index = [f"scope_cell_{i}" for i in range(len(obs))]
    out = AnnData(
        np.vstack([adata.X, adata.X[positions]]), obs=obs, var=adata.var.copy()
    )
    for key, values in adata.obsm.items():
        appended = values[positions].copy()
        if key == "velocity_display":
            appended[:] = [-100, -100]
        out.obsm[key] = np.vstack([values, appended])
    out.uns["paga_fixture"] = deepcopy(adata.uns["paga_fixture"])
    return out


def rendering_benchmark(adata, repeats=3):
    """Record full/budgeted drawing with identical full-data selection values."""
    indices = np.tile(np.arange(adata.n_obs), 20)
    obs = adata.obs.iloc[indices].copy()
    obs.index = [f"benchmark_cell_{i}" for i in range(len(indices))]
    large = AnnData(np.asarray(adata.X)[indices], obs=obs, var=adata.var.copy())
    large.obsm["X_display"] = adata.obsm["X_display"][indices]
    large.uns["state_colors"] = list(adata.uns["state_colors"])
    records, summaries = [], []
    for budget in (None, 1200):
        draws, exports, sizes = [], [], []
        displayed = None
        for _ in range(repeats):
            start = perf_counter()
            result = sg.pl.consensus_structure(
                large,
                score_key="marker_score",
                basis="display",
                groupby="state",
                topk=100,
                sample_key="sample",
                max_points=budget,
                rasterized=True,
                random_state=83,
                colorbar=False,
                print_legend=False,
                return_result=True,
                show=False,
            )
            result.fig.canvas.draw()
            draws.append(perf_counter() - start)
            stream = BytesIO()
            start = perf_counter()
            result.fig.savefig(stream, format="svg")
            exports.append(perf_counter() - start)
            sizes.append(stream.tell())
            if displayed is None:
                displayed = int(result.data["cells"].background_displayed.sum())
                summaries.append(result.data["summary"].copy())
            plt.close(result.fig)
        records.append(
            {
                "max_points": budget,
                "rasterized": True,
                "input_cells": large.n_obs,
                "displayed_background_cells": displayed,
                "median_create_and_draw_seconds": float(np.median(draws)),
                "median_svg_export_seconds": float(np.median(exports)),
                "median_svg_bytes": int(np.median(sizes)),
                "repeats": repeats,
            }
        )
    pd.testing.assert_frame_equal(summaries[0], summaries[1])
    return {
        "summary_invariant": True,
        "measurements": records,
        "interpretation": "Observed timings in this environment, with no pass/fail speed threshold.",
    }


def run_gallery(output_dir, formats=("png", "svg", "pdf"), *, seed=83, benchmark=False):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata, report, alignment = prepare_example(seed)
    before_report, before_alignment = report.copy(deep=True), alignment.copy(deep=True)
    manifests, checks = {}, {}
    with sg.pl.style_context("paper"):
        result = compass(
            adata,
            report,
            alignment,
            layout="summary",
            max_points=450,
            ring_mode="none",
            figsize=(15.5, 7.2),
            title="State response compass · stored evidence and individual samples",
        )
        evidence = result.data["evidence"].copy()
        enriched = evidence.set_index("state").loc["Enriched"]
        checks["sparse_source_effect_not_presented_as_supported"] = bool(
            not enriched.effect_available and np.isfinite(enriched.effect)
        )
        checks["supported_effects_remain_available"] = bool(
            evidence.set_index("state")
            .loc[["Shared", "Other"], "effect_available"]
            .all()
        )
        checks["pooled_sparse_node_remains_visible"] = (
            "Enriched" in result.artists["nodes"]
        )
        sample_counts = result.data["sample_abundance"]
        c1 = sample_counts.loc[
            sample_counts.state.eq("Enriched") & sample_counts["sample"].eq("C1")
        ].iloc[0]
        checks["zero_state_sample_retains_denominator"] = bool(
            c1["count"] == 0 and c1.denominator == 600 and c1.state_fraction == 0
        )
        manifests["stored_compass_summary"] = save_result(
            result, out, "stored_compass_summary", formats
        )

        changed = adata.copy()
        changed.obsm["X_display"] = changed.obsm["X_display"] @ np.array(
            [[0, 2], [-0.3, 0]]
        ) + [4, -2]
        changed.obsm["velocity_display"] = -changed.obsm["velocity_display"]
        alternate = compass(
            changed,
            report,
            alignment,
            layout="summary",
            max_points=50,
            ring_mode="none",
            legend=False,
        )
        pd.testing.assert_frame_equal(evidence, alternate.data["evidence"])
        checks["stored_evidence_invariant_to_display_change"] = True
        plt.close(alternate.fig)

        ood = compass(
            adata,
            report,
            alignment,
            ring_mode="ood_frac",
            title="Optional OOD overlay · explicit threshold and missing score coverage",
        )
        checks["missing_ood_is_unavailable"] = not bool(
            ood.data["glyphs"].set_index("node").loc["Other", "ring_available"]
        )
        manifests["compass_ood_coverage"] = save_result(
            ood, out, "compass_ood_coverage", formats
        )

        plain = compass(adata, ring_mode="none", legend=False)
        unrelated = compass(
            add_unrelated_condition(adata), ring_mode="none", legend=False
        )
        left = plain.data["nodes"].set_index("state").alignment_cosine
        right = unrelated.data["nodes"].set_index("state").alignment_cosine
        np.testing.assert_allclose(left.loc[STATES], right.loc[STATES], equal_nan=True)
        checks["unrelated_condition_does_not_change_default_alignment"] = True
        plt.close(plain.fig)
        plt.close(unrelated.fig)

        concordance = sg.pl.geometry_velocity_concordance(
            report,
            alignment,
            alignment_state_key="node",
            effect_metric="delta_norm",
            analysis_rep="X_pca",
            label_key="state_name",
            min_cells=10,
            min_samples=2,
            sample_count_keys=("n_samples0", "n_samples1"),
            highlight_states=["Shared"],
            max_labels=1,
            return_result=True,
            show=False,
            title="Stored state evidence · supported values and unavailable states",
        )
        manifests["named_concordance"] = save_result(
            concordance, out, "named_concordance", formats
        )

        selection_options = {
            "score_key": "marker_score",
            "basis": "display",
            "groupby": "state",
            "selection": "threshold",
            "threshold": 4.5,
            "sample_key": "sample",
            "random_state": seed,
            "hi_alpha": 0.9,
            "outline_alpha": 0.35,
            "norm": Normalize(2.5, 5.8),
            "print_legend": False,
            "return_result": True,
            "show": False,
        }
        selected = sg.pl.consensus_structure(adata, max_points=300, **selection_options)
        full = sg.pl.consensus_structure(adata, **selection_options)
        pd.testing.assert_frame_equal(selected.data["summary"], full.data["summary"])
        checks["selection_prevalence_invariant_to_display_budget"] = True
        checks["highlight_face_and_outline_opacity_independent"] = bool(
            selected.artists["selected"].get_alpha() == 0.9
            and selected.artists["selected_outlines"].get_alpha() == 0.35
        )
        plt.close(full.fig)
        manifests["threshold_selection"] = save_result(
            selected, out, "threshold_selection", formats
        )

        atlas = sg.pl.geometry_dynamics_atlas(
            adata,
            basis="display",
            color="marker_score",
            node_key="state",
            condition_key="condition",
            group0="control",
            group1="treated",
            sample_key="sample",
            report=report,
            effect_metric="delta_norm",
            panels=("states", "velocity", "evidence"),
            velocity_mode="grid",
            velocity_options={
                "grid_size": 25,
                "n_neighbors": 12,
                "radius": 0.55,
                "min_neighbors": 3,
            },
            cmap="viridis",
            norm=Normalize(2.5, 5.8),
            max_points=450,
            highlight_states=["Enriched"],
            max_labels=2,
            layout="adaptive",
            title="Shared numeric scale · highlighted sparse state and source coverage",
            return_result=True,
            show=False,
        )
        checks["highlighted_atlas_state_is_labelled"] = bool(
            atlas.data["labels"].set_index("state").loc["Enriched", "displayed"]
        )
        checks["adaptive_atlas_resolves_evidence_strip"] = (
            atlas.metadata["resolved_layout"] == "adaptive_evidence_strip"
        )
        manifests["compact_atlas"] = save_result(atlas, out, "compact_atlas", formats)

    pd.testing.assert_frame_equal(report, before_report)
    pd.testing.assert_frame_equal(alignment, before_alignment)
    checks["stored_source_tables_unchanged"] = True
    if not all(checks.values()):
        raise AssertionError(
            f"Failed gallery checks: {[k for k, v in checks.items() if not v]}"
        )
    report.to_csv(out / "stored_state_report.csv", index=False)
    alignment.to_csv(out / "stored_alignment.csv", index=False)
    adata.obs.to_csv(out / "synthetic_observations.csv")
    record = {
        "schema_version": 1,
        "data_kind": "fully synthetic; no private input or external download",
        "seed": seed,
        "n_cells": adata.n_obs,
        "checks": checks,
        "source_report_provenance": report.attrs["provenance"],
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "implementation_sha256": {
            name: hashlib.sha256(
                (Path(sg.__file__).parent / "pl" / name).read_bytes()
            ).hexdigest()
            for name in (
                "_recovery_compass.py",
                "_compass_data.py",
                "_stored_evidence.py",
                "_geometry_velocity_concordance.py",
                "_consensus_structure.py",
                "_geometry_dynamics_atlas.py",
                "_paga_layout.py",
            )
        },
        "input_generator_sha256": hashlib.sha256(
            Path(__file__).with_name("condition_sparse_states_gallery.py").read_bytes()
        ).hexdigest(),
        "plot_manifests": {
            name: Path(value["metadata"]).name for name, value in manifests.items()
        },
        "benchmark": rendering_benchmark(adata) if benchmark else None,
        "reproduce": "python examples/compass_evidence_gallery.py --output-dir NEW_OUTPUT_DIRECTORY --formats png svg pdf",
        "interpretation": "Stored effects and eligible intervals are read from real ScGeo outputs on synthetic input; no inference about biological recovery.",
    }
    (out / "validation.json").write_text(
        json.dumps(record, indent=2, default=str) + "\n"
    )
    print(
        json.dumps(
            {
                "checks_passed": len(checks),
                "figures": len(manifests),
                "output_dir": str(out),
            },
            indent=2,
        )
    )
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/compass_upgrade")
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "svg", "pdf"],
        choices=["png", "svg", "pdf"],
    )
    parser.add_argument("--seed", type=int, default=83)
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    run_gallery(
        args.output_dir, tuple(args.formats), seed=args.seed, benchmark=args.benchmark
    )


if __name__ == "__main__":
    main()
