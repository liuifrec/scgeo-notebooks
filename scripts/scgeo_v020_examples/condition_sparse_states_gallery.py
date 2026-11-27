"""Reproduce useful plots for a synthetic state with zero or one control cell.

Run from a checkout containing the condition-sparse plotting extensions::

    python examples/condition_sparse_states_gallery.py \
        --output-dir results/condition_sparse_states

All expression, sample identities, counts, coordinates and graph edges are
generated here. The generic sample names C1/C2/T1/T2 are synthetic. Expression
values are continuous arbitrary units, not counts. The symmetric connectivity
matrix is an explicitly supplied topology fixture, not inferred PAGA evidence.
No external datasets, downloads or private study inputs are used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from copy import deepcopy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.decomposition import PCA

import scgeo as sg

STATES = ["Shared", "Enriched", "Other"]
GENES = ["Ieg1", "Response2", "Identity1", "Identity2", "Housekeeping", "Background"]
SAMPLE_COUNTS = {
    "C1": {"Shared": 540, "Enriched": 0, "Other": 60},
    "C2": {"Shared": 29, "Enriched": 1, "Other": 30},
    "T1": {"Shared": 20, "Enriched": 60, "Other": 20},
    "T2": {"Shared": 90, "Enriched": 120, "Other": 90},
}
CONDITION_COLORS = {"control": "#376D91", "treated": "#D97943"}
SAMPLE_COLORS = {"C1": "#376D91", "C2": "#8AB6CF", "T1": "#D97943", "T2": "#EAB777"}
STATE_COLORS = {"Shared": "#647E98", "Enriched": "#B9835F", "Other": "#7B9B8C"}


def make_example_data(seed=83):
    """Build three states with unequal capture and a declared affine PCA map."""
    rng = np.random.default_rng(seed)
    centers = {
        "Shared": np.array([-2.0, 0.0]),
        "Enriched": np.array([1.5, 1.65]),
        "Other": np.array([1.6, -1.85]),
    }
    matrices, displays, records = [], [], []
    for sample, state_counts in SAMPLE_COUNTS.items():
        condition = "control" if sample.startswith("C") else "treated"
        treated = float(condition == "treated")
        sample_offset = rng.normal(0, 0.10, len(GENES))
        for state, n_cells in state_counts.items():
            if not n_cells:
                continue
            display = centers[state] + rng.normal(0, [0.26, 0.22], (n_cells, 2))
            baseline = rng.normal(0, 0.45, n_cells)
            values = np.empty((n_cells, len(GENES)))
            values[:, 0] = (
                3.0 + 1.5 * (state == "Enriched") + 0.70 * treated + 0.25 * baseline
            )
            values[:, 1] = 3.0 + 0.8 * (state == "Enriched") + 0.55 * treated
            values[:, 2] = 4.0 + 1.1 * display[:, 0]
            values[:, 3] = 4.0 + 0.8 * display[:, 1]
            values[:, 4] = 3.0
            values[:, 5] = 2.0 + 0.15 * baseline
            values += sample_offset + rng.normal(0, 0.12, values.shape)
            matrices.append(values)
            displays.append(display)
            records.extend(
                {
                    "sample": sample,
                    "condition": condition,
                    "state": state,
                    "baseline_score": float(value),
                }
                for value in baseline
            )
    obs = pd.DataFrame(
        records, index=[f"synthetic_cell_{i:04d}" for i in range(len(records))]
    )
    obs["state"] = pd.Categorical(obs.state, categories=STATES, ordered=True)
    obs["condition"] = pd.Categorical(obs.condition, categories=["control", "treated"])
    obs["sample"] = pd.Categorical(obs["sample"], categories=list(SAMPLE_COUNTS))
    adata = AnnData(np.vstack(matrices), obs=obs, var=pd.DataFrame(index=GENES))
    pca = PCA(n_components=4, svd_solver="full")
    adata.obsm["X_pca"] = pca.fit_transform(adata.X)
    adata.varm["PCs"] = pca.components_.T.copy()
    adata.obsm["X_display"] = np.vstack(displays)
    adata.uns["condition_colors"] = list(CONDITION_COLORS.values())
    adata.uns["sample_colors"] = list(SAMPLE_COLORS.values())
    adata.uns["state_colors"] = list(STATE_COLORS.values())
    adata.uns["paga_fixture"] = {
        "groups": "state",
        "node_order": list(STATES),
        "connectivities": np.array(
            [[0.0, 0.80, 0.45], [0.80, 0.0, 0.20], [0.45, 0.20, 0.0]]
        ),
        "provenance": "Constructed symmetric topology fixture; no PAGA inference or transition inference was run.",
    }
    adata.uns["synthetic_provenance"] = {
        "kind": "fully synthetic",
        "seed": int(seed),
        "expression": "continuous arbitrary units; generic feature names; not raw counts",
        "representation": "four-component PCA fitted once; exact affine loadings supplied",
        "display": "known synthetic two-dimensional state coordinates; not UMAP",
        "topology": "supplied symmetric graph fixture; not inferred PAGA evidence",
        "sample_design": "two independent synthetic samples per condition with unequal capture",
        "biological_validation": False,
    }
    return adata


def sample_state_counts(adata):
    """Keep zero-state samples and declare the full-sample denominator."""
    rows = []
    for sample in SAMPLE_COUNTS:
        condition = "control" if sample.startswith("C") else "treated"
        selected = adata.obs.loc[adata.obs["sample"].eq(sample)]
        denominator = len(selected)
        for state in STATES:
            count = int(selected.state.eq(state).sum())
            rows.append(
                {
                    "sample": sample,
                    "condition": condition,
                    "state": state,
                    "count": count,
                    "denominator": denominator,
                    "fraction": count / denominator if denominator else np.nan,
                    "denominator_scope": "all synthetic Shared/Enriched/Other cells in the sample",
                }
            )
    return pd.DataFrame(rows)


def _explain(adata, seed):
    sg.tl.robust_shift(
        adata,
        rep="X_pca",
        condition_key="condition",
        group0="control",
        group1="treated",
        by="state",
        sample_key="sample",
        center="geometric_median",
        n_boot=0,
        seed=seed,
    )
    return sg.tl.explain_shift(
        adata,
        modules={"Response": ["Ieg1", "Response2"]},
        min_cells=10,
        min_samples=2,
        n_boot=0,
        seed=seed,
    )


def _paga(adata, *, pie_mode="cell_fraction", pie_key="condition", title=None):
    return sg.pl.paga_shift_map(
        adata,
        node_key="state",
        condition_key="condition",
        group0="control",
        group1="treated",
        basis="display",
        paga_key="paga_fixture",
        node_position="pooled",
        show_shift=False,
        show_endpoints=False,
        show_velocity=False,
        node_color_mode="constant",
        constant_node_color="#566777",
        palette=STATE_COLORS,
        pie_key=pie_key,
        pie_mode=pie_mode,
        pie_sample_key="sample" if pie_mode == "sample_abundance" else None,
        pie_palette=SAMPLE_COLORS if pie_key == "sample" else CONDITION_COLORS,
        node_scale_by_n=False,
        node_size=780,
        pie_size_scale=1.25,
        background_size=7,
        background_alpha=0.10,
        min_cells=10,
        connectivity_threshold=0.15,
        figsize=(8.6, 5.2),
        title=title,
        label_mode="on_data",
        label_offsets={"Shared": (0, 30), "Enriched": (0, 30), "Other": (0, -34)},
        label_fontsize=10,
        return_result=True,
        show=False,
    )


def _duplicate_capture(adata, sample="T2", copies=5):
    """Uniformly multiply one sample's capture without changing its composition."""
    selected = np.flatnonzero(adata.obs["sample"].eq(sample).to_numpy())
    indices = np.r_[np.arange(adata.n_obs), np.tile(selected, copies - 1)]
    obs = adata.obs.iloc[indices].copy()
    obs.index = [f"duplicated_cell_{i:04d}" for i in range(len(indices))]
    repeated = AnnData(adata.X[indices].copy(), obs=obs, var=adata.var.copy())
    for key, coordinates in adata.obsm.items():
        repeated.obsm[key] = coordinates[indices].copy()
    for key, loadings in adata.varm.items():
        repeated.varm[key] = loadings.copy()
    for key in [
        "paga_fixture",
        "synthetic_provenance",
        "condition_colors",
        "sample_colors",
        "state_colors",
    ]:
        repeated.uns[key] = deepcopy(adata.uns[key])
    return repeated


def _pie_matrix(result):
    return (
        result.data["composition"].set_index("node").loc[STATES, ["control", "treated"]]
    )


def _normalized_or_missing(result):
    table = result.data["composition"].set_index("node")
    values = table.to_numpy(dtype=float)
    return all(
        np.isnan(row).all() or (np.isfinite(row).all() and np.isclose(row.sum(), 1.0))
        for row in values
    )


def run_gallery(output_dir, formats=("png", "svg", "pdf"), *, seed=83, dpi=110):
    """Render four figures, export sources, and audit the declared sparse case."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = make_example_data(seed)
    explanation = _explain(adata, seed)
    count_table = sample_state_counts(adata)
    count_table.to_csv(output_dir / "synthetic_sample_state_counts.csv", index=False)
    adata.obs.to_csv(output_dir / "synthetic_cells.csv")
    pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names).to_csv(
        output_dir / "synthetic_expression.csv"
    )
    pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names).to_csv(
        output_dir / "synthetic_pca.csv"
    )
    pd.DataFrame(adata.varm["PCs"], index=adata.var_names).to_csv(
        output_dir / "synthetic_pca_loadings.csv"
    )
    for name, table in explanation.items():
        if isinstance(table, pd.DataFrame):
            table.to_csv(output_dir / f"explanation_{name}.csv", index=False)

    expression = sg.pl.shift_curves(
        adata,
        var_names=["Ieg1"],
        state="Enriched",
        basis="display",
        panels=("distribution", "expression"),
        palette=CONDITION_COLORS,
        title="Rare-control state: observed expression remains visible",
        figsize=(8.8, 4.3),
        size=14,
        alpha=0.7,
        return_result=True,
        show=False,
    )
    expression.fig.get_layout_engine().set(rect=(0, 0.06, 1, 0.94))
    raw = _paga(adata, title="Synthetic state graph · captured-cell fractions")
    balanced = _paga(
        adata,
        pie_mode="sample_abundance",
        title="Synthetic state graph · equal sample weighting",
    )
    sample_identity = _paga(
        adata, pie_key="sample", title="Synthetic state graph · sample identity"
    )

    duplicated = _duplicate_capture(adata)
    duplicated_raw = _paga(duplicated)
    duplicated_balanced = _paga(duplicated, pie_mode="sample_abundance")
    raw_difference = float(
        np.max(np.abs((_pie_matrix(duplicated_raw) - _pie_matrix(raw)).to_numpy()))
    )
    balanced_difference = float(
        np.max(
            np.abs(
                (_pie_matrix(duplicated_balanced) - _pie_matrix(balanced)).to_numpy()
            )
        )
    )
    audit = pd.concat(
        {
            "original_cell_fraction": _pie_matrix(raw),
            "duplicated_cell_fraction": _pie_matrix(duplicated_raw),
            "original_sample_abundance": _pie_matrix(balanced),
            "duplicated_sample_abundance": _pie_matrix(duplicated_balanced),
        },
        names=["case", "state"],
    ).reset_index()
    audit.to_csv(output_dir / "capture_duplication_pie_audit.csv", index=False)
    sample_state_counts(duplicated).to_csv(
        output_dir / "duplicated_sample_state_counts.csv", index=False
    )
    plt.close(duplicated_raw.fig)
    plt.close(duplicated_balanced.fig)

    status = explanation["summary"].set_index("state").loc["Enriched"]
    state_result = next(
        value
        for value in explanation["state_results"].values()
        if value["state"] == "Enriched"
    )
    enriched_genes = explanation["genes"].loc[explanation["genes"].state.eq("Enriched")]
    distributions = expression.data["sample_expression"].set_index("sample")
    coverage = expression.data["coverage"].set_index("sample")
    balanced_counts = balanced.data["pie_sample_counts"]
    zero_count = balanced_counts.loc[
        balanced_counts.state.eq("Enriched") & balanced_counts["sample"].eq("C1")
    ].iloc[0]
    checks = {
        "insufficient_cells_status_preserved": status["status"] == "insufficient_cells",
        "zero_and_one_control_cells_preserved": int(coverage.loc["C1", "n_cells"]) == 0
        and int(coverage.loc["C2", "n_cells"]) == 1,
        "no_invented_direction": bool(np.isnan(state_result["direction"]).all()),
        "no_invented_contributions": bool(enriched_genes.contribution.isna().all()),
        "absent_control_expression_is_missing": bool(
            distributions.loc["C1", ["mean", "median", "q25", "q75"]].isna().all()
        ),
        "singleton_is_a_point_without_interval": bool(
            np.isfinite(distributions.loc["C2", "median"])
            and distributions.loc["C2", ["q25", "q75"]].isna().all()
        ),
        "no_invented_shift_curve": bool(
            expression.data["curves"].expression.isna().all()
        ),
        "sample_balanced_zero_count_keeps_denominator": int(zero_count["count"]) == 0
        and int(zero_count["denominator"]) == 600
        and zero_count["state_fraction"] == 0,
        "sample_abundance_invariant_to_capture_duplication": balanced_difference
        < 1e-12,
        "raw_cell_pies_respond_to_capture_duplication": raw_difference > 0.05,
        "pie_slices_normalized_or_missing": all(
            _normalized_or_missing(result)
            for result in [raw, balanced, sample_identity]
        ),
        "enriched_pie_visible_in_all_modes": all(
            "Enriched" in result.artists["pies"]
            for result in [raw, balanced, sample_identity]
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"Synthetic sparse-state gallery checks failed: {checks}")

    results = {
        "sparse_state_expression": expression,
        "paga_condition_raw_cells": raw,
        "paga_condition_sample_abundance": balanced,
        "paga_sample_identity": sample_identity,
    }
    exports = {}
    for name, result in results.items():
        result.metadata["data_provenance"] = dict(adata.uns["synthetic_provenance"])
        result.fig.text(
            0.012,
            0.014,
            "Fully synthetic · supplied topology fixture · no inferred transitions or biological validation",
            fontsize=7,
            color="#5F6C78",
        )
        files = result.save(output_dir / name, formats=formats, dpi=dpi)
        exports[name] = str(Path(files["metadata"]).name)
        plt.close(result.fig)

    metrics = {
        "n_cells": adata.n_obs,
        "enriched_counts": {
            str(sample): int(value) for sample, value in coverage.n_cells.items()
        },
        "raw_pie_maximum_change_after_fivefold_T2_capture": raw_difference,
        "sample_abundance_maximum_change_after_fivefold_T2_capture": balanced_difference,
        "original_condition_cell_fractions": json.loads(
            _pie_matrix(raw).reset_index().to_json(orient="records")
        ),
        "original_condition_sample_abundance_slices": json.loads(
            _pie_matrix(balanced).reset_index().to_json(orient="records")
        ),
    }
    summary = {
        "schema_version": 1,
        "data_kind": "fully synthetic; generic samples; no private inputs or downloads",
        "seed": int(seed),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": {name: bool(value) for name, value in checks.items()},
        "metrics": metrics,
        "mapping": {
            **explanation["mapping"],
            "intercept": explanation["mapping"]["intercept"].tolist(),
        },
        "source_parameters": explanation["source_params"],
        "plot_manifests": exports,
        "interpretation": {
            "expression": "Observed cells only; absent control cells have missing expression, not zero expression.",
            "cell_fraction": "Fraction of captured cells in each state assigned to each displayed category.",
            "sample_abundance": "Mean per-sample state fractions within each condition, then normalized across conditions for pie slices.",
            "denominator": "All annotated synthetic Shared/Enriched/Other cells in each sample.",
            "graph": "Supplied symmetric fixture on pooled display anchors; no condition shift or velocity arrows.",
            "validation": "Software semantics on synthetic input; not evidence for any biological mechanism.",
        },
        "reproduce": f"python examples/condition_sparse_states_gallery.py --seed {seed} --output-dir NEW_OUTPUT_DIRECTORY",
    }
    (output_dir / "gallery_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/condition_sparse_states")
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "svg", "pdf"),
        default=("png", "svg", "pdf"),
    )
    parser.add_argument("--seed", type=int, default=83)
    parser.add_argument("--dpi", type=int, default=110)
    args = parser.parse_args()
    result = run_gallery(args.output_dir, args.formats, seed=args.seed, dpi=args.dpi)
    print(
        json.dumps(
            {
                "checks": result["checks"],
                "metrics": result["metrics"],
                "plot_manifests": result["plot_manifests"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
