"""Reproducible synthetic examples for ScGeo shift explanations.

The inputs are simulated continuous expression values in arbitrary units,
not raw counts, observed biology, or an RNA-velocity trajectory. PCA is fitted
once to these inputs, and its affine transformation is supplied explicitly.
No downloads or private study data are needed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.decomposition import PCA

import scgeo as sg


def make_example_data(seed=37):
    """Simulate samples with global and opposing local expression responses.

    Unequal cell counts exercise the sample-level estimator. ``position`` is
    a known simulation coordinate, not an inferred lineage or pseudotime.
    """
    rng = np.random.RandomState(seed)
    programs = {
        "Response": [f"Response_{i}" for i in range(1, 4)],
        "LocalResponse": [f"LocalResponse_{i}" for i in range(1, 4)],
        "Identity": [f"Identity_{i}" for i in range(1, 4)],
        "SampleVariation": [f"SampleVariation_{i}" for i in range(1, 4)],
        "Null": [f"Null_{i}" for i in range(1, 7)],
    }
    genes = [gene for names in programs.values() for gene in names]
    matrices, records, truth_records = [], [], []
    for condition_index, condition in enumerate(("control", "treated")):
        for sample_index in range(6):
            sample = f"{condition}_{sample_index + 1}"
            n = 90 + 9 * sample_index
            position = np.linspace(-1.0, 1.0, n)
            position += rng.normal(0, 0.012, n)
            sample_offset = rng.normal(0, 0.09, len(genes))
            values = 4 + rng.normal(0, 0.11, (n, len(genes)))
            values += sample_offset
            values[:, :3] += 0.8 * position[:, None]
            values[:, :3] += 0.8 * condition_index
            values[:, 3:6] += 0.25 * position[:, None]
            values[:, 3:6] += 1.15 * condition_index * position[:, None]
            values[:, 6:9] += 1.2 * position[:, None]
            values[:, 9:12] += rng.normal(0, 0.28, 3)
            matrices.append(values)
            for cell_index, t in enumerate(position):
                cell = f"{sample}_cell_{cell_index:03d}"
                records.append({
                    "cell": cell,
                    "state": "Simulated state",
                    "condition": condition,
                    "sample": sample,
                    "position": float(t),
                    "region": "early" if t < 0 else "late",
                })
                truth_records.append({
                    "cell": cell,
                    "position": float(t),
                    "response_effect": 0.8,
                    "local_response_effect": float(1.15 * t),
                    "identity_effect": 0.0,
                    "null_effect": 0.0,
                })
    obs = pd.DataFrame.from_records(records).set_index("cell")
    obs["condition"] = pd.Categorical(obs.condition, ["control", "treated"])
    obs["state"] = pd.Categorical(obs.state)
    obs["region"] = pd.Categorical(obs.region, ["early", "late"])
    adata = AnnData(
        X=np.vstack(matrices),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    pca = PCA(n_components=8, svd_solver="full")
    adata.obsm["X_pca"] = pca.fit_transform(adata.X)
    adata.obsm["X_display"] = adata.obsm["X_pca"][:, :2].copy()
    adata.varm["PCs"] = pca.components_.T.copy()
    adata.uns["pca"] = {"params": {"zero_center": True}}
    adata.uns["condition_colors"] = ["#4477AA", "#D67732"]
    adata.uns["synthetic_provenance"] = {
        "kind": "fully synthetic",
        "seed": int(seed),
        "expression": "continuous positive values in arbitrary units; not counts",
        "representation": "eight-component sklearn PCA fitted once",
        "display": "first two fitted PCs; not UMAP",
        "samples": "six simulated independent samples per condition",
        "local_effect": "1.15 times known simulation position; opposing signs",
        "biological_validation": False,
    }
    feature_truth = pd.DataFrame([
        {
            "gene": gene,
            "program": program,
            "planted_global_effect": 0.8 if program == "Response" else 0.0,
            "planted_local_effect": "1.15 * position" if program == "LocalResponse" else "0",
            "interpretation": (
                "condition response" if program in {"Response", "LocalResponse"}
                else "sample-level nuisance" if program == "SampleVariation"
                else "baseline identity" if program == "Identity"
                else "unperturbed noise"
            ),
        }
        for program, names in programs.items()
        for gene in names
    ])
    return adata, programs, pd.DataFrame(truth_records), feature_truth


def distributional_counterexample():
    """Identical gene marginals can coexist with altered joint geometry."""
    control = np.array([[0.0, 0.0], [2.0, 2.0]])
    treated = np.array([[0.0, 2.0], [2.0, 0.0]])

    def distances(left, right):
        return np.linalg.norm(left[:, None, :] - right[None, :, :], axis=2)

    energy = (
        2 * distances(control, treated).mean()
        - distances(control, control).mean()
        - distances(treated, treated).mean()
    )
    return {
        "control": control.tolist(),
        "treated": treated.tolist(),
        "gene_marginals_identical": bool(
            np.array_equal(np.sort(control, axis=0), np.sort(treated, axis=0))
        ),
        "mean_displacement": float(np.linalg.norm(treated.mean(0) - control.mean(0))),
        "wasserstein_1_euclidean": 2.0,
        "energy_statistic": float(energy),
        "interpretation": "distributional change without a centroid displacement axis",
    }


def _input_hash(adata):
    digest = hashlib.sha256(adata.obs.to_csv().encode("utf-8"))
    digest.update(np.asarray(adata.X, dtype="<f8").tobytes())
    digest.update("\n".join(adata.var_names).encode("utf-8"))
    return digest.hexdigest()


def _explain(adata, programs, *, n_boot=80, seed=37,
             rep="X_pca", loadings_key="PCs", center="geometric_median"):
    sg.tl.robust_shift(
        adata,
        rep=rep,
        condition_key="condition",
        group0="control",
        group1="treated",
        by="state",
        sample_key="sample",
        center=center,
        n_boot=0,
        seed=seed,
    )
    return sg.tl.explain_shift(
        adata,
        shift_key="robust_shift",
        loadings_key=loadings_key,
        modules=programs,
        local_key="region",
        min_cells=10,
        min_samples=2,
        n_boot=n_boot,
        seed=seed,
    )


def _gene_contributions(explanation):
    return explanation["genes"].set_index(["state", "gene"])["contribution"].sort_index()


def composition_only_example():
    """A known mixture change with exactly unchanged local expression."""
    matrices, records = [], []
    for condition, n_late in (("control", 12), ("treated", 36)):
        for sample_index in range(4):
            sample = f"{condition}_{sample_index + 1}"
            local = np.r_[np.zeros(48 - n_late), np.ones(n_late)]
            matrices.append(np.column_stack([1 + 4 * local, 2 + 2 * local, 3 + 0 * local]))
            records.extend({
                "condition": condition,
                "sample": sample,
                "state": "Composition example",
                "region": "late" if value else "early",
            } for value in local)
    obs = pd.DataFrame.from_records(records)
    obs.index = [f"mixture_cell_{i}" for i in range(len(obs))]
    adata = AnnData(
        X=np.vstack(matrices),
        obs=obs,
        var=pd.DataFrame(index=["Substate_A", "Substate_B", "Constant"]),
    )
    adata.obsm["X_linear"] = adata.X.copy()
    adata.varm["identity_loadings"] = np.eye(adata.n_vars)
    result = _explain(
        adata,
        {"Substate program": ["Substate_A", "Substate_B"]},
        rep="X_linear",
        loadings_key="identity_loadings",
        center="mean",
        n_boot=0,
    )
    local_delta = result["local_genes"]["expression_delta"].to_numpy(dtype=float)
    return adata, result, {
        "data_kind": "deterministic mixture counterexample",
        "source_estimator": "mean of sample means",
        "reference_late_fraction": 0.25,
        "comparison_late_fraction": 0.75,
        "global_delta_norm": float(result["summary"].iloc[0]["delta_norm"]),
        "maximum_absolute_local_expression_change": float(np.max(np.abs(local_delta))),
        "interpretation": "the broad-state displacement comes entirely from supplied-substate composition",
        "biological_validation": False,
    }


def _local_response_figure(explanation, feature_truth):
    local = explanation["local_genes"].merge(
        feature_truth[["gene", "program"]], on="gene", how="left"
    )
    shown_programs = ["Response", "LocalResponse", "Identity"]
    selected = local[local.program.isin(shown_programs)].copy()
    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    colors = {"early": "#4477AA", "late": "#D67732"}
    for j, region in enumerate(("early", "late")):
        for i, program in enumerate(shown_programs):
            values = selected.loc[
                selected.local_state.eq(region) & selected.program.eq(program),
                "expression_delta",
            ].to_numpy(dtype=float)
            x = i + (-0.15 if j == 0 else 0.15)
            jitter = np.linspace(-0.035, 0.035, len(values))
            ax.scatter(x + jitter, values, color=colors[region], s=34, alpha=0.75,
                       label=region.capitalize() if i == 0 else None)
            if len(values):
                ax.plot([x - 0.07, x + 0.07], [values.mean()] * 2,
                        color=colors[region], linewidth=2.5)
    ax.axhline(0, color="#7B8792", linewidth=0.8)
    ax.set_xticks(range(len(shown_programs)))
    ax.set_xticklabels(["Uniform response", "Opposing local response", "Baseline identity"])
    ax.set_ylabel("Treated − control expression\n(arbitrary input units)")
    ax.set_title("Expression response within supplied baseline regions")
    ax.legend(frameon=False, title="Baseline region", loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.01, 0.01,
             "Synthetic input · dots are genes; short lines are program means · no causal inference",
             fontsize=8, color="#5F6C78")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return sg.pl.PlotResult(
        fig=fig,
        axes={"local_expression": ax},
        data={"local_gene_contrasts": selected},
        metadata={
            "plot": "synthetic local expression contrasts",
            "metric": "expression_delta from explain_shift local_genes",
            "uncertainty": "none; dots show individual member genes, not independent replicates",
            "local_direction": "expression signs, not each local shift's own direction",
            "biological_validation": False,
        },
    )


def _write_tables(result, output_dir, prefix):
    paths = {}
    for name, table in result.items():
        if not isinstance(table, pd.DataFrame):
            continue
        filename = f"{prefix}_{name}.csv"
        table.to_csv(output_dir / filename, index=False)
        paths[name] = {"file": filename, "rows": len(table)}
    return paths


def run_gallery(output_dir, formats=("png", "svg", "pdf"), *, seed=37, n_boot=80):
    """Run a small scientific example and export every plotted source table."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata, programs, cell_truth, feature_truth = make_example_data(seed=seed)
    explanation = _explain(adata, programs, n_boot=n_boot, seed=seed)
    source_tables = _write_tables(explanation, output_dir, "explanation")
    adata.obs.to_csv(output_dir / "synthetic_cells.csv")
    pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names).to_csv(
        output_dir / "synthetic_expression.csv"
    )
    cell_truth.to_csv(output_dir / "synthetic_cell_truth.csv", index=False)
    feature_truth.to_csv(output_dir / "synthetic_feature_truth.csv", index=False)

    main_contributions = _gene_contributions(explanation)
    rng = np.random.RandomState(seed + 100)
    rotation, _ = np.linalg.qr(rng.normal(size=(8, 8)))
    rotated = adata.copy()
    rotated.obsm["X_rotated"] = adata.obsm["X_pca"] @ rotation
    rotated.varm["PCs_rotated"] = adata.varm["PCs"] @ rotation
    rotated_explanation = _explain(
        rotated, programs, n_boot=0, seed=seed,
        rep="X_rotated", loadings_key="PCs_rotated",
    )
    rotation_error = float(np.max(np.abs(
        _gene_contributions(rotated_explanation) - main_contributions
    )))

    chosen = np.flatnonzero(adata.obs["sample"].eq("treated_1").to_numpy())
    indices = np.r_[np.arange(adata.n_obs), np.tile(chosen, 4)]
    repeated_obs = adata.obs.iloc[indices].copy()
    repeated_obs.index = [f"repeat_{i}" for i in range(len(indices))]
    repeated = AnnData(X=adata.X[indices].copy(), obs=repeated_obs, var=adata.var.copy())
    repeated.obsm["X_pca"] = adata.obsm["X_pca"][indices].copy()
    repeated.varm["PCs"] = adata.varm["PCs"].copy()
    repeated_explanation = _explain(repeated, programs, n_boot=0, seed=seed)
    duplication_error = float(np.max(np.abs(
        _gene_contributions(repeated_explanation) - main_contributions
    )))

    mixture_data, mixture_result, mixture_metrics = composition_only_example()
    mixture_tables = _write_tables(mixture_result, output_dir, "composition_only")
    mixture_data.obs.to_csv(output_dir / "composition_only_cells.csv")
    pd.DataFrame(mixture_data.X, index=mixture_data.obs_names,
                 columns=mixture_data.var_names).to_csv(
                     output_dir / "composition_only_expression.csv")

    summary = explanation["summary"]
    reconstruction_error = float(summary["reconstruction_error"].abs().max())
    local = explanation["local_genes"]
    planted_local = local[local.gene.str.startswith("LocalResponse_")]
    local_signs = planted_local.groupby("local_state", observed=True)["expression_delta"].mean()
    checks = {
        "affine_reconstruction": reconstruction_error < 1e-5,
        "rotation_invariance": rotation_error < 1e-5,
        "sample_duplication_invariance": duplication_error < 1e-5,
        "opposing_local_expression_signs": bool(local_signs["early"] < 0 < local_signs["late"]),
        "composition_only_local_expression_unchanged": (
            mixture_metrics["maximum_absolute_local_expression_change"] < 1e-10
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"Synthetic gallery checks failed: {checks}")

    exports = {}
    condition_palette = {"control": "#4477AA", "treated": "#D67732"}
    with sg.pl.style_context("paper"):
        figures = [
            ("gene_shift_curves", sg.pl.shift_curves(
                adata,
                var_names=["Response_1", "LocalResponse_1", "Identity_1"],
                state="Simulated state",
                basis="display",
                palette=condition_palette,
                panels=("curve", "expression", "contribution"),
                title="ScGeo shift explanations · simulated expression",
                return_result=True,
                show=False,
            )),
            ("program_shift_curves", sg.pl.shift_curves(
                adata,
                program_names=["Response", "LocalResponse"],
                state="Simulated state",
                basis="display",
                palette=condition_palette,
                panels=("curve",),
                title="Program responses · simulated input",
                figsize=(6.4, 7.2),
                return_result=True,
                show=False,
            )),
            ("state_shift_embedding", sg.pl.shift_embedding(
                adata,
                states=["Simulated state"],
                basis="display",
                color="condition",
                palette=condition_palette,
                title="State contrast on fitted PCA display coordinates",
                return_result=True,
                show=False,
            )),
            ("local_expression_contrasts", _local_response_figure(explanation, feature_truth)),
        ]
        for name, result in figures:
            result.metadata["data_provenance"] = dict(adata.uns["synthetic_provenance"])
            files = result.save(output_dir / name, formats=formats, dpi=180)
            exports[name] = str(Path(files["metadata"]).name)
            plt.close(result.fig)

    benchmark = {
        "schema_version": 1,
        "data_kind": "fully synthetic; no private inputs or downloads",
        "seed": int(seed),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "samples_per_condition": 6,
        "n_boot": int(n_boot),
        "input_sha256": _input_hash(adata),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_params": explanation["source_params"],
        "mapping": {
            **explanation["mapping"],
            "intercept": np.asarray(explanation["mapping"]["intercept"]).tolist(),
        },
        "checks": checks,
        "metrics": {
            "maximum_reconstruction_error": reconstruction_error,
            "maximum_gene_contribution_rotation_error": rotation_error,
            "maximum_gene_contribution_sample_duplication_error": duplication_error,
            "local_program_expression_contrasts": {str(k): float(v) for k, v in local_signs.items()},
        },
        "summary": json.loads(summary.to_json(orient="records")),
        "local_summary": json.loads(explanation["local_summary"].to_json(orient="records")),
        "composition_only": mixture_metrics,
        "distributional_counterexample": distributional_counterexample(),
        "explanation_source_tables": source_tables,
        "composition_source_tables": mixture_tables,
        "plot_manifests": exports,
        "interpretation": {
            "descriptive_accounting_validated": True,
            "biological_driver_validation": False,
            "causal_recovery_validation": False,
            "scope": "fixed affine PCA mapping and stated synthetic examples",
            "curves": "expression along a geometric projection, not inferred time or RNA velocity",
            "arrows": "schematic condition connectors, not reconstructed cell paths",
        },
        "reproduce": (
            "python examples/explain_shift_gallery.py "
            f"--seed {seed} --n-boot {n_boot} --output-dir NEW_OUTPUT_DIRECTORY"
        ),
    }
    (output_dir / "benchmark_summary.json").write_text(
        json.dumps(benchmark, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return benchmark


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/explain_shift"))
    parser.add_argument("--formats", nargs="+", choices=("png", "svg", "pdf"),
                        default=("png", "svg", "pdf"))
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument("--n-boot", type=int, default=80)
    args = parser.parse_args()
    result = run_gallery(args.output_dir, args.formats, seed=args.seed, n_boot=args.n_boot)
    print(json.dumps({"checks": result["checks"], "metrics": result["metrics"],
                      "plot_manifests": result["plot_manifests"]}, indent=2))


if __name__ == "__main__":
    main()
