#!/usr/bin/env python3
"""Assemble frozen GSE249479 Dataset B evidence into manuscript artifacts.

This module performs plotting and evidence collation only. It does not run
normalization, feature selection, representation learning, ScGeo, Augur, or
any other numerical model.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/public_validation/gse249479_dataset_b"
OUT = BASE / "manuscript"
SCGEO = ROOT.parent / "scgeo"
SCGEO_COMMIT = "9a0ed16cbaa57f935f9c9bc87d1643a25b51012c"
INFERENCE_STATUS = "descriptive_only"

EXPECTED_HASHES = {
    "figure_sources/03d_display_umap_coordinates.csv.gz": "563c1ca84d582d2d0dbb331ec811c0473bfa4eb2e43f04299638559785d54cff",
    "audit/01_qc_distribution_summary.csv": "ac0552a3d2864a03e8f79c34d4585bb43817376e364ea63daffd7616a8a1e616",
    "audit/02_annotation_confidence_table.csv": "d5e1615feeb5795556759dd58d527c4a7331357175d9a4968d0453c99dd966ff",
    "audit/02_annotation_summary.json": "a614f02de2068eb564a0b19c20deffe4bece04b4ff7d512f90e8a57c822f5547",
    "figure_sources/03d_representation_quality_metrics.csv": "6d5885fabc758d76326fbdad68d92d291db4329e1580562ee7c578b770dee0e3",
    "figure_sources/03d_rare_state_sensitivity.csv": "b1dca72b9e0a6cdc5a486ba5e558f64d1d4c3a62bc85685dcdf0ccebefaacdec",
    "scgeo/04_tnf_primary_consensus.csv": "ad43151a1d6124cf5383ebca606714b7ac3a9a048269cff43e69eceb1a255357",
    "scgeo/04_lps_primary_consensus.csv": "556aa43fed6ef20e874bec3780038287a7d9a580cfd434564e963c3cde68ee07",
    "scgeo/04_tnf_detailed_representation_state_evidence.csv": "b0b2227c8ff8e18730421905117e90807cc210a20ffc9bd5a8e033909868edb0",
    "scgeo/04_lps_detailed_representation_state_evidence.csv": "42fdade846a0e2665f48311bd4b38aa0f7adb4ea2b2298e7dc37e26a2622131f",
    "scgeo/04_local_geometry_representation_summary.csv": "86607d9fd5702ce883d17a31c4889ed8b6721c78129fbf229de86c47e01321c1",
    "controls/04_exclude_diffusion_sensitivity.csv": "cf47c4d718a18ee64d9280d50ee7de843d3fecd3c33f61076565ae2a46bad73f",
    "controls/04_tnf_leave_one_representation.csv": "9527edb5eb402e0326ae25598fc2ad41ddda214adf71698377f67eb9bb6771cc",
    "controls/04_lps_leave_one_representation.csv": "0ba17a969aa55c119fd2750ff9a90db9e8c457e3094e44cf826ee791811e2f44",
    "controls/04_balanced_state_displacement_subsampling.csv": "13cd9f7ceb53f58248e2552d1c6b70de929ad2c318f1c9c9929a76de42ec76bc",
    "controls/04_shuffled_condition_computational_null.csv": "70f467b3185fcc81fb649cc5f0f4b586ce282f637a283a1b9eba09b29b6c1db2",
    "comparator/06_scgeo_augur_full_comparator.csv": "c2b8c1782900ef684cfb67f97cbd6fca3e2f473d15ca27444afd49ed0f9ab004",
    "comparator/06_scgeo_augur_rank_correlations.csv": "2cc77335d007cf57dc480f3e5bdfb2b78dbc7962f5887c75f31b56f296c4d5be",
    "comparator/06_official_python_implementation_sensitivity.csv": "ced5dbd8887b3f35fe6890e27fc8626f2d467d7df5fd3cb5b111de0f443cd56e",
    "figure_sources/04_signature_displacement_associations.csv": "1d0d66f66c12fcb0177639b8ab82e8b09cbb0d85468679f886c8a63af8bf39f5",
    "signatures/signature_definitions.csv": "b8c9b2d874c27288f9f59609e7266d9af0865ba55674c1118204870f0053065c",
    "audit/replication_reassessment.json": "ce1afb1839dc916de29909bfedff781fc37f97fa73c3f9289d853a067d6242e8",
    "audit/06_comparator_documentation_correction.json": "8e1d432f55f37e66323a10bcf6f45ed0923901957dcf6479643bd01c0a937e08",
}

STATE_ORDER = [
    "HSC/quiescent", "Activated HSC", "MPP/progenitor", "Myeloid",
    "Megakaryocyte/erythroid", "Lymphoid", "Ambiguous HSPC",
]
STATE_MAP = {
    "HSC_quiescent": "HSC/quiescent", "activated_HSC": "Activated HSC",
    "MPP_progenitor": "MPP/progenitor", "myeloid": "Myeloid",
    "megakaryocyte_erythroid": "Megakaryocyte/erythroid",
    "lymphoid": "Lymphoid", "ambiguous_HSPC": "Ambiguous HSPC",
}
STATE_COLORS = {
    "HSC/quiescent": "#4C78A8", "Activated HSC": "#F58518",
    "MPP/progenitor": "#ECA82C", "Myeloid": "#54A24B",
    "Megakaryocyte/erythroid": "#B279A2", "Lymphoid": "#E45756",
    "Ambiguous HSPC": "#9D9D9D",
}
CONDITION_COLORS = {"PBS": "#777777", "TNF": "#D1495B", "LPS": "#3E7CB1"}
REP_ORDER = ["X_pca30", "X_pca50", "X_scvi", "X_pca20", "X_diffmap"]
REP_LABELS = ["PCA30\nprimary", "PCA50\nprimary", "scVI\nprimary", "PCA20\nsensitivity", "Diffusion\nexploratory"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_output(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def validate_frozen_sources() -> dict[str, str]:
    observed: dict[str, str] = {}
    for relative, expected in EXPECTED_HASHES.items():
        path = BASE / relative
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen manuscript source: {path}")
        observed[relative] = sha256(path)
        if observed[relative] != expected:
            raise RuntimeError(f"Frozen source checksum mismatch for {relative}")
    if git_output(SCGEO, "rev-parse", "HEAD") != SCGEO_COMMIT:
        raise RuntimeError("Frozen ScGeo commit mismatch")
    if git_output(SCGEO, "status", "--short"):
        raise RuntimeError("Frozen ScGeo checkout is not clean")
    return observed


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".csv.gz" if compression == "gzip" else ".csv"
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}{suffix}")
    frame.to_csv(tmp, index=False, compression=compression)
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text.rstrip() + "\n", encoding="utf-8")
    os.replace(tmp, path)


def save_figure(fig: plt.Figure, stem: str) -> None:
    directory = OUT / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    metadata = {"Creator": "GSE249479 frozen-artifact manuscript assembly", "Date": None}
    for extension in ["png", "svg"]:
        fig.savefig(directory / f"{stem}.{extension}", dpi=300, bbox_inches="tight", metadata=metadata)
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


def read_sources() -> dict[str, Any]:
    return {
        "umap": pd.read_csv(BASE / "figure_sources/03d_display_umap_coordinates.csv.gz"),
        "qc": pd.read_csv(BASE / "audit/01_qc_distribution_summary.csv"),
        "annotation": pd.read_csv(BASE / "audit/02_annotation_confidence_table.csv"),
        "annotation_summary": json.loads((BASE / "audit/02_annotation_summary.json").read_text()),
        "quality": pd.read_csv(BASE / "figure_sources/03d_representation_quality_metrics.csv"),
        "rare": pd.read_csv(BASE / "figure_sources/03d_rare_state_sensitivity.csv"),
        "tnf_primary": pd.read_csv(BASE / "scgeo/04_tnf_primary_consensus.csv"),
        "lps_primary": pd.read_csv(BASE / "scgeo/04_lps_primary_consensus.csv"),
        "tnf_rep": pd.read_csv(BASE / "scgeo/04_tnf_detailed_representation_state_evidence.csv"),
        "lps_rep": pd.read_csv(BASE / "scgeo/04_lps_detailed_representation_state_evidence.csv"),
        "geometry": pd.read_csv(BASE / "scgeo/04_local_geometry_representation_summary.csv"),
        "exclude_diffusion": pd.read_csv(BASE / "controls/04_exclude_diffusion_sensitivity.csv"),
        "loo": pd.concat([
            pd.read_csv(BASE / "controls/04_tnf_leave_one_representation.csv"),
            pd.read_csv(BASE / "controls/04_lps_leave_one_representation.csv"),
        ], ignore_index=True),
        "balanced": pd.read_csv(BASE / "controls/04_balanced_state_displacement_subsampling.csv"),
        "shuffled": pd.read_csv(BASE / "controls/04_shuffled_condition_computational_null.csv"),
        "comparator": pd.read_csv(BASE / "comparator/06_scgeo_augur_full_comparator.csv"),
        "correlations": pd.read_csv(BASE / "comparator/06_scgeo_augur_rank_correlations.csv"),
        "implementation": pd.read_csv(BASE / "comparator/06_official_python_implementation_sensitivity.csv"),
        "signature_assoc": pd.read_csv(BASE / "figure_sources/04_signature_displacement_associations.csv"),
        "signature_defs": pd.read_csv(BASE / "signatures/signature_definitions.csv"),
        "replication": json.loads((BASE / "audit/replication_reassessment.json").read_text()),
    }


def main_figure(data: dict[str, Any]) -> None:
    umap = data["umap"].copy()
    umap["state"] = umap["marker_inferred_label"].map(STATE_MAP)
    comparator = data["comparator"].copy()
    reps = pd.concat([data["tnf_rep"], data["lps_rep"]], ignore_index=True)
    reps = reps[reps["hierarchy"].eq("detailed") & reps["representation"].isin(REP_ORDER)]

    atomic_csv(umap[["obs_name", "condition", "state", "X_umap_pca50_1", "X_umap_pca50_2"]], OUT / "figure_sources/main_A_display_umap.csv.gz", "gzip")
    atomic_csv(comparator[["contrast", "state", "scgeo_normalized_displacement", "representation_consensus_status", "official_augur_auc", "abundance_change", "energy_distance", "mmd", "sliced_wasserstein"]], OUT / "figure_sources/main_B_D_E_evidence.csv")
    atomic_csv(reps[["contrast", "state", "representation", "normalized_effect", "status", "representation_role"]], OUT / "figure_sources/main_C_representation_robustness.csv")

    fig = plt.figure(figsize=(16, 20), constrained_layout=True)
    outer = fig.add_gridspec(4, 2, height_ratios=[1.05, 1.2, 1.15, 0.32])
    sub_a = outer[0, :].subgridspec(1, 2, wspace=0.05)
    ax_state, ax_condition = fig.add_subplot(sub_a[0, 0]), fig.add_subplot(sub_a[0, 1])
    for state in STATE_ORDER:
        frame = umap[umap["state"].eq(state)]
        ax_state.scatter(frame.X_umap_pca50_1, frame.X_umap_pca50_2, s=1, alpha=.35, color=STATE_COLORS[state], rasterized=True, label=state)
    for condition in ["PBS", "TNF", "LPS"]:
        frame = umap[umap["condition"].eq(condition)]
        ax_condition.scatter(frame.X_umap_pca50_1, frame.X_umap_pca50_2, s=1, alpha=.30, color=CONDITION_COLORS[condition], rasterized=True, label=condition)
    for ax, title in [(ax_state, "Conservative marker-inferred states"), (ax_condition, "Condition distributions")]:
        ax.set(title=title, xlabel="display UMAP 1", ylabel="display UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    ax_state.legend(loc="upper left", bbox_to_anchor=(1.00, 1.0), frameon=False, markerscale=4, fontsize=8)
    ax_condition.legend(frameon=False, markerscale=4)
    panel_label(ax_state, "A")
    ax_condition.text(.99, .01, "PCA50 UMAP: visualization only\nnot used for quantitative geometry", transform=ax_condition.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")

    ax_b = fig.add_subplot(outer[1, 0]); panel_label(ax_b, "B")
    offsets = {"PBS_vs_TNF": -.12, "PBS_vs_LPS": .12}; markers = {"PBS_vs_TNF": "o", "PBS_vs_LPS": "s"}
    for contrast in ["PBS_vs_TNF", "PBS_vs_LPS"]:
        frame = comparator[comparator.contrast.eq(contrast)].set_index("state").loc[STATE_ORDER].reset_index()
        y = np.arange(len(STATE_ORDER)) + offsets[contrast]
        stable = frame.representation_consensus_status.eq("stable_effect")
        color = CONDITION_COLORS[contrast.removeprefix("PBS_vs_")]
        ax_b.scatter(frame.scgeo_normalized_displacement, y, marker=markers[contrast], s=60, facecolors=np.where(stable, color, "white"), edgecolors=color, label=contrast.replace("PBS_vs_", "PBS → "))
    ax_b.set_yticks(np.arange(len(STATE_ORDER)), STATE_ORDER); ax_b.invert_yaxis(); ax_b.set_xlabel("Primary-consensus normalized displacement"); ax_b.set_title("State displacement and ScGeo consensus status")
    ax_b.grid(axis="x", color="#EEEEEE"); ax_b.legend(frameon=False, title="Filled = stable_effect")
    rho_tnf = data["tnf_primary"].pairwise_spearman_median.iloc[0]; rho_lps = data["lps_primary"].pairwise_spearman_median.iloc[0]
    ax_b.text(.99, .02, f"Primary rank agreement: TNF ρ={rho_tnf:.3f}; LPS ρ={rho_lps:.3f}", transform=ax_b.transAxes, ha="right", fontsize=8)

    ax_c = fig.add_subplot(outer[1, 1]); panel_label(ax_c, "C")
    row_order = [(contrast, state) for contrast in ["PBS_vs_TNF", "PBS_vs_LPS"] for state in STATE_ORDER]
    matrix = reps.pivot(index=["contrast", "state"], columns="representation", values="normalized_effect").reindex(index=row_order, columns=REP_ORDER)
    im = ax_c.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=np.nanmax(matrix.to_numpy()))
    ax_c.set_xticks(range(len(REP_ORDER)), REP_LABELS); ax_c.set_yticks(range(len(row_order)), [f"{c.split('_')[-1]}: {s}" for c, s in row_order], fontsize=7)
    ax_c.axvline(2.5, color="white", lw=2); ax_c.set_title("Representation robustness: primary and sensitivity views")
    fig.colorbar(im, ax=ax_c, fraction=.025, label="normalized displacement")

    sub_d = outer[2, 0].subgridspec(1, 3, wspace=.08)
    metrics = [("abundance_change", "Abundance Δ", "coolwarm"), ("scgeo_normalized_displacement", "Displacement", "YlGnBu"), ("energy_distance", "Distribution\nenergy", "magma")]
    rows = comparator.set_index(["contrast", "state"]).reindex(row_order)
    axes_d = []
    for index, (column, title, cmap) in enumerate(metrics):
        ax = fig.add_subplot(sub_d[0, index]); axes_d.append(ax)
        values = rows[column].to_numpy()[:, None]
        limit = np.nanmax(np.abs(values)) if column == "abundance_change" else None
        image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=-limit if limit else None, vmax=limit if limit else None)
        ax.set_xticks([0], [title]); ax.set_yticks(range(len(row_order)), [f"{c.split('_')[-1]}: {s}" for c, s in row_order] if index == 0 else [], fontsize=7)
        for y, value in enumerate(values[:, 0]): ax.text(0, y, f"{value:+.2f}" if column == "abundance_change" else f"{value:.2f}", ha="center", va="center", fontsize=6, color="black" if abs(value) < .7 * (limit or np.nanmax(values)) else "white")
        fig.colorbar(image, ax=ax, fraction=.08)
    panel_label(axes_d[0], "D"); axes_d[1].set_title("Separate descriptive axes; column-specific scales; no composite", pad=25)

    ax_e = fig.add_subplot(outer[2, 1]); panel_label(ax_e, "E")
    short_state = {"HSC/quiescent":"HSC/q", "Activated HSC":"Act HSC", "MPP/progenitor":"MPP", "Myeloid":"Myeloid", "Megakaryocyte/erythroid":"Mk/Ery", "Lymphoid":"Lymphoid", "Ambiguous HSPC":"Ambig"}
    label_offsets = {
        ("PBS_vs_TNF", "HSC/quiescent"):(-24, 7), ("PBS_vs_TNF", "Activated HSC"):(3, 7),
        ("PBS_vs_TNF", "MPP/progenitor"):(3, 7), ("PBS_vs_TNF", "Myeloid"):(3, -12),
        ("PBS_vs_TNF", "Megakaryocyte/erythroid"):(-32, -13), ("PBS_vs_TNF", "Lymphoid"):(3, 6),
        ("PBS_vs_TNF", "Ambiguous HSPC"):(3, 7), ("PBS_vs_LPS", "HSC/quiescent"):(-28, -12),
        ("PBS_vs_LPS", "Activated HSC"):(3, -12), ("PBS_vs_LPS", "MPP/progenitor"):(3, 7),
        ("PBS_vs_LPS", "Myeloid"):(3, -12), ("PBS_vs_LPS", "Megakaryocyte/erythroid"):(-36, -12),
        ("PBS_vs_LPS", "Lymphoid"):(3, 6), ("PBS_vs_LPS", "Ambiguous HSPC"):(3, 7),
    }
    for contrast, color, marker in [("PBS_vs_TNF", CONDITION_COLORS["TNF"], "o"), ("PBS_vs_LPS", CONDITION_COLORS["LPS"], "s")]:
        frame = comparator[comparator.contrast.eq(contrast)]
        ax_e.scatter(frame.official_augur_auc, frame.scgeo_normalized_displacement, color=color, marker=marker, s=55, label=contrast.replace("PBS_vs_", "PBS → "))
        for row in frame.itertuples():
            ax_e.annotate(short_state[row.state], (row.official_augur_auc, row.scgeo_normalized_displacement), fontsize=6, xytext=label_offsets[(contrast, row.state)], textcoords="offset points")
    ax_e.set(xlabel="Official R Augur AUC (treatment separability)", ylabel="ScGeo primary normalized displacement", title="Official comparator versus representation-stable displacement\nDifferent estimands; n=7 marker-inferred states per contrast")
    ax_e.grid(color="#EEEEEE"); ax_e.legend(frameon=False)

    ax_f = fig.add_subplot(outer[3, :]); panel_label(ax_f, "F"); ax_f.axis("off")
    ax_f.text(.5, .55, "LIMITATION  •  descriptive_only  •  no independent biological replicate identity  •  cells, conditions, libraries, and SouporCell clades are not replicates", ha="center", va="center", fontsize=11, fontweight="bold", bbox=dict(boxstyle="round,pad=.55", fc="#F5F5F5", ec="#666666"))
    fig.suptitle("Dataset B (GSE249479): descriptive treatment geometry across conservative HSPC states", fontsize=16, fontweight="bold")
    save_figure(fig, "dataset_b_main_figure")


def supplementary_figures(data: dict[str, Any]) -> None:
    # S1: QC and annotation confidence.
    qc = data["qc"]; annotation = data["annotation"].copy(); annotation["state"] = annotation.marker_inferred_label.map(STATE_MAP)
    atomic_csv(qc, OUT / "figure_sources/supp1_qc_distribution_summary.csv"); atomic_csv(annotation, OUT / "figure_sources/supp1_annotation_confidence.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for metric, color in [("total_counts", "#4C78A8"), ("detected_genes", "#59A14F"), ("pct_counts_mito", "#E15759")]:
        frame = qc[(qc.condition == "ALL") & (qc.metric == metric)].iloc[0]
        axes[0].errorbar(metric.replace("pct_counts_", "% ").replace("_", "\n"), frame["median"], yerr=[[frame["median"]-frame.q05], [frame.q95-frame["median"]]], fmt="o", color=color, capsize=4)
    axes[0].set_title("QC median and 5th–95th percentiles"); axes[0].set_yscale("symlog"); axes[0].set_ylabel("metric-specific scale")
    confidence = annotation.pivot(index="state", columns="annotation_confidence", values="n_cells").fillna(0).reindex(STATE_ORDER)
    confidence.plot(kind="barh", stacked=True, ax=axes[1], color={"low":"#F2CF5B", "moderate":"#59A14F", "ambiguous":"#9D9D9D"}); axes[1].set(title="Marker-label confidence", xlabel="cells", ylabel="")
    counts = pd.Series(data["annotation_summary"]["confidence_counts"]); axes[2].bar(counts.index, counts.values, color=["#9D9D9D", "#F2CF5B", "#59A14F"]); axes[2].set(title="Overall annotation confidence", ylabel="cells"); axes[2].tick_params(axis="x", rotation=30)
    fig.suptitle("Supplementary 1. Sparse QC and conservative annotation confidence")
    save_figure(fig, "dataset_b_supp1_qc_annotation")

    # S2: representation quality.
    quality = data["quality"].copy(); atomic_csv(quality, OUT / "figure_sources/supp2_representation_quality.csv")
    metrics = ["marker_label_neighbor_preservation", "souporcell_clade_neighbor_preservation_mapped", "neighborhood_overlap_vs_pca50", "local_distance_spearman", "largest_graph_component_fraction"]
    fig, ax = plt.subplots(figsize=(10, 4.8)); matrix = quality.set_index("representation").reindex(REP_ORDER)[metrics]
    im = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1); ax.set_xticks(range(len(metrics)), [m.replace("_", "\n") for m in metrics], fontsize=8); ax.set_yticks(range(len(REP_ORDER)), REP_LABELS); fig.colorbar(im, ax=ax, label="metric value")
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]): ax.text(x, y, f"{matrix.iloc[y,x]:.2f}", ha="center", va="center", fontsize=7, color="white" if matrix.iloc[y,x] < .45 else "black")
    ax.set_title("Supplementary 2. Representation-quality metrics (condition mixing is descriptive, not an optimization target)")
    save_figure(fig, "dataset_b_supp2_representation_quality")

    # S3: diffusion distortion and leave-one-representation sensitivity.
    geometry = data["geometry"].copy(); exclude = data["exclude_diffusion"].copy(); loo = data["loo"].copy()
    loo_summary = loo.groupby(["contrast", "state"], observed=False).relative_change.max().reset_index(name="max_leave_one_relative_change")
    atomic_csv(geometry, OUT / "figure_sources/supp3_geometry_distortion.csv"); atomic_csv(exclude, OUT / "figure_sources/supp3_exclude_diffusion.csv"); atomic_csv(loo_summary, OUT / "figure_sources/supp3_leave_one_representation.csv")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    axes[0].bar(geometry.rep, geometry.median_global_distortion, color=["#4C78A8" if r != "X_diffmap" else "#E15759" for r in geometry.rep]); axes[0].tick_params(axis="x", rotation=35); axes[0].set(title="Global distortion", ylabel="median distortion")
    for contrast, color in [("PBS_vs_TNF", CONDITION_COLORS["TNF"]), ("PBS_vs_LPS", CONDITION_COLORS["LPS"])]:
        frame=exclude[exclude.contrast.eq(contrast)]; axes[1].scatter(frame.absolute_change_after_excluding_diffusion, frame.state, label=contrast, color=color)
    axes[1].set(title="Effect of excluding diffusion map", xlabel="absolute change in median normalized effect"); axes[1].legend(frameon=False)
    pivot=loo_summary.pivot(index="state",columns="contrast",values="max_leave_one_relative_change").reindex(STATE_ORDER); pivot.plot(kind="barh",ax=axes[2],color=[CONDITION_COLORS["LPS"],CONDITION_COLORS["TNF"]]); axes[2].set(title="Leave-one-representation sensitivity",xlabel="maximum relative change",ylabel="")
    fig.suptitle("Supplementary 3. Diffusion-map sensitivity and robustness of primary conclusions")
    save_figure(fig, "dataset_b_supp3_diffusion_leave_one")

    # S4: cell-balanced and shuffled-label controls.
    balanced = data["balanced"]; shuffled = data["shuffled"]
    balanced_summary = balanced.groupby(["contrast","state"],observed=False).normalized_effect.agg(median="median",q025=lambda x:x.quantile(.025),q975=lambda x:x.quantile(.975)).reset_index()
    shuffled_summary = shuffled.groupby(["contrast","state"],observed=False).normalized_effect.agg(median="median",q025=lambda x:x.quantile(.025),q975=lambda x:x.quantile(.975)).reset_index()
    atomic_csv(balanced_summary, OUT / "figure_sources/supp4_balanced_cell_control.csv"); atomic_csv(shuffled_summary, OUT / "figure_sources/supp4_shuffled_label_control.csv")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    for ax, frame, title in [(axes[0],balanced_summary,"Repeated balanced cell subsampling"),(axes[1],shuffled_summary,"Shuffled cell-label computational null")]:
        for contrast,color,offset in [("PBS_vs_TNF",CONDITION_COLORS["TNF"],-.1),("PBS_vs_LPS",CONDITION_COLORS["LPS"],.1)]:
            part=frame[frame.contrast.eq(contrast)].set_index("state").reindex(STATE_ORDER); y=np.arange(len(STATE_ORDER))+offset
            ax.errorbar(part["median"],y,xerr=[part["median"]-part.q025,part.q975-part["median"]],fmt="o",color=color,capsize=2,label=contrast)
        ax.set_yticks(range(len(STATE_ORDER)),STATE_ORDER); ax.invert_yaxis(); ax.set(title=title,xlabel="normalized displacement"); ax.legend(frameon=False)
    fig.suptitle("Supplementary 4. Cell-level computational controls only; no biological uncertainty or p-values")
    save_figure(fig, "dataset_b_supp4_computational_controls")

    # S5: official versus Python implementation sensitivity.
    implementation=data["implementation"].copy(); atomic_csv(implementation,OUT/"figure_sources/supp5_official_python_augur.csv")
    fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    for contrast,color,marker in [("PBS_vs_TNF",CONDITION_COLORS["TNF"],"o"),("PBS_vs_LPS",CONDITION_COLORS["LPS"],"s")]:
        frame=implementation[implementation.contrast.eq(contrast)]; axes[0].scatter(frame.official_augur_auc,frame.python_augur_inspired_auc,color=color,marker=marker,label=contrast)
        for row in frame.itertuples(): axes[0].annotate(row.state,(row.official_augur_auc,row.python_augur_inspired_auc),fontsize=6)
    axes[0].plot([.48,.72],[.48,.72],ls="--",color="#777777"); axes[0].set(xlabel="Official R Augur AUC",ylabel="Augur-inspired Python AUC",title="AUC implementation sensitivity"); axes[0].legend(frameon=False)
    ordered=implementation.sort_values(["contrast","jaccard_overlap"]); axes[1].barh(range(len(ordered)),ordered.jaccard_overlap,color="#76B7B2"); axes[1].set_yticks(range(len(ordered)),[f"{r.contrast.split('_')[-1]}: {r.state}" for r in ordered.itertuples()]); axes[1].set(xlabel="feature-set Jaccard",title="Official LOESS residuals vs Python raw variance",xlim=(0,1))
    fig.suptitle("Supplementary 5. Official Augur is primary; Python is implementation sensitivity only")
    save_figure(fig,"dataset_b_supp5_official_python_augur")

    # S6: signature associations and HSC-I sparsity.
    associations=data["signature_assoc"].copy(); definitions=data["signature_defs"].copy(); atomic_csv(associations,OUT/"figure_sources/supp6_signature_associations.csv"); atomic_csv(definitions,OUT/"figure_sources/supp6_signature_definitions.csv")
    pivot=associations.pivot(index="signature",columns="contrast",values="spearman_across_states")
    fig,axes=plt.subplots(1,2,figsize=(15,6),gridspec_kw={"width_ratios":[1.4,1]},constrained_layout=True)
    im=axes[0].imshow(pivot,aspect="auto",cmap="coolwarm",vmin=-1,vmax=1); axes[0].set_xticks(range(len(pivot.columns)),[c.replace("PBS_vs_","PBS → ") for c in pivot.columns]); axes[0].set_yticks(range(len(pivot.index)),pivot.index); fig.colorbar(im,ax=axes[0],label="Spearman across 7 states")
    for y in range(pivot.shape[0]):
        for x in range(pivot.shape[1]): axes[0].text(x,y,f"{pivot.iloc[y,x]:.2f}",ha="center",va="center",fontsize=7)
    hsc_i=definitions[definitions.signature.eq("study_hsc_i_top50_subset")].iloc[0]; axes[1].axis("off"); axes[1].text(.03,.92,"HSC-I signature limitation",fontsize=13,fontweight="bold",va="top"); axes[1].text(.03,.78,f"Available genes: {int(hsc_i.n_available)} / {int(hsc_i.n_requested)}",fontsize=12); axes[1].text(.03,.64,f"Missing: {hsc_i.missing_genes}",fontsize=9,wrap=True); axes[1].text(.03,.40,"Scores are strongly zero-inflated and retained only as descriptive context. Signature differences are not independent replication evidence.",fontsize=10,wrap=True,va="top")
    fig.suptitle("Supplementary 6. Signature–displacement associations and HSC-I sparsity")
    save_figure(fig,"dataset_b_supp6_signature_associations")


def evidence_ledger(data: dict[str, Any]) -> pd.DataFrame:
    comparator=data["comparator"].set_index(["contrast","state"]); tnf=data["tnf_primary"].set_index("state"); lps=data["lps_primary"].set_index("state")
    geometry=data["geometry"].set_index("rep"); exclude=data["exclude_diffusion"]
    correlations=data["correlations"]
    corr_tnf=correlations[(correlations.contrast=="PBS_vs_TNF")&(correlations.comparison_metric=="scgeo_normalized_displacement")].iloc[0]
    corr_lps=correlations[(correlations.contrast=="PBS_vs_LPS")&(correlations.comparison_metric=="scgeo_normalized_displacement")].iloc[0]
    hsc_lps=comparator.loc[("PBS_vs_LPS","HSC/quiescent")]
    rows=[
        ["TNF-associated geometry is stable for Activated HSC and Lymphoid states.",f"Activated HSC normalized displacement={tnf.loc['Activated HSC','normalized_magnitude_median']:.6f}; Lymphoid={tnf.loc['Lymphoid','normalized_magnitude_median']:.6f}; both primary_consensus_status=stable_effect; primary rank agreement ρ={tnf.pairwise_spearman_median.iloc[0]:.6f}.","PCA30, PCA50 and scVI primary consensus","ScGeo only",INFERENCE_STATUS,"Marker-inferred states; no independent biological replicates.","Dataset B Results / main panel B","Whether key TNF findings are representation-stable"],
        ["LPS shows broader state displacement than TNF in this descriptive dataset.",f"LPS: 7/7 states stable_effect, normalized displacement range {lps.normalized_magnitude_median.min():.6f}–{lps.normalized_magnitude_median.max():.6f}; TNF: 2/7 states stable_effect.","PCA30, PCA50 and scVI primary consensus","ScGeo only",INFERENCE_STATUS,"Breadth is descriptive and cannot establish a population-level treatment effect.","Dataset B Results / main panel B","Whether LPS response extends beyond the TNF-shifted states"],
        ["LPS HSC/quiescent is primarily notable for abundance change relative to its displacement rank.",f"Abundance proportion change={hsc_lps.abundance_change:+.6f}; normalized displacement={hsc_lps.scgeo_normalized_displacement:.6f} (rank {int(hsc_lps.scgeo_rank)}/7); energy distance={hsc_lps.energy_distance:.6f}.","Primary displacement; separate abundance and distribution axes","ScGeo geometry plus descriptive abundance/distribution",INFERENCE_STATUS,"Cells are not replicates; abundance is from one pooled condition library.","Dataset B Results / main panel D","Separation of composition, displacement and distribution"],
        ["Primary representation rankings are highly concordant.",f"Median pairwise state-rank Spearman: TNF={tnf.pairwise_spearman_median.iloc[0]:.6f}; LPS={lps.pairwise_spearman_median.iloc[0]:.6f}.","PCA30, PCA50 and scVI","ScGeo only",INFERENCE_STATUS,"Three representations are sensitivity views of the same cells, not independent confirmations.","Dataset B Results / main panels B–C","Robustness to representation choice"],
        ["Diffusion-map geometry is more distorted, with limited influence on primary scalar conclusions.",f"Diffusion median global distortion={geometry.loc['X_diffmap','median_global_distortion']:.6f}; state_graph_outlier={bool(geometry.loc['X_diffmap','state_graph_outlier'])}; maximum absolute change after excluding diffusion={exclude.absolute_change_after_excluding_diffusion.max():.6f}.","Diffusion exploratory sensitivity versus primary representations","ScGeo sensitivity only",INFERENCE_STATUS,"Scalar agreement does not validate diffusion coordinates; diffusion remains exploratory.","Supplementary representation-sensitivity figure","Concern that diffusion-map instability drives conclusions"],
        ["Official R Augur provides complementary treatment-separability evidence rather than the same estimand as ScGeo.",f"Official Augur–ScGeo Spearman: TNF={corr_tnf.spearman_rank_correlation:.6f}; LPS={corr_lps.spearman_rank_correlation:.6f}; n=7 marker-inferred states per contrast.","ScGeo primary consensus","Official R Augur 1.0.3 only in main figure",INFERENCE_STATUS,"Seven states only; rank association is not equivalence and has no biological-replicate inference.","Dataset B Results / main panel E","Comparator validity and complementarity"],
        ["Dataset B cannot support biological-replicate inference.","valid_biological_replicate_unit=false; three RNA matrices are pooled condition-level libraries; donor, cord-blood pool, recipient and mouse identities are unavailable per cell.","All representations","All comparators",INFERENCE_STATUS,"Cells, conditions, libraries, samples, condition_batch and SouporCell clades are not biological replicates.","Main limitation box and Methods","Pseudoreplication and inferential overclaiming"],
        ["Conservative annotations retain substantial uncertainty.",f"Annotation counts: ambiguous={data['annotation_summary']['confidence_counts']['ambiguous']}; low={data['annotation_summary']['confidence_counts']['low']}; moderate={data['annotation_summary']['confidence_counts']['moderate']}.","Marker-inferred labels","All state-level methods",INFERENCE_STATUS,"Rare or ambiguous states must not be over-annotated; labels are not reference-mapped ground truth.","Supplementary QC/annotation figure","Annotation confidence and state-definition sensitivity"],
        ["Computational controls assess cell-level stability only.",f"Balanced-control rows={len(data['balanced'])}; shuffled-label rows={len(data['shuffled'])}; shuffled control explicitly reports no p-value.","Primary representations and fixed computational resamples","ScGeo computational controls",INFERENCE_STATUS,"No computational null or cell resampling substitutes for biological replication.","Supplementary controls figure","Interpretation of resampling and null controls"],
    ]
    columns=["proposed_claim","supporting_output_and_exact_value","representation_scope","comparator_scope","descriptive_or_inferential_status","limitation","appropriate_manuscript_location","reviewer_concern_addressed"]
    return pd.DataFrame(rows,columns=columns)


def write_ledgers_and_text(data: dict[str, Any], hashes: dict[str, str]) -> None:
    ledger=evidence_ledger(data); atomic_csv(ledger,OUT/"dataset_b_evidence_ledger.csv")
    limitations={
        "dataset":"GSE249479 Dataset B","inference_status":INFERENCE_STATUS,"valid_biological_replicate_unit":False,
        "allowed_claim_scope":"descriptive state geometry, abundance, distribution and treatment separability within the analyzed cells",
        "prohibited_interpretations":["causal treatment effects","population-level inference","biological-replicate p-values","independent confirmation across PCA dimensions","UMAP geometry"],
        "limitations":[
            "Three pooled condition-level RNA libraries lack cell-level donor, cord-blood pool, recipient or mouse identity.",
            "Condition, sample, condition_batch, library and SouporCell clade are not biological replicate identifiers.",
            "Annotations are marker-inferred: 10,154 cells are ambiguous and 10,272 have low confidence.",
            "Official Augur and cell-resampling uncertainty measure computational stability, not biological uncertainty.",
            "Rank correlations contain only seven marker-inferred states per contrast.",
            "Diffusion-map geometry is an exploratory sensitivity with elevated distortion.",
            "HSC-I contains 45/50 requested genes and strongly zero-inflated scores.",
            "UMAP is display-only and was not used for quantitative geometry."
        ],
        "official_comparator":{"name":"Official R Augur","version":"1.0.3","commit":"b252b84e4af687d9817813b1db409267eb44ec3f"},
        "python_comparator_role":"supplementary implementation sensitivity only; never averaged or merged with official AUC",
        "frozen_source_hashes":hashes,
    }
    atomic_json(limitations,OUT/"dataset_b_claims_and_limitations.json")
    main_caption="""**Dataset B (GSE249479), descriptive treatment geometry.** (A) PCA50-derived UMAPs display conservative marker-inferred states and condition distributions; UMAP is visualization only and is not used for quantitative geometry. (B) Primary-consensus normalized displacement for PBS→TNF and PBS→LPS, with filled symbols indicating `stable_effect` across PCA30, PCA50 and scVI. (C) Representation sensitivity across the three primary representations, PCA20 dimensional sensitivity and exploratory diffusion map. (D) Abundance change, robust displacement and within-state energy distance are shown as separate columns with independent scales; no composite score is formed. (E) Official R Augur 1.0.3 treatment separability versus ScGeo representation-stable displacement; the quantities are complementary, not equivalent. (F) Dataset B is `descriptive_only` because no independent biological replicate identity is available. Cells, conditions, libraries and SouporCell clades are not replicates."""
    atomic_text(main_caption,OUT/"dataset_b_main_figure_caption.md")
    supplementary="""**Supplementary 1 — QC and annotation confidence.** Sparse-QC summaries and conservative marker-label confidence; annotation remains marker-inferred and includes ambiguous and low-confidence cells.

**Supplementary 2 — Representation quality.** Marker-label and mapped-clade neighbourhood preservation, PCA50 neighbourhood overlap, local-distance agreement and graph connectivity. Condition mixing is descriptive because condition is biological.

**Supplementary 3 — Diffusion and leave-one-representation sensitivity.** Diffusion map has elevated distortion and state-graph outlier status, while excluding it changes median normalized effects by at most 0.016432. Primary conclusions use PCA30, PCA50 and scVI.

**Supplementary 4 — Computational controls.** Repeated cell-balanced subsampling and shuffled-cell-label controls quantify computational stability only; they provide no biological uncertainty or inferential p-values.

**Supplementary 5 — Augur implementation sensitivity.** Official R Augur is the primary comparator. The Augur-inspired Python approximation is a supplementary implementation sensitivity using raw sparse variance rather than official LOESS-residual variable-feature selection; AUCs are not averaged or merged.

**Supplementary 6 — Signature associations.** Descriptive across-state signature–displacement rank associations (seven marker-inferred states per contrast). The study HSC-I score retains only 45/50 requested genes and is strongly zero-inflated."""
    atomic_text(supplementary,OUT/"dataset_b_supplementary_captions.md")

    methods=pd.DataFrame([
        ["Sparse QC and annotation","Phase 2 source","288fa73","compact H5AD and audit tables","marker-inferred; descriptive_only"],
        ["Representation ensemble","PCA20/30/50, diffusion, scVI","a669852","frozen representation object and Phase 3A tables","UMAP display-only"],
        ["Treatment geometry","Frozen ScGeo defaults","bca6ace",f"ScGeo package {SCGEO_COMMIT}","PCA30/PCA50/scVI primary; no replicate inference"],
        ["Primary comparator","Official R Augur 1.0.3","e307a1b","commit b252b84e4af687d9817813b1db409267eb44ec3f","treatment separability; cell-level CV only"],
        ["Implementation sensitivity","Augur-inspired Python approximation","76bfd4a","sklearn RF; raw sparse-variance features","supplementary only; no retrospective tuning"],
        ["Manuscript assembly","Frozen-artifact plotting and evidence collation","working tree source","all source hashes validated before and after","no numerical analysis recomputed"],
    ],columns=["stage","implementation","repository_commit","numerical_source","scope_and_limitation"])
    methods["inference_status"]=INFERENCE_STATUS; atomic_csv(methods,OUT/"dataset_b_methods_provenance.csv")
    reviewer=pd.DataFrame([
        ["Biological replication","Replication reassessment and limitation box","Explicitly prevents pseudoreplicate inference","No biological replicate identity can be recovered"],
        ["Representation robustness","Primary rank agreement, five-representation matrix, leave-one analyses","Separates primary consensus from sensitivity views","Representations share the same cells"],
        ["Comparator validity","Official R Augur primary; Python supplementary","Uses the official package and preserves implementation differences","Only seven marker-inferred states"],
        ["Annotation uncertainty","Confidence counts and ambiguous-state visibility","Avoids over-annotation","Labels remain marker-inferred"],
        ["Abundance versus geometry","Separate abundance, displacement and distribution columns","Avoids opaque composite score","All axes remain descriptive"],
        ["Computational controls","Balanced and shuffled-cell-label panels","Shows cell-level stability and negative controls","Controls do not create biological uncertainty"],
        ["Reproducibility","Exact frozen hashes, versions and commits","Auditable manuscript assembly without recomputation","Depends on preserved ignored result artifacts"],
    ],columns=["reviewer_concern","relevant_output","how_addressed","residual_limitation"])
    reviewer["inference_status"]=INFERENCE_STATUS; atomic_csv(reviewer,OUT/"dataset_b_reviewer_relevance.csv")

    alt={
        "dataset_b_main_figure":"Six-panel Dataset B figure. Panel A shows display-only PCA50 UMAPs colored by seven conservative marker-inferred states and by PBS, TNF and LPS. Panel B shows normalized state displacement for PBS-to-TNF and PBS-to-LPS, with stable effects filled; TNF stable effects occur for Activated HSC and Lymphoid, whereas all seven LPS states are stable effects. Panel C is a state-by-representation heatmap separating PCA30, PCA50 and scVI primary views from PCA20 and diffusion sensitivities. Panel D shows abundance change, displacement and energy distance in separate columns with no composite. Panel E plots official R Augur AUC against ScGeo displacement for seven states per contrast. Panel F states descriptive_only and no independent biological replicate identity.",
        "dataset_b_supp1_qc_annotation":"QC percentile summaries, marker-label confidence by state, and overall ambiguous, low and moderate annotation counts.",
        "dataset_b_supp2_representation_quality":"Heatmap of five representation-quality metrics across PCA20, PCA30, PCA50, diffusion map and scVI; diffusion has weaker local preservation than the primary representations.",
        "dataset_b_supp3_diffusion_leave_one":"Bars and points show diffusion-map distortion, small changes after excluding diffusion, and leave-one-representation sensitivity by state and contrast.",
        "dataset_b_supp4_computational_controls":"Balanced cell-subsampling and shuffled-cell-label normalized displacement ranges, explicitly labeled computational and non-biological.",
        "dataset_b_supp5_official_python_augur":"Official R Augur AUC versus the Augur-inspired Python approximation and feature-set Jaccard overlap; official output is primary and Python is supplementary.",
        "dataset_b_supp6_signature_associations":"Heatmap of signature-displacement Spearman correlations across seven states and a text panel noting 45 of 50 HSC-I genes available with zero-inflation limitation.",
    }
    for name,text in alt.items(): atomic_text(text,OUT/"alt_text"/f"{name}.txt")


def assemble() -> dict[str, Any]:
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titlesize":11,"axes.labelsize":9,"svg.hashsalt":"gse249479-dataset-b-manuscript-v1"})
    before=validate_frozen_sources(); data=read_sources()
    main_figure(data); supplementary_figures(data); write_ledgers_and_text(data,before)
    after=validate_frozen_sources()
    if before != after: raise RuntimeError("A frozen numerical source changed during manuscript assembly")
    generated=sorted(str(path.relative_to(ROOT)) for path in OUT.rglob("*") if path.is_file())
    report={"status":"passed","inference_status":INFERENCE_STATUS,"numerical_analysis_recomputed":False,"scgeo_rerun":False,"comparator_rerun":False,"normalization_run":False,"feature_selection_run":False,"source_hashes_before":before,"source_hashes_after":after,"source_hashes_unchanged":before==after,"frozen_scgeo_commit":SCGEO_COMMIT,"frozen_scgeo_clean":True,"generated_files":generated}
    atomic_json(report,OUT/"dataset_b_manuscript_assembly_validation.json")
    return report


if __name__ == "__main__":
    print(json.dumps(assemble(),indent=2,sort_keys=True))
