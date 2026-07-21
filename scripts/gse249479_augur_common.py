"""Official Augur comparator plus a descriptive-only Python sensitivity analysis."""

from __future__ import annotations

import json
import os
import platform
import time
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from scripts.gse249479_representation_common import (
    ROOT, ResourceLog, active_branch, atomic_csv, atomic_json, close_and_collect,
    ensure_tree, git_commit, package_versions, sha256, utc_now,
)


CONFIG_PATH = ROOT / "configs/gse249479_augur_v1.json"


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def paths(config: dict[str, Any] | None = None) -> dict[str, Path]:
    config = load_config() if config is None else config
    out = Path(os.environ.get("SCGEO_GSE249479_OUTPUT_DIR", config["default_output_dir"]))
    if not out.is_absolute(): out = ROOT / out
    return {
        "compact": Path(os.environ.get("SCGEO_GSE249479_COMPACT_H5AD", config["default_compact_h5ad"])).resolve(),
        "representation": Path(os.environ.get("SCGEO_GSE249479_REPRESENTATION_H5AD", config["default_representation_h5ad"])).resolve(),
        "output": out.resolve(),
        "scgeo": Path(os.environ.get("SCGEO_SOURCE_REPO", ROOT.parent / "scgeo")).resolve(),
    }


def ensure_dirs(out: Path) -> None:
    ensure_tree(out)
    for name in ["comparator", "controls"]: (out / name).mkdir(parents=True, exist_ok=True)


def validate(config: dict[str, Any], p: dict[str, Path]) -> None:
    if active_branch() != config["required_branch"]: raise RuntimeError("Wrong branch")
    if git_commit(p["scgeo"]) != config["frozen_scgeo_commit"]: raise RuntimeError("Frozen ScGeo commit mismatch")
    if sha256(p["compact"]) != config["expected_compact_sha256"]: raise RuntimeError("Compact checksum mismatch")
    if sha256(p["representation"]) != config["expected_representation_sha256"]: raise RuntimeError("Representation checksum mismatch")


def write_resources(resources: ResourceLog, out: Path) -> None:
    frame = pd.DataFrame(resources.rows)
    frame["inference_status"] = "descriptive_only"
    frame["uncertainty_status"] = "computational_cell_resampling_stability_only_where_applicable"
    atomic_csv(frame, out / "audit" / f"{resources.stage}_resource_log.csv")
    combined = out / "audit" / "05_06_comparator_resource_log.csv"
    if combined.exists(): frame = pd.concat([pd.read_csv(combined), frame], ignore_index=True)
    atomic_csv(frame, combined)


def state_labels(obs: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    obs = obs.copy()
    obs["state_detailed"] = pd.Categorical(obs["marker_inferred_label"].astype(str).map(config["detailed_state_mapping"]), categories=config["detailed_state_order"], ordered=True)
    obs["state_coarse"] = pd.Categorical(obs["state_detailed"].astype(str).map(config["coarse_state_mapping"]), categories=config["coarse_state_order"], ordered=True)
    return obs


def load_sparse_hvg(config: dict[str, Any], p: dict[str, Path]) -> tuple[sparse.csr_matrix, pd.DataFrame, np.ndarray]:
    backed = ad.read_h5ad(p["compact"], backed="r")
    hvg = backed.var["highly_variable"].to_numpy(bool)
    if int(hvg.sum()) != int(config["expression"]["n_hvg"]): raise RuntimeError("Expected exactly 3,000 HVGs")
    matrix = backed[:, hvg].X
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    X = matrix.tocsr().astype(np.float32)
    obs = state_labels(backed.obs.copy(), config)
    genes = backed.var_names[hvg].astype(str).to_numpy()
    backed.file.close()
    totals = obs["total_counts"].to_numpy(np.float32)
    factors = np.divide(np.float32(1e4), totals, out=np.zeros_like(totals), where=totals > 0)
    X = X.multiply(factors[:, None]).tocsr().astype(np.float32)
    np.log1p(X.data, out=X.data)
    if not sparse.isspmatrix_csr(X): raise RuntimeError("Augur input lost CSR encoding")
    return X, obs, genes


def variable_features(X: sparse.csr_matrix, idx: np.ndarray, retain_fraction: float) -> np.ndarray:
    block = X[idx]
    mean = np.asarray(block.mean(axis=0)).ravel()
    sqmean = np.asarray(block.power(2).mean(axis=0)).ravel()
    variance = np.maximum(sqmean - mean * mean, 0)
    n = max(2, int(np.ceil(variance.size * retain_fraction)))
    return np.argsort(variance, kind="mergesort")[-n:]


def one_subsample_auc(X: sparse.csr_matrix, idx0: np.ndarray, idx1: np.ndarray, features: np.ndarray, seed: int, cfg: dict[str, Any]) -> tuple[float, list[float], int]:
    rng = np.random.RandomState(seed)
    n = int(cfg["subsample_size_per_condition"])
    chosen = np.concatenate([rng.choice(idx0, n, replace=False), rng.choice(idx1, n, replace=False)])
    y = np.concatenate([np.zeros(n, dtype=np.int8), np.ones(n, dtype=np.int8)])
    nfeat = max(2, int(np.ceil(features.size * float(cfg["random_feature_fraction"]))))
    picked = np.sort(rng.choice(features, nfeat, replace=False))
    block = X[chosen][:, picked]
    mean = np.asarray(block.mean(axis=0)).ravel(); sq = np.asarray(block.power(2).mean(axis=0)).ravel()
    keep = (sq - mean * mean) > 0
    block = block[:, keep]
    if block.shape[1] == 0: return np.nan, [], 0
    cv = StratifiedKFold(n_splits=int(cfg["cv_folds"]), shuffle=True, random_state=seed)
    fold_auc = []
    for fold, (train, test) in enumerate(cv.split(np.zeros(y.size), y)):
        model = RandomForestClassifier(
            n_estimators=int(cfg["n_trees"]), max_features=int(cfg["max_features_mtry"]),
            min_samples_split=int(cfg["min_samples_split"]), random_state=seed + fold,
            n_jobs=int(cfg["n_jobs"]), class_weight=None,
        )
        model.fit(block[train], y[train])
        fold_auc.append(float(roc_auc_score(y[test], model.predict_proba(block[test])[:, 1])))
    return float(np.mean(fold_auc)), fold_auc, int(block.shape[1])


def run_hierarchy(X: sparse.csr_matrix, obs: pd.DataFrame, config: dict[str, Any], hierarchy: str, shuffled: bool = False) -> pd.DataFrame:
    key = f"state_{hierarchy}"; states = config[f"{hierarchy}_state_order"]
    cond = obs["condition"].astype(str).to_numpy(); labels = obs[key].astype(str).to_numpy()
    aug = config["augur_defaults"]
    blocks = config["shuffled_null"]["fixed_seed_blocks"] if shuffled else config["fixed_seed_blocks"]
    n_inner = int(config["shuffled_null"]["subsamples_per_seed_block"] if shuffled else config["subsamples_per_seed_block"])
    rows = []
    for contrast_i, treated in enumerate(config["contrasts"]):
        for state_i, state in enumerate(states):
            pool = np.flatnonzero(np.isin(cond, ["PBS", treated]) & (labels == state))
            ypool = cond[pool].copy()
            n0_orig, n1_orig = int(np.sum(ypool == "PBS")), int(np.sum(ypool == treated))
            if min(n0_orig, n1_orig) < int(aug["subsample_size_per_condition"]):
                rows.append({"contrast": f"PBS_vs_{treated}", "hierarchy": hierarchy, "state": state, "status": "insufficient_coverage", "n_cells_PBS": n0_orig, "n_cells_treated": n1_orig, "inference_status": "descriptive_only"})
                continue
            features = variable_features(X, pool, float(aug["variance_quantile_retained"]))
            for block_i, block_seed in enumerate(blocks):
                if shuffled:
                    rng = np.random.RandomState(int(block_seed) + 100 * contrast_i + state_i)
                    assigned = ypool.copy(); rng.shuffle(assigned)
                    idx0, idx1 = pool[assigned == "PBS"], pool[assigned == treated]
                else:
                    idx0, idx1 = pool[ypool == "PBS"], pool[ypool == treated]
                for inner in range(n_inner):
                    seed = int(block_seed) + 10000 * contrast_i + 100 * state_i + inner
                    auc, folds, nfeat = one_subsample_auc(X, idx0, idx1, features, seed, aug)
                    rows.append({
                        "contrast": f"PBS_vs_{treated}", "hierarchy": hierarchy, "state": state,
                        "seed_block": int(block_seed), "subsample_in_block": inner, "subsample_seed": seed,
                        "auc": auc, "fold_auc_1": folds[0] if folds else np.nan, "fold_auc_2": folds[1] if folds else np.nan,
                        "fold_auc_3": folds[2] if folds else np.nan, "n_cells_PBS_available": n0_orig,
                        "n_cells_treated_available": n1_orig, "n_cells_used_per_condition": int(aug["subsample_size_per_condition"]),
                        "n_variable_features": int(features.size), "n_features_used": nfeat, "status": "usable",
                        "control_status": config["shuffled_null"]["status"] if shuffled else "observed_labels",
                        "uncertainty_status": "computational_cell_resampling_stability_only", "inference_status": "descriptive_only",
                    })
    return pd.DataFrame(rows)


def summarize_auc(rows: pd.DataFrame) -> pd.DataFrame:
    usable = rows[rows["status"] == "usable"]
    summary = usable.groupby(["contrast", "hierarchy", "state"], observed=False).agg(
        augur_auc=("auc", "mean"), auc_sd=("auc", "std"), auc_median=("auc", "median"),
        auc_q025=("auc", lambda x: x.quantile(.025)), auc_q975=("auc", lambda x: x.quantile(.975)),
        n_computational_subsamples=("auc", "size"), n_cells_PBS_available=("n_cells_PBS_available", "first"),
        n_cells_treated_available=("n_cells_treated_available", "first"), n_cells_used_per_condition=("n_cells_used_per_condition", "first"),
    ).reset_index()
    summary["status"] = "usable"; summary["metric_interpretation"] = "treatment_separability_not_effect_magnitude_or_direction"
    summary["uncertainty_status"] = "computational_cell_resampling_stability_only"; summary["inference_status"] = "descriptive_only"
    return summary


def seed_summary(rows: pd.DataFrame) -> pd.DataFrame:
    out = rows[rows["status"] == "usable"].groupby(["contrast", "hierarchy", "state", "seed_block"], observed=False).agg(seed_block_auc=("auc", "mean"), seed_block_auc_sd=("auc", "std"), n_subsamples=("auc", "size")).reset_index()
    out["uncertainty_status"] = "computational_fixed_seed_block_sensitivity"; out["inference_status"] = "descriptive_only"
    return out


def metadata(stage: str, config: dict[str, Any], p: dict[str, Path], resources: ResourceLog, payload: dict[str, Any], implementation: dict[str, Any] | None = None) -> None:
    implementation = config["implementation"] if implementation is None else implementation
    record = {**payload, "stage": stage, "timestamp_utc": utc_now(), "inference_status": "descriptive_only",
              "implementation": implementation,
              "compact_sha256": sha256(p["compact"]), "representation_sha256": sha256(p["representation"]),
              "frozen_scgeo_commit": git_commit(p["scgeo"]), "scgeo_modified": False, "biological_replicates": False,
              "packages": package_versions(), "peak_cpu_rss_gib": resources.peak_rss_gib}
    if implementation.get("role") == "primary_external_comparator" or implementation.get("manuscript_role") == "primary_external_comparator":
        record["official_augur_settings"] = config["canonical_comparator"]["settings"]
    else:
        record["python_augur_inspired_headline_settings"] = config["augur_defaults"]
    atomic_json(record, p["output"] / "metadata" / f"{stage}_metadata.json")
    atomic_json({"timestamp_utc": utc_now(), "python": platform.python_version(), "packages": package_versions(), "implementation": implementation, "inference_status": "descriptive_only"}, p["output"] / "version_records" / f"{stage}_versions.json")


def python_sensitivity_export(frame: pd.DataFrame) -> pd.DataFrame:
    """Label Python results precisely without changing any numeric value."""
    rename = {
        "auc": "python_augur_inspired_auc",
        "augur_auc": "python_augur_inspired_auc",
        "seed_block_auc": "python_augur_inspired_seed_block_auc",
        "seed_block_auc_sd": "python_augur_inspired_seed_block_auc_sd",
        "detailed_augur_auc": "detailed_python_augur_inspired_auc",
        "coarse_augur_auc": "coarse_python_augur_inspired_auc",
        "observed_augur_auc": "observed_python_augur_inspired_auc",
    }
    out = frame.rename(columns={key: value for key, value in rename.items() if key in frame.columns}).copy()
    out["implementation"] = "Augur-inspired Python approximation"
    out["manuscript_role"] = "supplementary_implementation_sensitivity"
    return out


def run_augur() -> dict[str, Any]:
    config, p = load_config(), paths(); validate(config, p); ensure_dirs(p["output"])
    stage = "05_augur_comparator"; resources = ResourceLog(stage, config); started = time.perf_counter()
    before = {k: sha256(p[k]) for k in ["compact", "representation"]}
    with resources.operation("load_sparse_normalized_3000_hvgs"):
        X, obs, genes = load_sparse_hvg(config, p)
    with resources.operation("augur_observed_detailed_and_coarse"):
        observed = pd.concat([run_hierarchy(X, obs, config, h, False) for h in ["detailed", "coarse"]], ignore_index=True, sort=False)
        summary = summarize_auc(observed); seeds = seed_summary(observed)
        seeds["seed_block_rank"] = seeds.groupby(["contrast", "hierarchy", "seed_block"], observed=False)["seed_block_auc"].rank(method="average", ascending=False)
        rank_stability = seeds.groupby(["contrast", "hierarchy", "state"], observed=False).agg(
            seed_rank_sd=("seed_block_rank", "std"), seed_rank_min=("seed_block_rank", "min"), seed_rank_max=("seed_block_rank", "max")
        ).reset_index()
        summary = summary.merge(rank_stability, on=["contrast", "hierarchy", "state"], how="left")
        summary["ranking_stability_status"] = np.where(summary["seed_rank_sd"] > 1.0, "ranking_sensitive_across_fixed_seed_blocks", "ranking_stable_across_fixed_seed_blocks")
    with resources.operation("label_shuffled_computational_null"):
        shuffled = pd.concat([run_hierarchy(X, obs, config, h, True) for h in ["detailed", "coarse"]], ignore_index=True, sort=False)
        shuffled_summary = summarize_auc(shuffled)
        shuffled_summary["control_status"] = config["shuffled_null"]["status"]
    with resources.operation("write_augur_outputs"):
        out = p["output"]
        atomic_csv(python_sensitivity_export(observed), out / "comparator" / "05_augur_subsample_auc.csv.gz")
        atomic_csv(python_sensitivity_export(summary), out / "comparator" / "05_augur_state_summary.csv")
        atomic_csv(python_sensitivity_export(summary), out / "comparator" / "05_python_augur_inspired_state_summary.csv")
        atomic_csv(python_sensitivity_export(summary), out / "figure_sources" / "05_augur_state_summary.csv")
        atomic_csv(python_sensitivity_export(seeds), out / "controls" / "05_augur_fixed_seed_sensitivity.csv")
        atomic_csv(python_sensitivity_export(shuffled), out / "controls" / "05_augur_shuffled_subsample_auc.csv.gz")
        atomic_csv(python_sensitivity_export(shuffled_summary), out / "controls" / "05_augur_shuffled_null_summary.csv")
        null_figure_source = shuffled_summary.merge(
            summary[["contrast", "hierarchy", "state", "augur_auc"]].rename(columns={"augur_auc": "observed_augur_auc"}),
            on=["contrast", "hierarchy", "state"], how="left",
        )
        atomic_csv(python_sensitivity_export(null_figure_source), out / "figure_sources" / "05_augur_sensitivity_and_null.csv")
        exclusion = summary[(summary["hierarchy"] == "detailed") & (summary["state"] != "Ambiguous HSPC")].copy()
        exclusion["sensitivity_status"] = "ambiguous_state_excluded_other_statewise_scores_unchanged"
        atomic_csv(python_sensitivity_export(exclusion), out / "controls" / "05_augur_exclude_ambiguous_hspc.csv")
        detailed_h = summary[summary["hierarchy"] == "detailed"][["contrast", "state", "augur_auc", "auc_sd", "seed_rank_sd"]].rename(columns={"state": "detailed_state", "augur_auc": "detailed_augur_auc", "auc_sd": "detailed_auc_sd", "seed_rank_sd": "detailed_seed_rank_sd"})
        detailed_h["coarse_state"] = detailed_h["detailed_state"].map(config["coarse_state_mapping"])
        coarse_h = summary[summary["hierarchy"] == "coarse"][["contrast", "state", "augur_auc", "auc_sd", "seed_rank_sd"]].rename(columns={"state": "coarse_state", "augur_auc": "coarse_augur_auc", "auc_sd": "coarse_auc_sd", "seed_rank_sd": "coarse_seed_rank_sd"})
        hierarchy = detailed_h.merge(coarse_h, on=["contrast", "coarse_state"], how="left")
        hierarchy["coarse_minus_detailed_auc"] = hierarchy["coarse_augur_auc"] - hierarchy["detailed_augur_auc"]
        hierarchy["sensitivity_status"] = "detailed_vs_coarse_state_definition_no_combined_score"; hierarchy["inference_status"] = "descriptive_only"
        atomic_csv(python_sensitivity_export(hierarchy), out / "controls" / "05_augur_hierarchy_sensitivity.csv")
    after = {k: sha256(p[k]) for k in ["compact", "representation"]}
    if before != after: raise RuntimeError("Immutable input changed")
    payload = {"status": "passed", "implementation": "Augur-inspired Python approximation", "manuscript_role": "supplementary_implementation_sensitivity", "n_summary_rows": int(summary.shape[0]), "n_observed_subsamples": int((observed.status == 'usable').sum()),
               "n_shuffled_subsamples": int((shuffled.status == 'usable').sum()), "runtime_seconds": time.perf_counter()-started,
               "input_sha256_before": before, "input_sha256_after": after, "full_matrix_densified": False}
    atomic_json(payload, p["output"] / "audit" / "05_augur_summary.json"); write_resources(resources, p["output"]); metadata(stage, config, p, resources, payload)
    del X; close_and_collect(); return payload


def _configured_result_path(path: str) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else ROOT / candidate


def official_augur_summary(config: dict[str, Any], p: dict[str, Path]) -> pd.DataFrame:
    canonical = config["canonical_comparator"]
    auc = pd.read_csv(_configured_result_path(canonical["state_auc_path"]))
    subs = pd.read_csv(_configured_result_path(canonical["subsample_auc_path"]))
    per_subsample = subs.groupby(["contrast", "state", "subsample_idx"], observed=False)["estimate"].mean().reset_index()
    spread = per_subsample.groupby(["contrast", "state"], observed=False)["estimate"].agg(
        official_auc_sd="std", official_auc_median="median", n_computational_subsamples="size"
    ).reset_index()
    auc = auc.merge(spread, on=["contrast", "state"], validate="one_to_one")
    auc = auc.rename(columns={"official_augur_auc": "official_augur_auc"})

    python_summary = pd.read_csv(p["output"] / "comparator" / "05_augur_state_summary.csv")
    python_summary = python_summary[python_summary["hierarchy"] == "detailed"].copy()
    counts = python_summary[["contrast", "state", "n_cells_PBS_available", "n_cells_treated_available", "n_cells_used_per_condition", "status"]]
    auc = auc.merge(counts, on=["contrast", "state"], validate="one_to_one")
    auc["implementation"] = "Official R Augur"
    auc["package_version"] = "Augur 1.0.3"
    auc["package_commit"] = canonical["commit"]
    auc["manuscript_role"] = "primary_external_comparator"
    auc["metric_interpretation"] = "treatment_separability_not_effect_magnitude_or_direction"
    auc["uncertainty_status"] = "computational_cell_resampling_stability_only"
    auc["inference_status"] = "descriptive_only"
    auc["rank_scope"] = "seven_marker_inferred_states_per_contrast"
    return auc


def comparator_methods_table() -> pd.DataFrame:
    rows = [
        {
            "implementation": "Official R Augur",
            "package/engine": "Augur 1.0.3; R randomForest 4.7-1.2",
            "subsampling": "50 balanced subsamples; 20 cells per condition",
            "CV": "3-fold stratified rsample::vfold_cv",
            "variable-feature selection": "official LOESS residuals of mean/SD after 1% tail trimming; upper 50% residuals",
            "random feature selection": "random 50% of selected features per subsample using floor(n_features*0.5)",
            "forest parameters": "100 trees; mtry=2; R randomForest engine",
            "role in manuscript": "primary external comparator",
        },
        {
            "implementation": "Augur-inspired Python approximation",
            "package/engine": "scikit-learn RandomForestClassifier",
            "subsampling": "50 balanced subsamples; 20 cells per condition; different fixed-seed schedule",
            "CV": "3-fold sklearn StratifiedKFold; different split construction and seeds",
            "variable-feature selection": "upper half of genes by raw sparse variance; not official LOESS residuals",
            "random feature selection": "random 50% of selected features using ceil(n_features*0.5)",
            "forest parameters": "100 trees; max_features/mtry=2; sklearn engine",
            "role in manuscript": "supplementary implementation sensitivity",
        },
    ]
    table = pd.DataFrame(rows)
    table["inference_status"] = "descriptive_only"
    return table


def implementation_sensitivity(config: dict[str, Any]) -> pd.DataFrame:
    table = pd.read_csv(_configured_result_path(config["canonical_comparator"]["implementation_comparison_path"]))
    required = {"official_augur_auc", "python_augur_inspired_auc", "absolute_auc_difference", "jaccard_overlap"}
    if missing := required.difference(table.columns):
        raise RuntimeError(f"Official/Python sensitivity table lacks corrected columns: {sorted(missing)}")
    table["official_role"] = "primary_external_comparator"
    table["python_role"] = "supplementary_implementation_sensitivity"
    table["rank_scope"] = "seven_marker_inferred_states_per_contrast"
    table["no_auc_combination"] = True
    table["inference_status"] = "descriptive_only"
    return table


def build_comparator(config: dict[str, Any], p: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    aug = official_augur_summary(config, p)
    sc = pd.read_csv(p["output"] / "scgeo" / "04_full_state_evidence.csv")
    sc = sc[(sc["hierarchy"] == "detailed") & (sc["representation_role"] == "primary")].copy()
    metrics = sc.groupby(["contrast", "state"], observed=False).agg(
        scgeo_normalized_displacement=("normalized_effect", "median"),
        representation_consensus_status=("consensus_status", "first"),
        abundance_change=("abundance_proportion_change_treated_minus_PBS", "first"),
        energy_distance=("energy_distance_biased_full_value", "median"), mmd=("mmd_rbf_squared_biased_full_value", "median"),
        sliced_wasserstein=("sliced_wasserstein_full_value", "median"),
        leave_one_representation_max_relative_change=("leave_one_representation_max_relative_change", "first"),
    ).reset_index()
    signatures = pd.read_csv(p["output"] / "scgeo" / "04_state_signature_changes.csv")
    signatures = signatures[signatures["signature"].isin(config["principal_signatures"])].pivot(index=["contrast", "state"], columns="signature", values="mean_difference_treated_minus_PBS").add_prefix("signature_change__").reset_index()
    rep = ad.read_h5ad(p["representation"], backed="r")
    label_map = rep.obs["marker_inferred_label"].astype(str).map(config["detailed_state_mapping"])
    confidence = rep.obs.assign(state=label_map).groupby("state", observed=False)["annotation_confidence"].value_counts(normalize=True).unstack(fill_value=0).add_prefix("annotation_fraction_").reset_index()
    rep.file.close()
    table = aug.merge(metrics, on=["contrast", "state"], how="left").merge(signatures, on=["contrast", "state"], how="left").merge(confidence, on="state", how="left")
    coarse_parts = []
    for slug in ["tnf", "lps"]:
        part = pd.read_csv(p["output"] / "scgeo" / f"04_{slug}_coarse_consensus.csv")
        coarse_parts.append(part[["contrast", "state", "all_representation_consensus_status"]].rename(columns={"state": "coarse_state", "all_representation_consensus_status": "coarse_all_representation_consensus_status"}))
    coarse = pd.concat(coarse_parts, ignore_index=True)
    table["coarse_state"] = table["state"].map(config["coarse_state_mapping"])
    table = table.merge(coarse, on=["contrast", "coarse_state"], how="left")
    table["study_hsc_i_limitation"] = "45_of_50_genes_available_and_scores_strongly_zero_inflated_descriptive_context_only"
    table["official_augur_rank"] = table.groupby("contrast")["official_augur_auc"].rank(method="average", ascending=False)
    table["scgeo_rank"] = table.groupby("contrast")["scgeo_normalized_displacement"].rank(method="average", ascending=False)
    rules = config["comparison_pattern_rules"]
    def label(row: pd.Series) -> str:
        a, s = row.official_augur_rank, row.scgeo_rank
        if a <= rules["high_rank_max"] and s <= rules["high_rank_max"]: return "high_official_Augur_high_ScGeo"
        if a <= rules["high_rank_max"] and s >= rules["moderate_or_low_rank_min"]: return "high_official_Augur_moderate_or_low_ScGeo"
        if s <= rules["high_rank_max"] and a >= rules["moderate_or_low_rank_min"]: return "high_ScGeo_moderate_or_low_official_Augur"
        if a >= rules["low_rank_min"] and s >= rules["low_rank_min"]: return "low_both"
        return "intermediate_or_mixed"
    table["transparent_rank_pattern"] = table.apply(label, axis=1)
    table["method_relation"] = "complementary_non_equivalent_quantities_no_combined_score"
    table["inference_status"] = "descriptive_only"
    correlations = []
    for contrast, group in table.groupby("contrast", observed=False):
        for metric in ["scgeo_normalized_displacement", "abundance_change", "energy_distance", "mmd", "sliced_wasserstein"]:
            correlations.append({"contrast": contrast, "augur_metric": "official_R_Augur_auc_treatment_separability", "comparison_metric": metric,
                                 "spearman_rank_correlation": group["official_augur_auc"].corr(group[metric], method="spearman"), "n_states": int(group.shape[0]),
                                 "rank_scope": "seven_marker_inferred_states", "interpretation": "descriptive_rank_association_not_method_equivalence", "inference_status": "descriptive_only"})
    return table, pd.DataFrame(correlations)


def figures(out: Path, table: pd.DataFrame, correlations: pd.DataFrame, seed: pd.DataFrame, null: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 180, "font.size": 8})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, contrast in zip(axes, ["PBS_vs_TNF", "PBS_vs_LPS"]):
        x = table[table.contrast == contrast].sort_values("official_augur_auc")
        ax.barh(x.state, x.official_augur_auc, xerr=x.official_auc_sd, color="#4C78A8"); ax.axvline(.5, color="grey", ls="--", lw=.8)
        ax.set(title=contrast.replace("_", " "), xlabel="Official R Augur AUC (separability)", xlim=(0.45, 1.0))
    fig.suptitle("A. Official R Augur treatment-separability ranking (primary comparator)"); fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(out / "figures" / f"05_augur_state_ranking.{ext}")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    for contrast, g in table.groupby("contrast"):
        ax.scatter(g.scgeo_normalized_displacement, g.official_augur_auc, label=contrast)
        for r in g.itertuples(): ax.annotate(r.state, (r.scgeo_normalized_displacement, r.official_augur_auc), fontsize=6)
    ax.set(xlabel="ScGeo primary normalized displacement", ylabel="Official R Augur AUC", title="B. Official separability versus displacement geometry"); ax.legend(); fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(out / "figures" / f"06_augur_vs_scgeo_scatter.{ext}")
    plt.close(fig)
    cols=["official_augur_auc","scgeo_normalized_displacement","abundance_change","energy_distance","mmd","sliced_wasserstein"]
    mat=table.set_index(["contrast","state"])[cols].copy(); shown=mat.copy()
    for c in cols:
        v=shown[c].abs() if c=="abundance_change" else shown[c]; lo,hi=v.min(),v.max(); shown[c]=(v-lo)/(hi-lo) if hi>lo else 0
    fig,ax=plt.subplots(figsize=(8.2,6)); im=ax.imshow(shown,aspect="auto",cmap="viridis")
    ax.set_xticks(range(len(cols)),["Official R Augur AUC","ScGeo displacement","|abundance change|","energy","MMD","sliced Wasserstein"],rotation=30,ha="right")
    ax.set_yticks(range(shown.shape[0]),[f"{a}: {b}" for a,b in shown.index]); ax.set_title("C. Method-evidence matrix (columnwise display scaling; no composite)"); fig.colorbar(im,ax=ax,label="Within-column display scale"); fig.tight_layout()
    for ext in ["png","svg"]: fig.savefig(out/"figures"/f"06_method_evidence_matrix.{ext}")
    plt.close(fig)
    seed_col = "python_augur_inspired_seed_block_auc" if "python_augur_inspired_seed_block_auc" in seed else "seed_block_auc"
    null_col = "python_augur_inspired_auc" if "python_augur_inspired_auc" in null else "augur_auc"
    null_d=null[null.hierarchy=="detailed"][["contrast","state",null_col]].rename(columns={null_col:"null_auc"})
    python_auc = sensitivity[["contrast","state","python_augur_inspired_auc"]]
    sens=seed[seed.hierarchy=="detailed"].groupby(["contrast","state"],observed=False)[seed_col].agg(["min","max"]).reset_index().merge(python_auc,on=["contrast","state"]).merge(null_d,on=["contrast","state"])
    fig,ax=plt.subplots(figsize=(8,4.5)); x=np.arange(sens.shape[0]); ax.vlines(x,sens["min"],sens["max"],color="#4C78A8"); ax.scatter(x,sens.python_augur_inspired_auc,label="Python approximation mean",s=18); ax.scatter(x,sens.null_auc,label="Python shuffled-null mean",s=18)
    ax.set_xticks(x,[f"{r.contrast.split('_')[-1]}:{r.state}" for r in sens.itertuples()],rotation=70,ha="right"); ax.set_ylabel("Augur-inspired Python AUC"); ax.set_title("Supplementary: Python implementation sensitivity and non-biological shuffled null"); ax.legend(); fig.tight_layout()
    for ext in ["png","svg"]: fig.savefig(out/"figures"/f"05_augur_sensitivity_and_null.{ext}")
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,2.8)); ax.axis("off")
    ax.text(.2,.55,"Official R Augur\nTreatment separability\nAUC from balanced cell CV",ha="center",va="center",bbox=dict(boxstyle="round",fc="#DCEAF7"))
    ax.text(.8,.55,"ScGeo\nDisplacement geometry\nRepresentation stability",ha="center",va="center",bbox=dict(boxstyle="round",fc="#F8E1C4"))
    ax.annotate("Complementary evidence; not the same estimand",xy=(.7,.55),xytext=(.3,.55),arrowprops=dict(arrowstyle="<->"),ha="center",va="bottom")
    ax.text(.5,.12,"Abundance and within-state distribution remain separate descriptive axes",ha="center")
    fig.tight_layout()
    for ext in ["png","svg"]: fig.savefig(out/"figures"/f"06_augur_scgeo_explanatory_schematic.{ext}")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for contrast, marker in [("PBS_vs_TNF", "o"), ("PBS_vs_LPS", "s")]:
        g = sensitivity[sensitivity.contrast == contrast]
        axes[0].scatter(g.official_augur_auc, g.python_augur_inspired_auc, marker=marker, label=contrast)
        for r in g.itertuples(): axes[0].annotate(r.state, (r.official_augur_auc, r.python_augur_inspired_auc), fontsize=6)
    limits = [0.48, 0.72]; axes[0].plot(limits, limits, color="grey", ls="--", lw=.8)
    axes[0].set(xlim=limits, ylim=limits, xlabel="Official R Augur AUC", ylabel="Augur-inspired Python AUC", title="AUC implementation sensitivity")
    axes[0].legend()
    ordered = sensitivity.sort_values(["contrast", "jaccard_overlap"])
    axes[1].barh(np.arange(len(ordered)), ordered.jaccard_overlap, color="#72B7B2")
    axes[1].set_yticks(np.arange(len(ordered)), [f"{r.contrast.split('_')[-1]}: {r.state}" for r in ordered.itertuples()])
    axes[1].set(xlabel="Variable-feature Jaccard overlap", xlim=(0, 1), title="Official LOESS residuals vs raw sparse variance")
    fig.suptitle("Supplementary implementation sensitivity; seven marker-inferred states per contrast")
    fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(out / "figures" / f"06_python_augur_inspired_sensitivity.{ext}")
    plt.close(fig)


def run_comparison() -> dict[str, Any]:
    config, p = load_config(), paths(); validate(config, p); ensure_dirs(p["output"])
    stage="06_scgeo_augur_comparison"; resources=ResourceLog(stage,config); started=time.perf_counter(); before={k:sha256(p[k]) for k in ["compact","representation"]}
    with resources.operation("assemble_full_comparator_table"):
        table, correlations=build_comparator(config,p)
        official_summary=official_augur_summary(config,p)
        sensitivity=implementation_sensitivity(config)
        methods=comparator_methods_table()
        seed=pd.read_csv(p["output"] / "controls" / "05_augur_fixed_seed_sensitivity.csv")
        null=pd.read_csv(p["output"] / "controls" / "05_augur_shuffled_null_summary.csv")
    with resources.operation("write_comparison_tables_and_figures"):
        atomic_csv(official_summary,p["output"] / "comparator" / "05_official_r_augur_state_summary.csv")
        atomic_csv(table,p["output"] / "comparator" / "06_scgeo_augur_full_comparator.csv")
        atomic_csv(correlations,p["output"] / "comparator" / "06_scgeo_augur_rank_correlations.csv")
        atomic_csv(sensitivity,p["output"] / "comparator" / "06_official_python_implementation_sensitivity.csv")
        atomic_csv(methods,p["output"] / "comparator" / "06_augur_comparator_methods.csv")
        atomic_csv(official_summary,p["output"] / "figure_sources" / "05_official_r_augur_state_summary.csv")
        atomic_csv(table,p["output"] / "figure_sources" / "06_scgeo_augur_full_comparator.csv")
        atomic_csv(correlations,p["output"] / "figure_sources" / "06_scgeo_augur_rank_correlations.csv")
        atomic_csv(sensitivity,p["output"] / "figure_sources" / "06_official_python_implementation_sensitivity.csv")
        atomic_csv(methods,p["output"] / "figure_sources" / "06_augur_comparator_methods.csv")
        figures(p["output"],table,correlations,seed,null,sensitivity)
        alt=("Deterministic Phase 3C comparator alt text. The main ranking, ScGeo scatter, method-evidence matrix, and explanatory schematic use official R Augur 1.0.3 at commit b252b84e4af687d9817813b1db409267eb44ec3f as the primary external comparator. Official Augur AUC denotes treatment separability; ScGeo denotes displacement geometry and representation stability, and abundance and distribution remain separate descriptive axes. A supplementary implementation-sensitivity panel compares official AUC with the Augur-inspired Python approximation without averaging, merging, or selecting between them, and reports feature-set Jaccard overlap for official LOESS-residual versus Python raw-sparse-variance selection. The fixed-seed and shuffled-null panel is explicitly Python-only. Rank correlations involve seven marker-inferred states per contrast. All uncertainty is computational cell-resampling stability and all evidence is descriptive_only.")
        (p["output"] / "alt_text" / "05_06_augur_scgeo_comparison.txt").write_text(alt+"\n",encoding="utf-8")
        captions=p["output"] / "captions"; captions.mkdir(parents=True,exist_ok=True)
        (captions / "05_06_augur_scgeo_comparison.md").write_text(
            "Official R Augur 1.0.3 (commit b252b84e4af687d9817813b1db409267eb44ec3f) is the primary external comparator. "
            "The Augur-inspired Python approximation matched headline subsampling and random-forest settings but used raw sparse variance instead of official LOESS-residual variable-feature selection; its forest engine, seed schedule, and CV construction also differed. "
            "Python results are shown only as supplementary implementation sensitivity and are never averaged or merged with official AUC. Rank correlations involve seven marker-inferred states per contrast. All results are descriptive_only.\n",
            encoding="utf-8",
        )
    after={k:sha256(p[k]) for k in ["compact","representation"]}
    if before!=after: raise RuntimeError("Immutable input changed")
    payload={"status":"passed","canonical_comparator":"Official R Augur 1.0.3","canonical_commit":config["canonical_comparator"]["commit"],"python_role":"supplementary_implementation_sensitivity","rank_scope":"seven_marker_inferred_states_per_contrast","n_comparator_rows":int(table.shape[0]),"n_rank_correlations":int(correlations.shape[0]),"runtime_seconds":time.perf_counter()-started,"input_sha256_before":before,"input_sha256_after":after,"combined_score_created":False,"official_python_auc_averaged_or_merged":False,"inferential_p_values":False}
    official_implementation={**config["canonical_comparator"],"implementation":"Official R Augur","role":"primary_external_comparator"}
    atomic_json(payload,p["output"] / "audit" / "06_scgeo_augur_comparison_summary.json"); write_resources(resources,p["output"]); metadata(stage,config,p,resources,payload,implementation=official_implementation); close_and_collect(); return payload


def relabel_existing_python_outputs() -> dict[str, Any]:
    """Apply documentation-only labels to existing Python outputs in place."""
    config, p = load_config(), paths(); ensure_dirs(p["output"])
    relpaths = [
        "comparator/05_augur_subsample_auc.csv.gz",
        "comparator/05_augur_state_summary.csv",
        "figure_sources/05_augur_state_summary.csv",
        "figure_sources/05_augur_sensitivity_and_null.csv",
        "controls/05_augur_fixed_seed_sensitivity.csv",
        "controls/05_augur_shuffled_subsample_auc.csv.gz",
        "controls/05_augur_shuffled_null_summary.csv",
        "controls/05_augur_exclude_ambiguous_hspc.csv",
        "controls/05_augur_hierarchy_sensitivity.csv",
    ]
    checked=[]
    for relpath in relpaths:
        path=p["output"] / relpath
        frame=pd.read_csv(path)
        numeric_before=frame.select_dtypes(include=[np.number]).to_numpy(copy=True)
        labelled=python_sensitivity_export(frame)
        numeric_after=labelled.select_dtypes(include=[np.number]).to_numpy(copy=True)
        if numeric_before.shape != numeric_after.shape or not np.allclose(numeric_before,numeric_after,equal_nan=True):
            raise RuntimeError(f"Numeric Python output changed while relabeling {relpath}")
        atomic_csv(labelled,path)
        checked.append(relpath)
    summary=pd.read_csv(p["output"] / "comparator" / "05_augur_state_summary.csv")
    atomic_csv(summary,p["output"] / "comparator" / "05_python_augur_inspired_state_summary.csv")
    for folder in ["metadata","version_records"]:
        for stem in ["05_augur_comparator"]:
            path=p["output"] / folder / f"{stem}_{'metadata' if folder == 'metadata' else 'versions'}.json"
            payload=json.loads(path.read_text(encoding="utf-8")); payload["implementation"]=config["implementation"]
            payload.pop("augur_parameters",None)
            if folder == "metadata": payload["python_augur_inspired_headline_settings"]=config["augur_defaults"]
            payload["manuscript_role"]="supplementary_implementation_sensitivity"; payload["inference_status"]="descriptive_only"
            atomic_json(payload,path)
    audit_path=p["output"] / "audit" / "05_augur_summary.json"
    audit=json.loads(audit_path.read_text(encoding="utf-8")); audit["implementation"]="Augur-inspired Python approximation"; audit["manuscript_role"]="supplementary_implementation_sensitivity"
    atomic_json(audit,audit_path)
    return {"status":"passed","tables_relabelled":checked,"numeric_values_changed":False,"inference_status":"descriptive_only"}


def regenerate_comparator_documentation() -> dict[str, Any]:
    config, p = load_config(), paths()
    immutable = {
        "compact_h5ad": p["compact"],
        "representation_h5ad": p["representation"],
        "scgeo_state_evidence": p["output"] / "scgeo" / "04_full_state_evidence.csv",
        "official_augur_state_auc": _configured_result_path(config["canonical_comparator"]["state_auc_path"]),
        "official_augur_subsample_auc": _configured_result_path(config["canonical_comparator"]["subsample_auc_path"]),
    }
    before={name:sha256(path) for name,path in immutable.items()}
    python_summary=p["output"] / "comparator" / "05_augur_state_summary.csv"
    python_before=sha256(python_summary)
    relabel=relabel_existing_python_outputs()
    from scripts.compare_gse249479_official_augur import main as refresh_official_python_comparison
    refresh_official_python_comparison()
    comparison=run_comparison()
    after={name:sha256(path) for name,path in immutable.items()}
    if before != after: raise RuntimeError("A numerical input changed during documentation regeneration")
    report={
        "status":"passed",
        "canonical_comparator":"Official R Augur 1.0.3",
        "canonical_commit":config["canonical_comparator"]["commit"],
        "python_implementation":"Augur-inspired Python approximation",
        "python_role":"supplementary_implementation_sensitivity",
        "python_headline_parameters_matched":config["implementation"]["headline_parameters_matched"],
        "python_documented_differences":config["implementation"]["documented_differences"],
        "rank_scope":"seven_marker_inferred_states_per_contrast",
        "python_output_relabeling":relabel,
        "canonical_comparison":comparison,
        "immutable_numerical_input_sha256_before":before,
        "immutable_numerical_input_sha256_after":after,
        "python_summary_sha256_before_documentation_relabel":python_before,
        "python_summary_sha256_after_documentation_relabel":sha256(python_summary),
        "python_numeric_values_changed":False,
        "official_python_auc_averaged_merged_or_selected":False,
        "scgeo_rerun":False,
        "python_retroactive_tuning":False,
        "inference_status":"descriptive_only",
    }
    atomic_json(report,p["output"] / "audit" / "06_comparator_documentation_correction.json")
    return report
