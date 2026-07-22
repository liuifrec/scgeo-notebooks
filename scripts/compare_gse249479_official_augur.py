"""Compare official R Augur with immutable Python and ScGeo Phase 3C outputs."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.gse249479_augur_common import load_config as load_phase3c_config
from scripts.gse249479_augur_common import load_sparse_hvg, paths as phase3c_paths, variable_features


CONFIG = ROOT / "configs/gse249479_official_augur_validation_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(payload: dict, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def agreement_label(rho: float, thresholds: dict) -> str:
    if not np.isfinite(rho):
        return "unavailable"
    if rho >= float(thresholds["high_spearman_min"]):
        return "high"
    if rho >= float(thresholds["moderate_spearman_min"]):
        return "moderate"
    return "low"


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    output = ROOT / config["validation_output"]
    official_path = output / "official_augur_state_auc.csv"
    python_path = ROOT / config["existing_python_summary"]
    scgeo_path = ROOT / config["existing_scgeo_comparator"]
    source_paths = [python_path, scgeo_path]
    baseline = {str(path.relative_to(ROOT)): sha256(path) for path in source_paths}

    official = pd.read_csv(official_path)
    python = pd.read_csv(python_path)
    python_auc_column = "python_augur_inspired_auc" if "python_augur_inspired_auc" in python.columns else "augur_auc"
    python = python[python["hierarchy"].eq("detailed")][["contrast", "state", python_auc_column]].rename(columns={python_auc_column: "python_augur_inspired_auc"})
    scgeo = pd.read_csv(scgeo_path)
    if "hierarchy" in scgeo.columns:
        scgeo = scgeo[scgeo["hierarchy"].eq("detailed")]
    scgeo = scgeo[["contrast", "state", "scgeo_normalized_displacement"]].drop_duplicates(["contrast", "state"])
    merged = official.merge(python, on=["contrast", "state"], validate="one_to_one").merge(scgeo, on=["contrast", "state"], validate="one_to_one")
    expected = len(config["contrasts"]) * len(config["detailed_state_order"])
    if len(merged) != expected:
        raise RuntimeError(f"Expected {expected} state/contrast rows, found {len(merged)}")
    merged["official_augur_rank"] = merged.groupby("contrast")["official_augur_auc"].rank(method="average", ascending=False)
    merged["python_augur_inspired_rank"] = merged.groupby("contrast")["python_augur_inspired_auc"].rank(method="average", ascending=False)
    merged["scgeo_displacement_rank"] = merged.groupby("contrast")["scgeo_normalized_displacement"].rank(method="average", ascending=False)
    merged["absolute_auc_difference"] = (merged["official_augur_auc"] - merged["python_augur_inspired_auc"]).abs()
    merged["auc_difference_official_minus_python"] = merged["official_augur_auc"] - merged["python_augur_inspired_auc"]
    merged["official_role"] = "primary_external_comparator"
    merged["python_role"] = "supplementary_implementation_sensitivity"
    merged["rank_scope"] = "seven_marker_inferred_states_per_contrast"
    merged["inference_status"] = "descriptive_only"
    merged["uncertainty_status"] = "computational_cell_resampling_stability_only"

    official_features = pd.read_csv(output / "official_augur_variance_selected_features.csv")
    phase3c_config = load_phase3c_config()
    X, obs, genes = load_sparse_hvg(phase3c_config, phase3c_paths(phase3c_config))
    condition = obs["condition"].astype(str).to_numpy()
    detailed = obs["state_detailed"].astype(str).to_numpy()
    feature_rows = []
    for treated in config["contrasts"]:
        contrast = f"PBS_vs_{treated}"
        for state in config["detailed_state_order"]:
            pool = np.flatnonzero(np.isin(condition, ["PBS", treated]) & (detailed == state))
            selected_idx = variable_features(X, pool, 0.5)
            python_set = set(genes[selected_idx].astype(str))
            official_set = set(official_features.loc[
                official_features["contrast"].eq(contrast) & official_features["state"].eq(state), "gene_id"
            ].astype(str))
            intersection = len(python_set & official_set)
            union = len(python_set | official_set)
            feature_rows.append({
                "contrast": contrast,
                "state": state,
                "official_selected_features": len(official_set),
                "python_selected_features": len(python_set),
                "intersection_features": intersection,
                "jaccard_overlap": intersection / union if union else np.nan,
                "official_selection_method": "LOESS_residual_of_mean_over_sd_after_1pct_tail_trim",
                "python_selection_method": "upper_half_by_raw_sparse_variance",
                "inference_status": "descriptive_only",
            })
    feature_overlap = pd.DataFrame(feature_rows)
    merged = merged.merge(feature_overlap, on=["contrast", "state", "inference_status"], validate="one_to_one")

    rows = []
    thresholds = config["agreement_decision_thresholds"]
    for contrast, frame in merged.groupby("contrast", sort=False):
        official_python = float(spearmanr(frame["official_augur_auc"], frame["python_augur_inspired_auc"]).statistic)
        official_scgeo = float(spearmanr(frame["official_augur_auc"], frame["scgeo_normalized_displacement"]).statistic)
        rows.append({
            "contrast": contrast,
            "official_vs_python_spearman": official_python,
            "official_vs_python_agreement": agreement_label(official_python, thresholds),
            "official_vs_scgeo_spearman": official_scgeo,
            "mean_absolute_auc_difference": float(frame["absolute_auc_difference"].mean()),
            "median_absolute_auc_difference": float(frame["absolute_auc_difference"].median()),
            "max_absolute_auc_difference": float(frame["absolute_auc_difference"].max()),
            "state_at_max_absolute_auc_difference": str(frame.loc[frame["absolute_auc_difference"].idxmax(), "state"]),
            "mean_feature_selection_jaccard": float(frame["jaccard_overlap"].mean()),
            "min_feature_selection_jaccard": float(frame["jaccard_overlap"].min()),
            "max_feature_selection_jaccard": float(frame["jaccard_overlap"].max()),
            "n_states": int(len(frame)),
            "inference_status": "descriptive_only",
        })
    summary = pd.DataFrame(rows)
    all_high = bool(summary["official_vs_python_agreement"].eq("high").all())
    any_low = bool(summary["official_vs_python_agreement"].eq("low").any())
    if all_high:
        decision = "official_primary_python_supplementary_reproducibility_sensitivity"
    elif not any_low:
        decision = "official_primary_python_supplementary_implementation_sensitivity_moderate_agreement"
    else:
        decision = "official_primary_python_supplementary_implementation_sensitivity_nonuniform_agreement"

    after = {str(path.relative_to(ROOT)): sha256(path) for path in source_paths}
    if baseline != after:
        raise RuntimeError("An existing Phase 3C output changed during validation")
    output.mkdir(parents=True, exist_ok=True)
    atomic_csv(merged.sort_values(["contrast", "official_augur_rank"]), output / "official_python_scgeo_state_comparison.csv")
    atomic_csv(summary, output / "official_python_scgeo_agreement_summary.csv")
    atomic_csv(feature_overlap, output / "official_python_feature_selection_overlap.csv")
    report = {
        "inference_status": "descriptive_only",
        "biological_replicates": False,
        "official_augur_status": "installed_and_run",
        "official_source": config["official_source"],
        "official_commit": config["official_commit"],
        "official_version": config["official_version"],
        "agreement_thresholds": thresholds,
        "decision": decision,
        "existing_phase3c_outputs_unchanged": True,
        "existing_phase3c_checksums": after,
        "documented_implementation_differences": [
            "official uses R randomForest while the Augur-inspired Python approximation uses sklearn RandomForestClassifier",
            "official resets each subsample seed to 1 through 50 while the Python sensitivity implementation uses prespecified fixed seed blocks",
            "official uses rsample stratified v-fold assignment while the Python sensitivity implementation uses sklearn StratifiedKFold",
            "official variance selection uses LOESS residuals of mean-over-SD after 1% tail trimming; the Python approximation uses the upper half by raw sparse variance",
            "official uses floor(n_features * 0.5) for random feature sampling while the Python approximation uses ceil(n_features * 0.5)"
        ],
        "canonical_comparator": "Official R Augur 1.0.3",
        "python_role": "supplementary_implementation_sensitivity",
        "rank_scope": "seven_marker_inferred_states_per_contrast",
        "interpretation": "Official Augur AUC is treatment separability; ScGeo is displacement geometry. The Augur-inspired Python approximation is a separate implementation sensitivity and is not averaged or merged with official AUC. Rank correlations involve seven marker-inferred states and are descriptive; they do not establish method equivalence or biological replication.",
        "summary": summary.to_dict(orient="records"),
    }
    atomic_json(report, output / "official_augur_validation_report.json")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
