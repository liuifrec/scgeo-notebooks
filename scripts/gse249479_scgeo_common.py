"""Descriptive-only ScGeo treatment geometry for GSE249479."""

from __future__ import annotations

import inspect
import json
import math
import os
import platform
import time
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist
from scipy.stats import wasserstein_distance

from scripts.gse249479_representation_common import (
    ROOT,
    ResourceLog,
    active_branch,
    atomic_csv,
    atomic_json,
    close_and_collect,
    ensure_tree,
    git_commit,
    package_versions,
    sha256,
    utc_now,
)


CONFIG_PATH = ROOT / "configs/gse249479_scgeo_v1.json"


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def configured_paths(config: dict[str, Any] | None = None) -> dict[str, Path]:
    config = load_config() if config is None else config
    output = Path(os.environ.get("SCGEO_GSE249479_OUTPUT_DIR", config["default_output_dir"]))
    if not output.is_absolute():
        output = ROOT / output
    return {
        "representation": Path(os.environ.get("SCGEO_GSE249479_REPRESENTATION_H5AD", config["default_representation_h5ad"])).resolve(),
        "compact": Path(os.environ.get("SCGEO_GSE249479_COMPACT_H5AD", config["default_compact_h5ad"])).resolve(),
        "output": output.resolve(),
        "source_repo": Path(os.environ.get("SCGEO_SOURCE_REPO", ROOT.parent / "scgeo")).resolve(),
    }


def ensure_scgeo_tree(output: Path) -> None:
    ensure_tree(output)
    for name in ["scgeo", "controls"]:
        (output / name).mkdir(parents=True, exist_ok=True)


def write_resources(resources: ResourceLog, output: Path) -> None:
    frame = pd.DataFrame(resources.rows)
    frame["inference_status"] = "descriptive_only"
    frame["resampling_status"] = "descriptive_cell_resampling_only_where_applicable"
    atomic_csv(frame, output / "audit" / f"{resources.stage}_resource_log.csv")
    combined = output / "audit" / "04_scgeo_resource_log.csv"
    if combined.exists():
        frame = pd.concat([pd.read_csv(combined), frame], ignore_index=True)
    atomic_csv(frame, combined)


def validate_environment(config: dict[str, Any], paths: dict[str, Path]) -> None:
    if active_branch() != config["required_branch"]:
        raise RuntimeError(f"Required branch {config['required_branch']!r}; observed {active_branch()!r}")
    if git_commit(paths["source_repo"]) != config["frozen_scgeo_commit"]:
        raise RuntimeError("Frozen ScGeo commit mismatch")
    for key in ["representation", "compact"]:
        if not paths[key].is_file():
            raise FileNotFoundError(paths[key])
    if sha256(paths["representation"]) != config["expected_representation_sha256"]:
        raise RuntimeError("Representation H5AD checksum mismatch")
    if sha256(paths["compact"]) != config["expected_compact_sha256"]:
        raise RuntimeError("Compact H5AD checksum mismatch")


def load_representation(config: dict[str, Any], paths: dict[str, Path]) -> ad.AnnData:
    obj = ad.read_h5ad(paths["representation"])
    if obj.n_obs != 34432 or obj.n_vars != 0:
        raise RuntimeError(f"Unexpected representation object shape {obj.shape}")
    required = config["representations"]["all_quantitative"]
    missing = [rep for rep in required if rep not in obj.obsm]
    if missing:
        raise RuntimeError(f"Missing quantitative representations: {missing}")
    if set(obj.obs["inference_status"].astype(str)) != {"descriptive_only"}:
        raise RuntimeError("Representation object does not retain descriptive_only")
    obj.obs["state_detailed"] = pd.Categorical(
        obj.obs["marker_inferred_label"].astype(str).map(config["detailed_state_mapping"]),
        categories=config["detailed_state_order"], ordered=True,
    )
    obj.obs["state_coarse"] = pd.Categorical(
        obj.obs["state_detailed"].astype(str).map(config["coarse_state_mapping"]),
        categories=config["coarse_state_order"], ordered=True,
    )
    return obj


def scgeo_defaults_snapshot(sg: Any) -> dict[str, Any]:
    def defaults(fn: Any) -> dict[str, Any]:
        out = {}
        for name, param in inspect.signature(fn).parameters.items():
            if param.default is not inspect.Parameter.empty:
                value = param.default
                out[name] = list(value) if isinstance(value, tuple) else value
        return out
    return {
        "robust_shift": defaults(sg.tl.robust_shift),
        "representation_stability": defaults(sg.tl.representation_stability),
        "local_geometry_stability": defaults(sg.tl.local_geometry_stability),
    }


def robust_output_rows(out: dict[str, Any], rep: str, states: list[str], hierarchy: str, contrast: str) -> pd.DataFrame:
    rows = []
    for state in states:
        rec = out.get("by", {}).get(state, {})
        ci = rec.get("bootstrap_magnitude_ci95", [np.nan, np.nan])
        n0, n1 = int(rec.get("n_cells0", 0)), int(rec.get("n_cells1", 0))
        rows.append({
            "contrast": f"PBS_vs_{contrast}", "hierarchy": hierarchy, "representation": rep, "state": state,
            "status": "usable" if min(n0, n1) >= 20 else "insufficient_coverage",
            "coverage_status_reason": "at_least_20_cells_per_condition" if min(n0, n1) >= 20 else "fewer_than_20_cells_in_at_least_one_condition",
            "n_cells_PBS": n0, f"n_cells_{contrast}": n1,
            "robust_displacement_magnitude": rec.get("delta_norm", np.nan),
            "normalized_effect": rec.get("normalized_delta_norm", np.nan),
            "cell_bootstrap_magnitude_low": ci[0], "cell_bootstrap_magnitude_high": ci[1],
            "direction_stability": rec.get("bootstrap_directional_resultant_length", np.nan),
            "direction_sign_stability": rec.get("direction_stability", np.nan),
            "outlier_cosine_to_mean": (rec.get("outlier_sensitivity") or {}).get("cosine_to_mean", np.nan),
            "resampling_status": "descriptive_cell_resampling_only", "inference_status": "descriptive_only",
        })
    return pd.DataFrame(rows)


def stability_rows(out: dict[str, Any], contrast: str, hierarchy: str) -> pd.DataFrame:
    frame = out["per_rep_state"].copy().rename(columns={
        "rep": "representation", "node": "state", "n_cells0": "n_cells_PBS",
        "delta_norm": "robust_displacement_magnitude", "normalized_delta_norm": "normalized_effect",
        "magnitude_ci95_low": "cell_bootstrap_magnitude_low", "magnitude_ci95_high": "cell_bootstrap_magnitude_high",
        "directional_resultant_length": "direction_stability",
    })
    frame[f"n_cells_{contrast}"] = frame.pop("n_cells1")
    frame["contrast"] = f"PBS_vs_{contrast}"
    frame["hierarchy"] = hierarchy
    frame["coverage_status_reason"] = np.where(frame["status"].eq("usable"), "at_least_20_cells_per_condition", frame["status"])
    frame["resampling_status"] = "descriptive_cell_resampling_only"
    frame["inference_status"] = "descriptive_only"
    return frame


def add_role(frame: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    role = {}
    for key in ["primary", "dimensional_sensitivity", "exploratory_sensitivity"]:
        role.update({rep: key for rep in config["representations"][key]})
    frame = frame.copy()
    frame["representation_role"] = frame["representation"].map(role)
    return frame


def representation_agreement(frame: pd.DataFrame, contrast: str) -> pd.DataFrame:
    usable = frame[(frame["hierarchy"] == "detailed") & frame["status"].eq("usable")].copy()
    usable["magnitude_rank"] = usable.groupby("representation")["normalized_effect"].rank(method="average", ascending=False)
    rows = []
    for state, group in usable.groupby("state", observed=False):
        vals = group["normalized_effect"].dropna().to_numpy(float)
        ranks = group["magnitude_rank"].dropna().to_numpy(float)
        rows.append({
            "contrast": f"PBS_vs_{contrast}", "state": state, "n_usable_representations": int(group.shape[0]),
            "normalized_effect_median_all_representations": float(np.median(vals)) if vals.size else np.nan,
            "normalized_effect_iqr_all_representations": float(np.subtract(*np.percentile(vals, [75, 25]))) if vals.size else np.nan,
            "magnitude_rank_std_all_representations": float(np.std(ranks)) if ranks.size else np.nan,
            "representation_agreement_status": "coordinate_safe_scalar_summary_only",
            "inference_status": "descriptive_only",
        })
    return pd.DataFrame(rows)


def leave_one_rows(frame: pd.DataFrame, contrast: str) -> pd.DataFrame:
    usable = frame[(frame["hierarchy"] == "detailed") & frame["status"].eq("usable")].copy()
    rows = []
    for state, group in usable.groupby("state", observed=False):
        by = group.set_index("representation")["normalized_effect"].to_dict()
        full = float(np.nanmedian(list(by.values())))
        for omitted in by:
            remaining = [v for k, v in by.items() if k != omitted and np.isfinite(v)]
            loo = float(np.median(remaining)) if remaining else np.nan
            rows.append({
                "contrast": f"PBS_vs_{contrast}", "state": state, "omitted_representation": omitted,
                "full_median_normalized_effect": full, "leave_one_median_normalized_effect": loo,
                "absolute_change": abs(loo - full) if np.isfinite(loo) else np.nan,
                "relative_change": abs(loo - full) / max(abs(full), 0.25) if np.isfinite(loo) else np.nan,
                "inference_status": "descriptive_only",
            })
    return pd.DataFrame(rows)


def run_contrast(treated: str) -> dict[str, Any]:
    import scgeo as sg

    config, paths = load_config(), configured_paths()
    validate_environment(config, paths)
    ensure_scgeo_tree(paths["output"])
    if treated not in config["contrasts"]:
        raise ValueError(treated)
    slug = treated.lower()
    stage = f"04{'a' if treated == 'TNF' else 'b'}_scgeo_pbs_{slug}"
    resources = ResourceLog(stage, config)
    before = {key: sha256(paths[key]) for key in ["representation", "compact"]}
    started = time.perf_counter()

    with resources.operation("load_expression_free_representation_object"):
        obj = load_representation(config, paths)

    with resources.operation("primary_detailed_frozen_scgeo_stability"):
        primary = sg.tl.representation_stability(
            obj, reps=config["representations"]["primary"], node_key="state_detailed",
            condition_key=config["condition_key"], group0=config["group0"], group1=treated,
            seed=int(config["random_seed"]), store_key=f"primary_{slug}",
        )
        detailed = stability_rows(primary, treated, "detailed")

    with resources.operation("pca20_and_diffusion_frozen_robust_shift"):
        for offset, rep in enumerate(config["representations"]["dimensional_sensitivity"] + config["representations"]["exploratory_sensitivity"]):
            out = sg.tl.robust_shift(
                obj, rep=rep, condition_key=config["condition_key"], group0=config["group0"], group1=treated,
                by="state_detailed", seed=int(config["random_seed"]) + offset + 1, store_key=f"{slug}_{rep}",
            )
            detailed = pd.concat([detailed, robust_output_rows(out, rep, config["detailed_state_order"], "detailed", treated)], ignore_index=True, sort=False)

    with resources.operation("coarse_hierarchy_frozen_scgeo_stability"):
        coarse = sg.tl.representation_stability(
            obj, reps=config["representations"]["all_quantitative"], node_key="state_coarse",
            condition_key=config["condition_key"], group0=config["group0"], group1=treated,
            seed=int(config["random_seed"]) + 20, store_key=f"coarse_{slug}",
        )
        coarse_rows = stability_rows(coarse, treated, "coarse")

    detailed = add_role(detailed, config)
    coarse_rows = add_role(coarse_rows, config)
    primary_consensus = primary["consensus_state"].copy().rename(columns={"node": "state", "consensus_label": "primary_consensus_status"})
    primary_consensus["contrast"] = f"PBS_vs_{treated}"
    primary_consensus["inference_status"] = "descriptive_only"
    coarse_consensus = coarse["consensus_state"].copy().rename(columns={"node": "state", "consensus_label": "all_representation_consensus_status"})
    coarse_consensus["contrast"] = f"PBS_vs_{treated}"
    coarse_consensus["hierarchy"] = "coarse"
    coarse_consensus["inference_status"] = "descriptive_only"
    agreement = representation_agreement(detailed, treated)
    loo = leave_one_rows(detailed, treated)

    with resources.operation("write_contrast_evidence"):
        base = paths["output"] / "scgeo"
        atomic_csv(detailed, base / f"04_{slug}_detailed_representation_state_evidence.csv")
        atomic_csv(coarse_rows, base / f"04_{slug}_coarse_representation_state_evidence.csv")
        atomic_csv(primary_consensus, base / f"04_{slug}_primary_consensus.csv")
        atomic_csv(coarse_consensus, base / f"04_{slug}_coarse_consensus.csv")
        atomic_csv(agreement, base / f"04_{slug}_representation_agreement.csv")
        atomic_csv(loo, paths["output"] / "controls" / f"04_{slug}_leave_one_representation.csv")
        atomic_json({
            "contrast": f"PBS_vs_{treated}", "inference_status": "descriptive_only",
            "sample_key": None, "biological_replicates": False,
            "resampling_status": "descriptive_cell_resampling_only",
            "scgeo_defaults_observed": scgeo_defaults_snapshot(sg),
            "primary_coverage": primary["coverage_summary"], "coarse_coverage": coarse["coverage_summary"],
            "primary_warnings": primary["warnings"], "coarse_warnings": coarse["warnings"],
        }, paths["output"] / "audit" / f"04_{slug}_scgeo_summary.json")

    after = {key: sha256(paths[key]) for key in ["representation", "compact"]}
    if after != before:
        raise RuntimeError("An immutable input changed")
    summary = {
        "status": "passed", "contrast": f"PBS_vs_{treated}", "inference_status": "descriptive_only",
        "n_detailed_rows": int(detailed.shape[0]), "n_coarse_rows": int(coarse_rows.shape[0]),
        "runtime_seconds": time.perf_counter() - started, "input_sha256_before": before, "input_sha256_after": after,
        "peak_cpu_rss_gib": resources.peak_rss_gib,
    }
    write_resources(resources, paths["output"])
    write_metadata(stage, config, paths, resources, summary, sg)
    del obj
    close_and_collect()
    return summary


def geometric_median(X: np.ndarray, tol: float = 1e-7, max_iter: int = 512) -> np.ndarray:
    y = np.median(X, axis=0)
    for _ in range(max_iter):
        d = np.linalg.norm(X - y, axis=1)
        if np.all(d <= 1e-12):
            return y
        w = 1.0 / np.maximum(d, 1e-12)
        nxt = (X * w[:, None]).sum(axis=0) / w.sum()
        if np.linalg.norm(nxt - y) <= tol * max(1.0, np.linalg.norm(y)):
            return nxt
        y = nxt
    return y


def normalized_shift(X0: np.ndarray, X1: np.ndarray) -> tuple[float, float]:
    c0, c1 = geometric_median(X0), geometric_median(X1)
    magnitude = float(np.linalg.norm(c1 - c0))
    r0 = float(np.median(np.linalg.norm(X0 - c0, axis=1)))
    r1 = float(np.median(np.linalg.norm(X1 - c1, axis=1)))
    scale = math.sqrt((r0 * r0 + r1 * r1) / 2.0)
    return magnitude, magnitude / scale if scale > 1e-12 else np.nan


def balanced_subsampling(obj: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    rows = []
    repeats = int(config["balanced_subsampling"]["repeats"])
    for treated in config["contrasts"]:
        for state in config["detailed_state_order"]:
            idx0 = np.flatnonzero((obj.obs["condition"].astype(str).to_numpy() == "PBS") & (obj.obs["state_detailed"].astype(str).to_numpy() == state))
            idx1 = np.flatnonzero((obj.obs["condition"].astype(str).to_numpy() == treated) & (obj.obs["state_detailed"].astype(str).to_numpy() == state))
            take = min(idx0.size, idx1.size)
            if take < 20:
                rows.append({"contrast": f"PBS_vs_{treated}", "state": state, "status": "insufficient_coverage", "n_balanced_each": take, "inference_status": "descriptive_only"})
                continue
            for repeat in range(repeats):
                rng = np.random.RandomState(int(config["random_seed"]) + repeat * int(config["balanced_subsampling"]["seed_stride"]))
                s0, s1 = rng.choice(idx0, take, replace=False), rng.choice(idx1, take, replace=False)
                for rep in config["representations"]["all_quantitative"]:
                    mag, norm = normalized_shift(np.asarray(obj.obsm[rep][s0], dtype=np.float64), np.asarray(obj.obsm[rep][s1], dtype=np.float64))
                    rows.append({
                        "contrast": f"PBS_vs_{treated}", "state": state, "representation": rep, "repeat": repeat,
                        "n_balanced_each": take, "robust_displacement_magnitude": mag, "normalized_effect": norm,
                        "status": "usable", "resampling_status": "descriptive_balanced_outer_cell_subsampling_point_estimate",
                        "inference_status": "descriptive_only",
                    })
    return pd.DataFrame(rows)


def abundance_tables(obj: ad.AnnData, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, sensitivity = [], []
    conditions = ["PBS", *config["contrasts"]]
    for hierarchy, key, states in [("detailed", "state_detailed", config["detailed_state_order"]), ("coarse", "state_coarse", config["coarse_state_order"])]:
        cond = obj.obs["condition"].astype(str).to_numpy()
        labels = obj.obs[key].astype(str).to_numpy()
        for condition in conditions:
            mask = cond == condition
            total = int(mask.sum())
            for state in states:
                n = int(np.sum(mask & (labels == state)))
                rows.append({"hierarchy": hierarchy, "condition": condition, "state": state, "n_cells": n, "condition_total": total, "proportion": n / total, "inference_status": "descriptive_only"})
        take = min(int(np.sum(cond == c)) for c in conditions)
        for repeat in range(int(config["balanced_subsampling"]["repeats"])):
            rng = np.random.RandomState(int(config["random_seed"]) + 5000 + repeat)
            for condition in conditions:
                idx = np.flatnonzero(cond == condition)
                chosen = rng.choice(idx, take, replace=False)
                for state in states:
                    prop = float(np.mean(labels[chosen] == state))
                    sensitivity.append({"hierarchy": hierarchy, "condition": condition, "state": state, "repeat": repeat, "balanced_total": take, "proportion": prop, "resampling_status": "descriptive_balanced_condition_cell_subsampling", "inference_status": "descriptive_only"})
    return pd.DataFrame(rows), pd.DataFrame(sensitivity)


def distribution_metrics(obj: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    rows = []
    cap = int(config["distribution"]["max_cells_per_condition_state"])
    n_proj = int(config["distribution"]["sliced_wasserstein_projections"])
    cond = obj.obs["condition"].astype(str).to_numpy()
    labels = obj.obs["state_detailed"].astype(str).to_numpy()
    for treated_i, treated in enumerate(config["contrasts"]):
        for state_i, state in enumerate(config["detailed_state_order"]):
            idx0 = np.flatnonzero((cond == "PBS") & (labels == state))
            idx1 = np.flatnonzero((cond == treated) & (labels == state))
            take = min(idx0.size, idx1.size, cap)
            if take < 20:
                rows.append({"contrast": f"PBS_vs_{treated}", "state": state, "status": "insufficient_coverage", "n_each": take, "inference_status": "descriptive_only"})
                continue
            rng = np.random.RandomState(int(config["random_seed"]) + 100 * treated_i + state_i)
            s0, s1 = rng.choice(idx0, take, replace=False), rng.choice(idx1, take, replace=False)
            for rep in config["distribution"]["representations"]:
                X0 = np.asarray(obj.obsm[rep][s0], dtype=np.float64)
                X1 = np.asarray(obj.obsm[rep][s1], dtype=np.float64)
                dxy = cdist(X0, X1)
                dxx, dyy = cdist(X0, X0), cdist(X1, X1)
                energy = float(2 * dxy.mean() - dxx.mean() - dyy.mean())
                pilot_n = min(int(config["distribution"]["mmd_bandwidth_pilot_cells"]), take)
                pilot = np.vstack([X0[:pilot_n], X1[:pilot_n]])
                bandwidth = float(np.median(pdist(pilot)))
                gamma = 1.0 / (2.0 * bandwidth * bandwidth) if bandwidth > 1e-12 else 1.0
                mmd2 = float(np.exp(-gamma * dxx**2).mean() + np.exp(-gamma * dyy**2).mean() - 2 * np.exp(-gamma * dxy**2).mean())
                projections = rng.normal(size=(X0.shape[1], n_proj))
                projections /= np.maximum(np.linalg.norm(projections, axis=0, keepdims=True), 1e-12)
                sw = float(np.mean([wasserstein_distance(X0 @ projections[:, j], X1 @ projections[:, j]) for j in range(n_proj)]))
                rows.append({
                    "contrast": f"PBS_vs_{treated}", "state": state, "representation": rep, "status": "usable", "n_each": take,
                    "energy_distance_biased_full_value": energy, "mmd_rbf_squared_biased_full_value": mmd2,
                    "mmd_bandwidth": bandwidth, "sliced_wasserstein_full_value": sw, "n_projections": n_proj,
                    "resampling_status": "descriptive_balanced_cell_subsample_for_bounded_computation", "inference_status": "descriptive_only",
                })
    return pd.DataFrame(rows)


def shuffled_controls(obj: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    rows = []
    cond = obj.obs["condition"].astype(str).to_numpy()
    labels = obj.obs["state_detailed"].astype(str).to_numpy()
    for treated_i, treated in enumerate(config["contrasts"]):
        pair_idx = np.flatnonzero(np.isin(cond, ["PBS", treated]))
        original = cond[pair_idx].copy()
        for repeat in range(int(config["shuffled_condition_control"]["repeats"])):
            rng = np.random.RandomState(int(config["random_seed"]) + 9000 + 100 * treated_i + repeat)
            shuffled = original.copy(); rng.shuffle(shuffled)
            for state in config["detailed_state_order"]:
                local = labels[pair_idx] == state
                i0 = pair_idx[local & (shuffled == "PBS")]
                i1 = pair_idx[local & (shuffled == treated)]
                if min(i0.size, i1.size) < 20:
                    continue
                for rep in config["shuffled_condition_control"]["representations"]:
                    mag, norm = normalized_shift(np.asarray(obj.obsm[rep][i0], dtype=np.float64), np.asarray(obj.obsm[rep][i1], dtype=np.float64))
                    rows.append({
                        "contrast": f"PBS_vs_{treated}", "state": state, "representation": rep, "repeat": repeat,
                        "robust_displacement_magnitude": mag, "normalized_effect": norm,
                        "control_status": "non_biological_shuffled_cell_label_computational_null_no_p_value",
                        "inference_status": "descriptive_only",
                    })
    return pd.DataFrame(rows)


def signature_tables(obj: ad.AnnData, config: dict[str, Any], evidence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    labels = obj.obs["state_detailed"].astype(str)
    for treated in config["contrasts"]:
        for state in config["detailed_state_order"]:
            for signature in config["signatures"]:
                col = f"score_{signature}"
                m0 = (obj.obs["condition"].astype(str) == "PBS") & (labels == state)
                m1 = (obj.obs["condition"].astype(str) == treated) & (labels == state)
                v0, v1 = obj.obs.loc[m0, col].astype(float), obj.obs.loc[m1, col].astype(float)
                rows.append({
                    "contrast": f"PBS_vs_{treated}", "state": state, "signature": signature,
                    "n_cells_PBS": int(v0.size), f"n_cells_{treated}": int(v1.size),
                    "mean_PBS": float(v0.mean()), "mean_treated": float(v1.mean()), "mean_difference_treated_minus_PBS": float(v1.mean() - v0.mean()),
                    "median_PBS": float(v0.median()), "median_treated": float(v1.median()), "median_difference_treated_minus_PBS": float(v1.median() - v0.median()),
                    "zero_fraction_PBS": float((v0 == 0).mean()), "zero_fraction_treated": float((v1 == 0).mean()),
                    "limitation": config["study_hsc_i_limitation"] if signature == "study_hsc_i_top50_subset" else "",
                    "inference_status": "descriptive_only",
                })
    sig = pd.DataFrame(rows)
    primary = evidence[(evidence["hierarchy"] == "detailed") & evidence["representation_role"].eq("primary") & evidence["status"].eq("usable")]
    disp = primary.groupby(["contrast", "state"], observed=False)["normalized_effect"].median().rename("primary_median_normalized_effect").reset_index()
    assoc = sig.merge(disp, on=["contrast", "state"], how="left")
    associations = []
    for (contrast, signature), group in assoc.groupby(["contrast", "signature"], observed=False):
        associations.append({
            "contrast": contrast, "signature": signature, "n_states": int(group.dropna(subset=["primary_median_normalized_effect", "mean_difference_treated_minus_PBS"]).shape[0]),
            "spearman_across_states": group["primary_median_normalized_effect"].corr(group["mean_difference_treated_minus_PBS"], method="spearman"),
            "association_status": "descriptive_across_marker_inferred_states_no_independent_replication",
            "inference_status": "descriptive_only",
        })
    return sig, pd.DataFrame(associations)


def local_geometry(obj: ad.AnnData, config: dict[str, Any], sg: Any) -> dict[str, Any]:
    return sg.tl.local_geometry_stability(
        obj, reps=config["representations"]["all_quantitative"], node_key="state_detailed",
        reference_rep=config["local_geometry"]["reference_rep"], pair_mode="reference",
        seed=int(config["random_seed"]), store_key="phase3b_local_geometry",
    )


def exclusion_sensitivities(evidence: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail = evidence[(evidence["hierarchy"] == "detailed") & evidence["status"].eq("usable")].copy()
    diffusion_rows = []
    ambiguous_rows = []
    for (contrast, state), group in detail.groupby(["contrast", "state"], observed=False):
        vals = group.set_index("representation")["normalized_effect"]
        all5 = float(vals.median())
        no_diff = float(vals.drop(index="X_diffmap", errors="ignore").median())
        primary = float(vals.reindex(config["representations"]["primary"]).median())
        diffusion_rows.append({
            "contrast": contrast, "state": state, "all_five_median_normalized_effect": all5,
            "exclude_diffusion_median_normalized_effect": no_diff, "primary_median_normalized_effect": primary,
            "absolute_change_after_excluding_diffusion": abs(no_diff - all5),
            "sensitivity_status": "coordinate_safe_scalar_sensitivity_not_independent_confirmation",
            "inference_status": "descriptive_only",
        })
    ranked = detail.groupby(["contrast", "state"], observed=False)["normalized_effect"].median().reset_index()
    for contrast, group in ranked.groupby("contrast", observed=False):
        group = group.copy()
        group["rank_with_ambiguous"] = group["normalized_effect"].rank(method="average", ascending=False)
        nonamb = group[group["state"] != "Ambiguous HSPC"].copy()
        nonamb["rank_excluding_ambiguous"] = nonamb["normalized_effect"].rank(method="average", ascending=False)
        rank_map = nonamb.set_index("state")["rank_excluding_ambiguous"].to_dict()
        for row in group.itertuples(index=False):
            ambiguous_rows.append({
                "contrast": contrast, "state": row.state, "median_normalized_effect": row.normalized_effect,
                "rank_with_ambiguous": row.rank_with_ambiguous,
                "rank_excluding_ambiguous": rank_map.get(row.state, np.nan),
                "state_point_estimate_changes_when_ambiguous_excluded": False,
                "sensitivity_status": "statewise_geometry_is_independent_of_other_label_levels_rank_context_changes_only",
                "inference_status": "descriptive_only",
            })
    return pd.DataFrame(diffusion_rows), pd.DataFrame(ambiguous_rows)


def plot_outputs(output: Path, evidence: pd.DataFrame, consensus: pd.DataFrame, abundance: pd.DataFrame, distribution: pd.DataFrame, associations: pd.DataFrame) -> None:
    plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 180, "font.size": 8})
    primary = evidence[(evidence["hierarchy"] == "detailed") & evidence["representation_role"].eq("primary")]
    for contrast in ["PBS_vs_TNF", "PBS_vs_LPS"]:
        sub = primary[primary["contrast"] == contrast]
        pivot = sub.pivot(index="state", columns="representation", values="normalized_effect")
        fig, ax = plt.subplots(figsize=(7.4, 3.8)); im = ax.imshow(pivot, aspect="auto", cmap="viridis")
        ax.set_xticks(range(pivot.shape[1]), pivot.columns, rotation=25, ha="right"); ax.set_yticks(range(pivot.shape[0]), pivot.index)
        ax.set_title(f"{contrast.replace('_', ' ')}: primary ScGeo normalized displacement"); fig.colorbar(im, ax=ax, label="Normalized displacement")
        fig.tight_layout()
        for ext in ["png", "svg"]: fig.savefig(output / "figures" / f"04_{contrast.lower()}_state_evidence_panel.{ext}")
        plt.close(fig)

    all_detail = evidence[evidence["hierarchy"] == "detailed"]
    pivot = all_detail.pivot_table(index=["contrast", "state"], columns="representation", values="normalized_effect")
    fig, ax = plt.subplots(figsize=(8.2, 6.2)); im = ax.imshow(pivot, aspect="auto", cmap="magma")
    ax.set_xticks(range(pivot.shape[1]), pivot.columns, rotation=30, ha="right"); ax.set_yticks(range(pivot.shape[0]), [f"{a}: {b}" for a,b in pivot.index])
    ax.set_title("Representation-by-state robust displacement"); fig.colorbar(im, ax=ax, label="Normalized displacement"); fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(output / "figures" / f"04_representation_state_displacement_heatmap.{ext}")
    plt.close(fig)

    comp = all_detail.pivot_table(index=["contrast", "state"], columns="representation", values="normalized_effect").reset_index()
    comp["primary_median"] = comp[["X_pca30", "X_pca50", "X_scvi"]].median(axis=1)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for contrast, group in comp.groupby("contrast"):
        ax.scatter(group["primary_median"], group["X_diffmap"], label=contrast, s=30)
    ax.set_xlabel("Primary median normalized displacement"); ax.set_ylabel("Diffusion-map normalized displacement"); ax.legend(); ax.set_title("Primary consensus versus diffusion sensitivity"); fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(output / "figures" / f"04_primary_vs_diffusion_sensitivity.{ext}")
    plt.close(fig)

    disp = primary.groupby(["contrast", "state"], observed=False)["normalized_effect"].median().reset_index()
    dist = distribution[distribution["representation"].isin(["X_pca30", "X_pca50", "X_scvi"])].groupby(["contrast", "state"], observed=False)["energy_distance_biased_full_value"].median().reset_index()
    abund = abundance[abundance["hierarchy"] == "detailed"].pivot(index="state", columns="condition", values="proportion").reset_index()
    frames=[]
    for contrast, treated in [("PBS_vs_TNF", "TNF"), ("PBS_vs_LPS", "LPS")]:
        x=disp[disp.contrast==contrast].merge(dist[dist.contrast==contrast],on=["contrast","state"]).merge(abund,on="state")
        x["absolute_abundance_change"]=(x[treated]-x["PBS"]).abs(); frames.append(x)
    tri=pd.concat(frames,ignore_index=True)
    fig, axes=plt.subplots(1,2,figsize=(9,3.8))
    for contrast,g in tri.groupby("contrast"):
        axes[0].scatter(g["absolute_abundance_change"],g["normalized_effect"],label=contrast)
        axes[1].scatter(g["energy_distance_biased_full_value"],g["normalized_effect"],label=contrast)
    axes[0].set(xlabel="Absolute abundance proportion change",ylabel="Primary median displacement"); axes[1].set(xlabel="Primary median energy distance",ylabel="Primary median displacement")
    axes[0].legend(); axes[1].legend(); fig.suptitle("Abundance, displacement, and within-state distribution are separate summaries"); fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(output / "figures" / f"04_abundance_displacement_distribution.{ext}")
    plt.close(fig)

    fig, ax=plt.subplots(figsize=(7.5,4.2))
    for contrast,g in associations.groupby("contrast"):
        ax.scatter(range(g.shape[0]),g["spearman_across_states"],label=contrast)
    ax.axhline(0,color="grey",lw=.7); ax.set_xticks(range(associations.signature.nunique()), associations.signature.unique(),rotation=60,ha="right")
    ax.set_ylabel("Descriptive Spearman across states"); ax.set_title("Signature change versus primary displacement"); ax.legend(); fig.tight_layout()
    for ext in ["png", "svg"]: fig.savefig(output / "figures" / f"04_signature_vs_displacement_association.{ext}")
    plt.close(fig)


def run_cross_summary() -> dict[str, Any]:
    import scgeo as sg

    config, paths = load_config(), configured_paths()
    validate_environment(config, paths); ensure_scgeo_tree(paths["output"])
    stage = "04c_cross_condition_state_summary"; resources = ResourceLog(stage, config)
    before = {key: sha256(paths[key]) for key in ["representation", "compact"]}; started = time.perf_counter()
    with resources.operation("load_expression_free_representation_object"):
        obj = load_representation(config, paths)
        evidence = pd.concat([pd.read_csv(paths["output"] / "scgeo" / f"04_{x}_detailed_representation_state_evidence.csv") for x in ["tnf", "lps"]] + [pd.read_csv(paths["output"] / "scgeo" / f"04_{x}_coarse_representation_state_evidence.csv") for x in ["tnf", "lps"]], ignore_index=True, sort=False)
        consensus = pd.concat([pd.read_csv(paths["output"] / "scgeo" / f"04_{x}_primary_consensus.csv") for x in ["tnf", "lps"]], ignore_index=True)
        coarse_consensus = pd.concat([pd.read_csv(paths["output"] / "scgeo" / f"04_{x}_coarse_consensus.csv") for x in ["tnf", "lps"]], ignore_index=True)
        agreement = pd.concat([pd.read_csv(paths["output"] / "scgeo" / f"04_{x}_representation_agreement.csv") for x in ["tnf", "lps"]], ignore_index=True)
        loo = pd.concat([pd.read_csv(paths["output"] / "controls" / f"04_{x}_leave_one_representation.csv") for x in ["tnf", "lps"]], ignore_index=True)
    with resources.operation("frozen_scgeo_local_geometry_stability"):
        geometry = local_geometry(obj, config, sg)
        for key in ["representation_summary", "state_pair_summary", "state_graph_summary"]:
            geometry[key] = geometry[key].copy()
            geometry[key]["resampling_status"] = "descriptive_cell_resampling_only_where_applicable"
            geometry[key]["inference_status"] = "descriptive_only"
        atomic_csv(geometry["representation_summary"], paths["output"] / "scgeo" / "04_local_geometry_representation_summary.csv")
        atomic_csv(geometry["state_pair_summary"], paths["output"] / "scgeo" / "04_local_geometry_state_pair_summary.csv")
        atomic_csv(geometry["state_graph_summary"], paths["output"] / "scgeo" / "04_local_geometry_state_graph_summary.csv")
    with resources.operation("abundance_and_balanced_subsampling"):
        abundance, abundance_sens = abundance_tables(obj, config)
        balanced = balanced_subsampling(obj, config)
    with resources.operation("within_state_distribution_metrics"):
        distribution = distribution_metrics(obj, config)
    with resources.operation("shuffled_condition_computational_null"):
        shuffled = shuffled_controls(obj, config)
    with resources.operation("signature_displacement_context"):
        signatures, associations = signature_tables(obj, config, evidence)
    with resources.operation("write_full_evidence_and_controls"):
        # Attach state-specific local geometry median across PCA50-reference pairs.
        geom_state = geometry["state_pair_summary"]
        geom_state = geom_state[geom_state["metric"] == "neighbor_overlap"].groupby("state", observed=False)["median"].median().rename("local_geometry_preservation_median").reset_index()
        evidence = evidence.merge(geom_state, on="state", how="left")
        evidence = evidence.merge(consensus[["contrast", "state", "primary_consensus_status"]], on=["contrast", "state"], how="left")
        evidence = evidence.merge(coarse_consensus[["contrast", "state", "all_representation_consensus_status"]], on=["contrast", "state"], how="left")
        evidence["consensus_status"] = np.where(evidence["hierarchy"].eq("detailed"), evidence["primary_consensus_status"], evidence["all_representation_consensus_status"])
        evidence = evidence.merge(agreement[["contrast", "state", "normalized_effect_iqr_all_representations", "magnitude_rank_std_all_representations"]], on=["contrast", "state"], how="left")
        loo_max = loo.groupby(["contrast", "state"], observed=False)["relative_change"].max().rename("leave_one_representation_max_relative_change").reset_index()
        evidence = evidence.merge(loo_max, on=["contrast", "state"], how="left")
        abundance_wide = abundance[abundance["hierarchy"] == "detailed"].pivot(index="state", columns="condition", values="proportion").reset_index()
        abundance_long = []
        for contrast, treated in [("PBS_vs_TNF", "TNF"), ("PBS_vs_LPS", "LPS")]:
            part = abundance_wide[["state", "PBS", treated]].copy(); part["contrast"] = contrast
            part["abundance_proportion_PBS"] = part["PBS"]; part["abundance_proportion_treated"] = part[treated]
            part["abundance_proportion_change_treated_minus_PBS"] = part[treated] - part["PBS"]
            abundance_long.append(part[["contrast", "state", "abundance_proportion_PBS", "abundance_proportion_treated", "abundance_proportion_change_treated_minus_PBS"]])
        evidence = evidence.merge(pd.concat(abundance_long, ignore_index=True), on=["contrast", "state"], how="left")
        evidence = evidence.merge(distribution[["contrast", "state", "representation", "energy_distance_biased_full_value", "mmd_rbf_squared_biased_full_value", "sliced_wasserstein_full_value"]], on=["contrast", "state", "representation"], how="left")
        diffusion_exclusion, ambiguous_exclusion = exclusion_sensitivities(evidence, config)
        atomic_csv(evidence, paths["output"] / "scgeo" / "04_full_state_evidence.csv")
        atomic_csv(abundance, paths["output"] / "scgeo" / "04_descriptive_abundance_proportions.csv")
        atomic_csv(distribution, paths["output"] / "scgeo" / "04_within_state_distribution_metrics.csv")
        atomic_csv(signatures, paths["output"] / "scgeo" / "04_state_signature_changes.csv")
        atomic_csv(associations, paths["output"] / "scgeo" / "04_signature_displacement_associations.csv")
        atomic_csv(balanced, paths["output"] / "controls" / "04_balanced_state_displacement_subsampling.csv")
        atomic_csv(abundance_sens, paths["output"] / "controls" / "04_balanced_abundance_subsampling.csv")
        atomic_csv(shuffled, paths["output"] / "controls" / "04_shuffled_condition_computational_null.csv")
        atomic_csv(diffusion_exclusion, paths["output"] / "controls" / "04_exclude_diffusion_sensitivity.csv")
        atomic_csv(ambiguous_exclusion, paths["output"] / "controls" / "04_exclude_ambiguous_hspc_sensitivity.csv")
        plot_outputs(paths["output"], evidence, consensus, abundance, distribution, associations)
        atomic_csv(evidence, paths["output"] / "figure_sources" / "04_full_state_evidence.csv")
        atomic_csv(distribution, paths["output"] / "figure_sources" / "04_distribution_metrics.csv")
        atomic_csv(associations, paths["output"] / "figure_sources" / "04_signature_displacement_associations.csv")
        alt = (
            "Phase 3B deterministic alt text. PBS-to-TNF and PBS-to-LPS panels show frozen-ScGeo normalized state displacement "
            "for PCA30, PCA50, and scVI. A heatmap adds nested PCA20 and exploratory diffusion map without counting PCA views "
            "as independent confirmation. Comparison plots keep abundance, directional displacement, and within-state distribution "
            "statistics separate. Signature associations are descriptive across seven marker-inferred states. No biological-replicate "
            "p-values are displayed; UMAP is not used for quantitative geometry."
        )
        (paths["output"] / "alt_text" / "04_scgeo_treatment_geometry.txt").write_text(alt + "\n", encoding="utf-8")
    after = {key: sha256(paths[key]) for key in ["representation", "compact"]}
    if after != before: raise RuntimeError("An immutable input changed")
    summary = {
        "status": "passed", "inference_status": "descriptive_only", "biological_replicates": False,
        "inferential_p_values": False, "comparators_run": False, "umap_used_quantitatively": False,
        "n_full_evidence_rows": int(evidence.shape[0]), "n_distribution_rows": int(distribution.shape[0]),
        "n_balanced_control_rows": int(balanced.shape[0]), "n_shuffled_control_rows": int(shuffled.shape[0]),
        "runtime_seconds": time.perf_counter() - started, "input_sha256_before": before, "input_sha256_after": after,
        "peak_cpu_rss_gib": resources.peak_rss_gib,
    }
    atomic_json(summary, paths["output"] / "audit" / "04c_cross_condition_summary.json")
    write_resources(resources, paths["output"]); write_metadata(stage, config, paths, resources, summary, sg)
    del obj; close_and_collect(); return summary


def write_metadata(stage: str, config: dict[str, Any], paths: dict[str, Path], resources: ResourceLog, summary: dict[str, Any], sg: Any) -> None:
    metadata = {
        **summary, "stage": stage, "timestamp_utc": utc_now(), "inference_status": "descriptive_only",
        "representation_h5ad": str(paths["representation"]), "representation_sha256": sha256(paths["representation"]),
        "compact_h5ad": str(paths["compact"]), "compact_sha256": sha256(paths["compact"]),
        "notebook_repository_commit": git_commit(ROOT), "frozen_scgeo_commit": git_commit(paths["source_repo"]),
        "scgeo_modified": False, "sample_key": None, "pseudo_replicates": False,
        "forbidden_methods_run": [], "package_versions": package_versions(["scgeo"]),
        "scgeo_defaults_observed": scgeo_defaults_snapshot(sg),
        "peak_cpu_rss_gib": resources.peak_rss_gib, "peak_gpu_allocated_gib": resources.peak_gpu_allocated_gib,
    }
    atomic_json(metadata, paths["output"] / "metadata" / f"{stage}_metadata.json")
    atomic_json({"timestamp_utc": utc_now(), "python": platform.python_version(), "packages": package_versions(["scgeo"]), "inference_status": "descriptive_only"}, paths["output"] / "version_records" / f"{stage}_versions.json")
