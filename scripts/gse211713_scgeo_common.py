"""Coverage-qualified, replicate-aware ScGeo workflow for GSE211713 C7 v2."""

from __future__ import annotations

import gc
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist

from scripts.gse211713_representation_common import (
    ROOT,
    ResourceMonitor,
    atomic_csv,
    atomic_json,
    configured_paths as representation_paths,
    git_commit,
    git_status,
    sha256,
    utc_now,
)


CONFIG_PATH = ROOT / "configs/gse211713_scgeo_v1.json"
SOURCE_REPO = Path(os.environ.get("SCGEO_SOURCE_REPO", ROOT.parent / "scgeo")).resolve()
if str(SOURCE_REPO) not in sys.path:
    sys.path.insert(0, str(SOURCE_REPO))

from scgeo.tl._robust_shift import _estimate_center, _robust_scale  # noqa: E402


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def configured_paths(config: dict[str, Any] | None = None) -> dict[str, Path]:
    config = load_config() if config is None else config
    output = Path(os.environ.get("SCGEO_GSE211713_C7_V2_OUTPUT_DIR", config["default_output_dir"]))
    if not output.is_absolute():
        output = ROOT / output
    return {
        "representation": Path(os.environ.get("SCGEO_GSE211713_REPRESENTATION_V2_H5AD", config["default_representation_h5ad"])).resolve(),
        "compact": Path(config["default_compact_h5ad"]).resolve(),
        "annotated": Path(config["default_annotated_h5ad"]).resolve(),
        "output": output.resolve(),
        "source_repo": SOURCE_REPO,
    }


def ensure_tree(output: Path) -> None:
    for directory in ["per_contrast", "metadata", "execution", "figures", "figure_sources", "alt_text"]:
        (output / directory).mkdir(parents=True, exist_ok=True)


def validate_environment(config: dict[str, Any], paths: dict[str, Path]) -> None:
    branch = subprocess.check_output(["git", "-C", str(ROOT), "branch", "--show-current"], text=True).strip()
    if branch != config["required_branch"]:
        raise RuntimeError(f"Branch mismatch: {branch}")
    if git_commit(paths["source_repo"]) != config["frozen_scgeo_commit"] or git_status(paths["source_repo"]):
        raise RuntimeError("Frozen ScGeo commit or cleanliness check failed")
    if sha256(paths["compact"]) != config["expected_compact_sha256"]:
        raise RuntimeError("Compact input checksum mismatch")
    if sha256(paths["annotated"]) != config["expected_annotated_sha256"]:
        raise RuntimeError("Annotated input checksum mismatch")
    c6_report_path = ROOT / config["representation_checksum_report"]
    if not c6_report_path.is_file():
        raise RuntimeError("C6 representation checksum report is absent")
    c6_report = json.loads(c6_report_path.read_text())
    if c6_report.get("status") != "passed" or sha256(paths["representation"]) != c6_report.get("representation_sha256"):
        raise RuntimeError("Representation v2 checksum validation failed")
    for path in config["coverage_tables"].values():
        if not (ROOT / path).is_file():
            raise RuntimeError(f"Frozen coverage table missing: {path}")


def load_representation(paths: dict[str, Path], config: dict[str, Any]) -> ad.AnnData:
    rep = ad.read_h5ad(paths["representation"])
    required = set(config["representations"]["all"])
    if required.difference(rep.obsm.keys()):
        raise RuntimeError("Representation v2 lacks a quantitative embedding")
    if "X_umap_display_only" not in rep.obsm:
        raise RuntimeError("Display-only UMAP is absent")
    if rep.n_vars != 0:
        raise RuntimeError("Representation object must remain zero-gene")
    return rep


def contrast_spec(config: dict[str, Any], contrast_id: str) -> dict[str, Any]:
    for item in config["primary_contrasts"] + config["secondary_contrasts"]:
        if item["contrast_id"] == contrast_id:
            return dict(item)
    raise KeyError(contrast_id)


def group_mask(obs: pd.DataFrame, label: str) -> np.ndarray:
    if label == "control":
        return obs["dose_gy"].eq(0).to_numpy()
    if label == "17Gy_early":
        return (obs["dose_gy"].eq(17) & obs["month_post_irradiation"].isin([1, 2])).to_numpy()
    if label == "17Gy_late":
        return (obs["dose_gy"].eq(17) & obs["month_post_irradiation"].isin([4, 5])).to_numpy()
    if label == "17Gy_all":
        return obs["dose_gy"].eq(17).to_numpy()
    if label == "10Gy_all":
        return obs["dose_gy"].eq(10).to_numpy()
    raise KeyError(label)


def _coverage_contrast_name(contrast_id: str) -> str:
    return {
        "control_vs_17gy_early": "control_vs_17Gy_early",
        "control_vs_17gy_late": "control_vs_17Gy_late",
        "17gy_early_vs_late": "17Gy_early_vs_late",
        "control_vs_all17gy": "control_vs_all_17Gy",
        "control_vs_all10gy": "control_vs_all_10Gy",
    }[contrast_id]


def eligibility_table(rep: ad.AnnData, config: dict[str, Any], contrast_id: str, hierarchy: str) -> pd.DataFrame:
    spec = contrast_spec(config, contrast_id)
    state_column = "state_major" if hierarchy == "major" else "fibroblast_state"
    states = config["states"][hierarchy]
    threshold = int(config["coverage"]["major_min_cells_per_mouse"] if hierarchy == "major" else config["coverage"]["fibroblast_min_cells_per_mouse"])
    minimum_mice = int(config["coverage"]["minimum_eligible_mice_per_group"])
    frozen = pd.read_csv(ROOT / config["coverage_tables"]["contrast_eligibility"])
    source_state = {state: state.lower() for state in states} if hierarchy == "major" else config["source_fibroblast_labels"]
    rows: list[dict[str, Any]] = []
    for group in [spec["group0"], spec["group1"]]:
        mask = group_mask(rep.obs, group)
        group_mice = sorted(rep.obs.loc[mask, "mouse_id"].astype(str).unique())
        for state in states:
            counts = rep.obs.loc[mask & rep.obs[state_column].astype(str).eq(state).to_numpy(), "mouse_id"].astype(str).value_counts().reindex(group_mice, fill_value=0)
            for mouse in group_mice:
                count = int(counts.loc[mouse])
                rows.append({
                    "contrast": contrast_id, "hierarchy": hierarchy, "state": state, "group": group, "mouse_id": mouse,
                    "state_cells": count, "minimum_cells_per_mouse": threshold, "eligible": bool(count >= threshold),
                    "exclusion_reason": "" if count >= threshold else f"fewer_than_{threshold}_cells",
                    "inference_status": spec["inference_status"], "coverage_source": "frozen_C1_C5_tables_plus_exact_obs_crosscheck",
                })
            expected = frozen.loc[
                frozen["contrast"].eq(_coverage_contrast_name(contrast_id))
                & frozen["hierarchy"].eq("major" if hierarchy == "major" else "refined")
                & frozen["state"].astype(str).eq(source_state[state])
            ]
            if expected.shape[0] != 1:
                raise RuntimeError(f"Frozen coverage row not uniquely recovered: {contrast_id}, {hierarchy}, {state}")
            expected_row = expected.iloc[0]
            if int(expected_row["minimum_eligible_mice_per_group"]) != minimum_mice:
                raise RuntimeError("Frozen minimum-mouse criterion differs from approved config")
            if hierarchy == "major" and int(expected_row["minimum_cells_per_mouse"]) != threshold:
                raise RuntimeError("Frozen major-state cell criterion differs from approved config")
            expected_group_counts = (
                int(expected_row["eligible_mice_group_a"]),
                int(expected_row["eligible_mice_group_b"]),
            ) if hierarchy == "major" else (
                int(expected_row["fibroblast_eligible_mice_group_a_at20"]),
                int(expected_row["fibroblast_eligible_mice_group_b_at20"]),
            )
            observed_group_count = int((counts >= threshold).sum())
            expected_group_count = expected_group_counts[0 if group == spec["group0"] else 1]
            if observed_group_count != expected_group_count:
                raise RuntimeError(
                    f"Frozen coverage count mismatch: {contrast_id}, {state}, {group}; "
                    f"observed={observed_group_count}, frozen={expected_group_count}"
                )
    frame = pd.DataFrame(rows)
    eligible_counts = frame.groupby(["state", "group"], observed=True)["eligible"].sum()
    for state in states:
        for group in [spec["group0"], spec["group1"]]:
            if int(eligible_counts.loc[(state, group)]) < minimum_mice:
                raise RuntimeError(f"Prespecified state lacks eligible mice: {contrast_id}, {state}, {group}")
    return frame


def qualified_subset(rep: ad.AnnData, eligibility: pd.DataFrame, config: dict[str, Any], contrast_id: str, hierarchy: str) -> ad.AnnData:
    spec = contrast_spec(config, contrast_id)
    state_column = "state_major" if hierarchy == "major" else "fibroblast_state"
    selected = np.zeros(rep.n_obs, dtype=bool)
    analysis_state = np.full(rep.n_obs, "", dtype=object)
    for state in config["states"][hierarchy]:
        allowed = set(eligibility.loc[eligibility["state"].eq(state) & eligibility["eligible"], "mouse_id"].astype(str))
        mask = rep.obs[state_column].astype(str).eq(state).to_numpy() & rep.obs["mouse_id"].astype(str).isin(allowed).to_numpy()
        selected |= mask
        analysis_state[mask] = state
    subset = rep[selected].copy()
    subset.obs["analysis_state"] = analysis_state[selected]
    m0 = group_mask(subset.obs, spec["group0"])
    m1 = group_mask(subset.obs, spec["group1"])
    if np.any(m0 & m1) or not np.all(m0 | m1):
        raise RuntimeError("Qualified contrast subset contains invalid group membership")
    subset.obs[config["condition_key"]] = np.where(m0, spec["group0"], spec["group1"])
    return subset


def _stable_seed(base: int, *parts: str) -> int:
    value = hashlib.sha256((str(base) + "|" + "|".join(parts)).encode()).digest()
    return int.from_bytes(value[:4], "little") & 0x7FFFFFFF


def mouse_centers(subset: ad.AnnData, state: str, representation: str, config: dict[str, Any], spec: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, list[str], list[str], dict[str, int]]:
    coords = np.asarray(subset.obsm[representation], dtype=np.float64)
    state_mask = subset.obs["analysis_state"].astype(str).eq(state).to_numpy()
    samples = subset.obs[config["sample_key"]].astype(str).to_numpy()
    counts: dict[str, int] = {}
    outputs: list[tuple[np.ndarray, list[str]]] = []
    for group in [spec["group0"], spec["group1"]]:
        group_state = state_mask & subset.obs[config["condition_key"]].eq(group).to_numpy()
        ids = list(dict.fromkeys(samples[group_state].tolist()))
        centers = []
        for mouse in ids:
            mask = group_state & (samples == mouse)
            counts[mouse] = int(mask.sum())
            centers.append(_estimate_center(coords[mask], "geometric_median", float(config["estimator"]["trim_fraction"])))
        outputs.append((np.vstack(centers), ids))
    return outputs[0][0], outputs[1][0], outputs[0][1], outputs[1][1], counts


def geometry_and_bootstrap(a: np.ndarray, b: np.ndarray, config: dict[str, Any], seed: int) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray]:
    center = "geometric_median"
    trim = float(config["estimator"]["trim_fraction"])
    c0 = _estimate_center(a, center, trim)
    c1 = _estimate_center(b, center, trim)
    delta = c1 - c0
    raw = float(np.linalg.norm(delta))
    scale, scale_value = _robust_scale(a, b, c0, c1, config["estimator"]["normalize_by"])
    normalized = float(raw / scale) if np.isfinite(scale) and scale > 0 else np.nan
    rng = np.random.RandomState(seed)
    n_boot = int(config["estimator"]["n_boot"])
    boot_delta = np.empty((n_boot, a.shape[1]), dtype=np.float64)
    for index in range(n_boot):
        aa = a[rng.choice(a.shape[0], size=a.shape[0], replace=True)]
        bb = b[rng.choice(b.shape[0], size=b.shape[0], replace=True)]
        boot_delta[index] = _estimate_center(bb, center, trim) - _estimate_center(aa, center, trim)
    boot_raw = np.linalg.norm(boot_delta, axis=1)
    boot_normalized = boot_raw / scale
    nonzero = np.isfinite(boot_raw) & (boot_raw > 1e-12)
    resultant = float(np.linalg.norm((boot_delta[nonzero] / boot_raw[nonzero, None]).mean(axis=0))) if np.any(nonzero) else np.nan
    stability = float(np.mean((boot_delta[nonzero] @ delta) >= 0)) if np.any(nonzero) and raw > 1e-12 else np.nan
    result = {
        "delta_norm": raw, "normalized_delta_norm": normalized, "normalization_scale": float(scale_value),
        "raw_ci95_low": float(np.percentile(boot_raw, 2.5)), "raw_ci95_high": float(np.percentile(boot_raw, 97.5)),
        "normalized_ci95_low": float(np.percentile(boot_normalized, 2.5)), "normalized_ci95_high": float(np.percentile(boot_normalized, 97.5)),
        "directional_resultant_length": resultant, "direction_stability": stability, "sign_stability": stability,
        "bootstrap_iterations": n_boot, "bootstrap_seed": seed, "resampling_status": config["estimator"]["resampling_label"],
    }
    return result, c0, c1, boot_raw


def normalized_statistic(a: np.ndarray, b: np.ndarray, config: dict[str, Any]) -> float:
    center, trim = "geometric_median", float(config["estimator"]["trim_fraction"])
    c0, c1 = _estimate_center(a, center, trim), _estimate_center(b, center, trim)
    raw = float(np.linalg.norm(c1 - c0))
    scale, _ = _robust_scale(a, b, c0, c1, config["estimator"]["normalize_by"])
    return float(raw / scale) if np.isfinite(scale) and scale > 0 else np.nan


def exact_permutation(a: np.ndarray, b: np.ndarray, config: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray]:
    observed = normalized_statistic(a, b, config)
    pooled = np.vstack([a, b])
    stats = []
    for chosen in combinations(range(pooled.shape[0]), a.shape[0]):
        mask = np.zeros(pooled.shape[0], dtype=bool)
        mask[list(chosen)] = True
        stats.append(normalized_statistic(pooled[mask], pooled[~mask], config))
    null = np.asarray(stats, dtype=float)
    tail = int(np.sum(null >= observed))
    assignments = int(null.size)
    equal_groups = bool(a.shape[0] == b.shape[0])
    return {
        "observed_normalized_displacement": observed, "assignment_count": assignments,
        "tail_count_greater_equal": tail, "empirical_upper_tail_fraction": float(tail / assignments),
        "minimum_attainable_grid": float(1 / assignments),
        "effective_minimum_with_complement_symmetry": float(2 / assignments) if equal_groups else float(1 / assignments),
        "ties": "included_statistic_greater_equal_observed", "observed_assignment_included": True,
        "plus_one_correction": False, "equal_group_complement_symmetry": equal_groups,
    }, null


def multivariate_energy(a: np.ndarray, b: np.ndarray) -> float:
    return float(2 * cdist(a, b).mean() - cdist(a, a).mean() - cdist(b, b).mean())


def mmd_rbf(a: np.ndarray, b: np.ndarray) -> float:
    pooled = np.vstack([a, b])
    distances = pdist(pooled)
    sigma = float(np.median(distances[distances > 0])) if np.any(distances > 0) else 1.0
    gamma = 1 / (2 * sigma**2)
    return float(np.exp(-gamma * cdist(a, a, "sqeuclidean")).mean() + np.exp(-gamma * cdist(b, b, "sqeuclidean")).mean() - 2 * np.exp(-gamma * cdist(a, b, "sqeuclidean")).mean())


def sliced_wasserstein(a: np.ndarray, b: np.ndarray, projections: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(projections):
        direction = rng.normal(size=a.shape[1])
        direction /= np.linalg.norm(direction)
        pa, pb = np.sort(a @ direction), np.sort(b @ direction)
        quantiles = np.linspace(0, 1, max(pa.size, pb.size))
        ia = np.interp(quantiles, np.linspace(0, 1, pa.size), pa)
        ib = np.interp(quantiles, np.linspace(0, 1, pb.size), pb)
        values.append(np.mean(np.abs(ia - ib)))
    return float(np.mean(values))


def bh_adjust(values: pd.Series) -> pd.Series:
    array = values.to_numpy(dtype=float)
    order = np.argsort(array)
    ranked = array[order]
    adjusted = np.minimum.accumulate((ranked * len(array) / np.arange(1, len(array) + 1))[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return pd.Series(output, index=values.index)


def run_contrast(contrast_id: str, hierarchy: str) -> dict[str, Any]:
    import scgeo as sg
    config, paths = load_config(), configured_paths()
    ensure_tree(paths["output"])
    validate_environment(config, paths)
    spec = contrast_spec(config, contrast_id)
    started = time.perf_counter()
    monitor = ResourceMonitor(f"{contrast_id}_{hierarchy}", config)
    rep = load_representation(paths, config)
    eligibility = eligibility_table(rep, config, contrast_id, hierarchy)
    subset = qualified_subset(rep, eligibility, config, contrast_id, hierarchy)
    states = config["states"][hierarchy]
    with monitor.operation("frozen_scgeo_primary_consensus"):
        primary = sg.tl.representation_stability(
            subset, reps=config["representations"]["primary"], node_key="analysis_state", condition_key=config["condition_key"],
            group0=spec["group0"], group1=spec["group1"], sample_key=config["sample_key"], center="geometric_median",
            trim_fraction=float(config["estimator"]["trim_fraction"]), n_boot=int(config["estimator"]["n_boot"]), min_cells=1,
            seed=int(config["random_seed"]), store_key=f"{contrast_id}_{hierarchy}_primary",
        )
    with monitor.operation("frozen_scgeo_all_representation_sensitivity"):
        all_reps = sg.tl.representation_stability(
            subset, reps=config["representations"]["all"], node_key="analysis_state", condition_key=config["condition_key"],
            group0=spec["group0"], group1=spec["group1"], sample_key=config["sample_key"], center="geometric_median",
            trim_fraction=float(config["estimator"]["trim_fraction"]), n_boot=int(config["estimator"]["n_boot"]), min_cells=1,
            seed=int(config["random_seed"]), store_key=f"{contrast_id}_{hierarchy}_all",
        )
    consensus = primary["consensus_state"].copy().rename(columns={"node": "state"})
    consensus["contrast"] = contrast_id
    consensus["hierarchy"] = hierarchy
    consensus["inference_status"] = spec["inference_status"]
    consensus["primary_representations"] = json.dumps(config["representations"]["primary"])
    all_per = all_reps["per_rep_state"].copy().rename(columns={"node": "state", "rep": "representation"})
    primary_consensus_columns = ["state", "consensus_label", "pairwise_spearman_median", "magnitude_rank_std", "loo_rep_magnitude_max_relative_deviation"]
    all_per = all_per.merge(consensus[primary_consensus_columns], on="state", how="left", suffixes=("", "_primary_consensus"))
    evidence_rows: list[dict[str, Any]] = []
    center_rows: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    bootstrap_draw_rows: list[dict[str, Any]] = []
    permutation_rows: list[dict[str, Any]] = []
    permutation_value_rows: list[dict[str, Any]] = []
    distribution_rows: list[dict[str, Any]] = []
    abundance_rows: list[dict[str, Any]] = []
    with monitor.operation("mouse_centers_bootstrap_permutation_distribution"):
        for state in states:
            state_elig = eligibility[eligibility["state"].eq(state)]
            eligible0 = state_elig.loc[state_elig["group"].eq(spec["group0"]) & state_elig["eligible"], "mouse_id"].astype(str).tolist()
            eligible1 = state_elig.loc[state_elig["group"].eq(spec["group1"]) & state_elig["eligible"], "mouse_id"].astype(str).tolist()
            excluded0 = state_elig.loc[state_elig["group"].eq(spec["group0"]) & ~state_elig["eligible"], "mouse_id"].astype(str).tolist()
            excluded1 = state_elig.loc[state_elig["group"].eq(spec["group1"]) & ~state_elig["eligible"], "mouse_id"].astype(str).tolist()
            for representation in config["representations"]["all"]:
                a, b, ids0, ids1, counts = mouse_centers(subset, state, representation, config, spec)
                seed = _stable_seed(int(config["random_seed"]), contrast_id, hierarchy, state, representation, "bootstrap")
                geometry, group_center0, group_center1, boot_raw = geometry_and_bootstrap(a, b, config, seed)
                frozen_row = all_per.loc[all_per["state"].astype(str).eq(state) & all_per["representation"].eq(representation)]
                if frozen_row.shape[0] != 1:
                    raise RuntimeError("Frozen ScGeo row was not uniquely recovered")
                frozen_row = frozen_row.iloc[0]
                if not np.isclose(float(frozen_row["normalized_delta_norm"]), geometry["normalized_delta_norm"], rtol=1e-6, atol=1e-7):
                    raise RuntimeError(f"Reconstructed normalized effect differs from frozen ScGeo: {contrast_id}, {state}, {representation}")
                role = "primary" if representation in config["representations"]["primary"] else "dimensional_sensitivity" if representation in config["representations"]["dimensional_sensitivity"] else "exploratory_sensitivity"
                row = {
                    "contrast": contrast_id, "hierarchy": hierarchy, "state": state, "representation": representation,
                    "representation_role": role, "inference_status": spec["inference_status"],
                    "eligible_mouse_ids_group0": json.dumps(ids0), "eligible_mouse_ids_group1": json.dumps(ids1),
                    "excluded_mouse_ids_group0": json.dumps(excluded0), "excluded_mouse_ids_group1": json.dumps(excluded1),
                    "n_mice_group0": len(ids0), "n_mice_group1": len(ids1),
                    "cells_per_mouse_group0": json.dumps({mouse: counts[mouse] for mouse in ids0}, sort_keys=True),
                    "cells_per_mouse_group1": json.dumps({mouse: counts[mouse] for mouse in ids1}, sort_keys=True),
                    "consensus_label": frozen_row["consensus_label"],
                    "primary_pairwise_spearman_median": frozen_row["pairwise_spearman_median"],
                    "primary_magnitude_rank_std": frozen_row["magnitude_rank_std"],
                    "leave_one_primary_max_relative_deviation": frozen_row["loo_rep_magnitude_max_relative_deviation"],
                    "coverage_status": "qualified_per_mouse_before_scgeo", **geometry,
                }
                evidence_rows.append(row)
                bootstrap_rows.append({key: row[key] for key in [
                    "contrast", "hierarchy", "state", "representation", "inference_status", "n_mice_group0", "n_mice_group1",
                    "delta_norm", "normalized_delta_norm", "normalization_scale", "raw_ci95_low", "raw_ci95_high",
                    "normalized_ci95_low", "normalized_ci95_high", "directional_resultant_length", "direction_stability", "sign_stability",
                    "bootstrap_iterations", "bootstrap_seed", "resampling_status",
                ]})
                for draw, value in enumerate(boot_raw, start=1):
                    bootstrap_draw_rows.append({"contrast": contrast_id, "hierarchy": hierarchy, "state": state, "representation": representation, "draw": draw, "raw_displacement": float(value), "normalized_displacement": float(value / geometry["normalization_scale"])})
                for group, ids, matrix in [(spec["group0"], ids0, a), (spec["group1"], ids1, b)]:
                    for mouse, vector in zip(ids, matrix, strict=True):
                        center_row = {"contrast": contrast_id, "hierarchy": hierarchy, "state": state, "representation": representation, "group": group, "mouse_id": mouse, "state_cells": counts[mouse], "dimensions": int(vector.size)}
                        center_row.update({f"coordinate_{i + 1}": float(value) for i, value in enumerate(vector)})
                        center_rows.append(center_row)
                if spec["inference_status"] == "replicate_aware_primary" and representation in config["representations"]["primary"]:
                    permutation_summary, permutation_values = exact_permutation(a, b, config)
                    permutation_rows.append({"contrast": contrast_id, "hierarchy": hierarchy, "state": state, "representation": representation, "n_mice_group0": len(ids0), "n_mice_group1": len(ids1), **permutation_summary})
                    for assignment_index, statistic in enumerate(permutation_values, start=1):
                        permutation_value_rows.append({
                            "contrast": contrast_id, "hierarchy": hierarchy, "state": state,
                            "representation": representation, "assignment_index": assignment_index,
                            "normalized_displacement": float(statistic),
                            "is_observed_assignment": bool(assignment_index == 1),
                        })
                cap = int(config["distribution"]["max_cells_per_state_group"])
                state_mask = subset.obs["analysis_state"].astype(str).eq(state).to_numpy()
                coords = np.asarray(subset.obsm[representation], dtype=np.float32)
                rng = np.random.default_rng(_stable_seed(int(config["random_seed"]), contrast_id, state, representation, "distribution"))
                pooled = []
                for group in [spec["group0"], spec["group1"]]:
                    mask = state_mask & subset.obs[config["condition_key"]].eq(group).to_numpy()
                    values = coords[mask]
                    if values.shape[0] > cap:
                        values = values[np.sort(rng.choice(values.shape[0], size=cap, replace=False))]
                    pooled.append(values)
                for unit, aa, bb in [("pooled_cells_descriptive", pooled[0], pooled[1]), ("biological_mouse_centers_descriptive", a, b)]:
                    distribution_rows.append({
                        "contrast": contrast_id, "hierarchy": hierarchy, "state": state, "representation": representation, "comparison_unit": unit,
                        "n_group0": int(aa.shape[0]), "n_group1": int(bb.shape[0]), "energy_distance": multivariate_energy(aa, bb),
                        "mmd_rbf": mmd_rbf(aa, bb), "sliced_wasserstein": sliced_wasserstein(aa, bb, int(config["distribution"]["sliced_wasserstein_projections"]), seed),
                        "inference_status": "descriptive_distribution_metric_not_replicate_test",
                    })
            state_col = "state_major" if hierarchy == "major" else "fibroblast_state"
            for group in [spec["group0"], spec["group1"]]:
                mask = group_mask(rep.obs, group)
                group_obs = rep.obs.loc[mask]
                totals = group_obs["mouse_id"].astype(str).value_counts()
                counts_by_mouse = group_obs.loc[group_obs[state_col].astype(str).eq(state), "mouse_id"].astype(str).value_counts()
                eligible_ids = set(state_elig.loc[state_elig["group"].eq(group) & state_elig["eligible"], "mouse_id"].astype(str))
                for mouse, total in totals.items():
                    count = int(counts_by_mouse.get(mouse, 0))
                    abundance_rows.append({"contrast": contrast_id, "hierarchy": hierarchy, "state": state, "group": group, "mouse_id": mouse, "state_cells": count, "total_retained_cells": int(total), "state_fraction": float(count / total), "coverage_eligible": mouse in eligible_ids, "comparison_unit": "biological_mouse", "inference_status": spec["inference_status"]})
    evidence = pd.DataFrame(evidence_rows)
    centers = pd.DataFrame(center_rows)
    bootstrap = pd.DataFrame(bootstrap_rows)
    draws = pd.DataFrame(bootstrap_draw_rows)
    permutation_columns = [
        "contrast", "hierarchy", "state", "representation", "n_mice_group0", "n_mice_group1",
        "observed_normalized_displacement", "assignment_count", "tail_count_greater_equal",
        "empirical_upper_tail_fraction", "minimum_attainable_grid",
        "effective_minimum_with_complement_symmetry", "ties", "observed_assignment_included",
        "plus_one_correction", "equal_group_complement_symmetry",
    ]
    permutation_value_columns = [
        "contrast", "hierarchy", "state", "representation", "assignment_index",
        "normalized_displacement", "is_observed_assignment",
    ]
    permutations = pd.DataFrame(permutation_rows, columns=permutation_columns)
    permutation_values = pd.DataFrame(permutation_value_rows, columns=permutation_value_columns)
    distributions = pd.DataFrame(distribution_rows)
    abundance = pd.DataFrame(abundance_rows)
    key = ["contrast", "hierarchy", "state", "representation"]
    if evidence.duplicated(key).any():
        raise RuntimeError("Duplicate per-contrast evidence keys")
    directory = paths["output"] / "per_contrast" / f"{contrast_id}__{hierarchy}"
    directory.mkdir(parents=True, exist_ok=False)
    atomic_csv(evidence, directory / "state_evidence.csv")
    atomic_csv(consensus, directory / "consensus.csv")
    atomic_csv(centers, directory / "individual_mouse_centers.csv")
    atomic_csv(bootstrap, directory / "bootstrap_intervals.csv")
    atomic_csv(draws, directory / "bootstrap_draws.csv")
    atomic_csv(permutations, directory / "exact_permutations.csv")
    atomic_csv(permutation_values, directory / "exact_permutation_values.csv")
    atomic_csv(distributions, directory / "distribution_metrics.csv")
    atomic_csv(abundance, directory / "mouse_level_abundance.csv")
    atomic_csv(eligibility, directory / "eligibility.csv")
    report = {
        "status": "passed", "contrast": contrast_id, "hierarchy": hierarchy, "inference_status": spec["inference_status"],
        "states": states, "evidence_rows": int(evidence.shape[0]), "consensus_rows": int(consensus.shape[0]),
        "mouse_center_rows": int(centers.shape[0]), "permutation_rows": int(permutations.shape[0]),
        "output_checksums": {path.name: sha256(path) for path in sorted(directory.iterdir()) if path.is_file()},
        **monitor.summary(time.perf_counter() - started),
    }
    atomic_csv(pd.DataFrame(monitor.rows), directory / "resource_log.csv")
    atomic_json(report, directory / "report.json")
    del rep, subset
    gc.collect()
    return report


def assemble_canonical_outputs() -> dict[str, Any]:
    config, paths = load_config(), configured_paths()
    validate_environment(config, paths)
    manifests = config["canonical_manifest"]
    frames: dict[str, list[pd.DataFrame]] = {name: [] for name in ["state_evidence", "consensus", "individual_mouse_centers", "bootstrap_intervals", "bootstrap_draws", "exact_permutations", "exact_permutation_values", "distribution_metrics", "mouse_level_abundance", "eligibility"]}
    source_paths = []
    for item in manifests:
        directory = paths["output"] / "per_contrast" / f"{item['contrast_id']}__{item['hierarchy']}"
        if not directory.is_dir():
            raise RuntimeError(f"Manifest source directory missing: {directory}")
        for name in frames:
            path = directory / f"{name}.csv"
            if not path.is_file() or "summary" in path.name or path.parent == paths["output"]:
                raise RuntimeError(f"Invalid canonical source: {path}")
            frames[name].append(pd.read_csv(path))
            source_paths.append(path)
    combined = {name: pd.concat(parts, ignore_index=True, sort=False) for name, parts in frames.items()}
    evidence = combined["state_evidence"]
    key = ["contrast", "hierarchy", "state", "representation"]
    if evidence.duplicated(key).any():
        raise RuntimeError("Canonical full evidence contains duplicate keys")
    if evidence["state"].astype(str).str.lower().isin(["nan", "ambiguous", "mast", "smoke_test"]).any() or evidence["contrast"].astype(str).str.contains("smoke", case=False).any():
        raise RuntimeError("Forbidden state or smoke-test result entered canonical evidence")
    observed_major = set(evidence.loc[evidence["hierarchy"].eq("major"), "state"])
    observed_fibro = set(evidence.loc[evidence["hierarchy"].eq("fibroblast"), "state"])
    if observed_major != set(config["states"]["major"]) or observed_fibro != set(config["states"]["fibroblast"]):
        raise RuntimeError("Canonical state sets differ from approved sets")
    primary_ids = {item["contrast_id"] for item in config["primary_contrasts"]}
    major_primary = evidence[evidence["contrast"].isin(primary_ids) & evidence["hierarchy"].eq("major")].copy()
    fibro_primary = evidence[evidence["contrast"].isin(primary_ids) & evidence["hierarchy"].eq("fibroblast")].copy()
    consensus = combined["consensus"]
    primary_consensus = consensus[consensus["contrast"].isin(primary_ids)].copy()
    consensus_key = ["contrast", "hierarchy", "state"]
    if primary_consensus.duplicated(consensus_key).any() or primary_consensus.shape[0] != 27:
        raise RuntimeError("Primary consensus keys are not unique and complete")
    sensitivity = evidence[evidence["contrast"].isin(primary_ids)].copy()
    negative = primary_consensus[primary_consensus["consensus_label"].ne("stable_effect")].copy()
    permutations = combined["exact_permutations"]
    if not permutations.empty:
        expected = {"control_vs_17gy_early": 126, "control_vs_17gy_late": 126, "17gy_early_vs_late": 70}
        for contrast, count in expected.items():
            if set(permutations.loc[permutations["contrast"].eq(contrast), "assignment_count"].astype(int)) != {count}:
                raise RuntimeError(f"Exact permutation count mismatch: {contrast}")
    bh = permutations[permutations["hierarchy"].eq("major")].copy()
    bh["bh_adjusted_fraction_within_contrast_representation_major_states"] = bh.groupby(["contrast", "representation"], observed=True)["empirical_upper_tail_fraction"].transform(bh_adjust)
    consistency_rows = []
    for (contrast, hierarchy, state), group in permutations.groupby(["contrast", "hierarchy", "state"], observed=True):
        by_rep = group.set_index("representation")["empirical_upper_tail_fraction"]
        consistency_rows.append({"contrast": contrast, "hierarchy": hierarchy, "state": state, "representations": json.dumps(sorted(by_rep.index)), "minimum_raw_fraction": float(by_rep.min()), "maximum_raw_fraction": float(by_rep.max()), "all_primary_fractions_below_0_05": bool((by_rep < 0.05).all()), "no_pvalue_averaging_or_combination": True})
    consistency = pd.DataFrame(consistency_rows)
    outputs = {
        "full_state_evidence.csv": evidence,
        "major_primary_state_evidence.csv": major_primary,
        "fibroblast_primary_state_evidence.csv": fibro_primary,
        "primary_consensus_state_evidence.csv": primary_consensus,
        "all_representation_sensitivity.csv": sensitivity,
        "negative_neutral_or_unstable_states.csv": negative,
        "individual_mouse_centers.csv": combined["individual_mouse_centers"],
        "biological_mouse_bootstrap_intervals.csv": combined["bootstrap_intervals"],
        "biological_mouse_bootstrap_draws.csv": combined["bootstrap_draws"],
        "exact_mouse_permutations.csv": permutations,
        "exact_mouse_permutation_values.csv": combined["exact_permutation_values"],
        "exact_mouse_permutation_bh_sensitivity.csv": bh,
        "exact_mouse_permutation_representation_consistency.csv": consistency,
        "mouse_level_abundance.csv": combined["mouse_level_abundance"],
        "distribution_metrics.csv": combined["distribution_metrics"],
        "eligibility_and_coverage_used.csv": combined["eligibility"],
    }
    for filename, frame in outputs.items():
        atomic_csv(frame, paths["output"] / filename)
    schema = {
        "schema_version": "gse211713_scgeo_v2_output_schema",
        "canonical_sources": [str(path.relative_to(paths["output"])) for path in source_paths],
        "files": {
            filename: {
                "rows": int(frame.shape[0]), "columns": frame.columns.tolist(),
                "key": key if filename in ["full_state_evidence.csv", "major_primary_state_evidence.csv", "fibroblast_primary_state_evidence.csv", "all_representation_sensitivity.csv"] else consensus_key if filename in ["primary_consensus_state_evidence.csv", "negative_neutral_or_unstable_states.csv"] else None,
            }
            for filename, frame in outputs.items()
        },
        "duplicate_full_evidence_keys": int(evidence.duplicated(key).sum()),
        "duplicate_primary_consensus_keys": int(primary_consensus.duplicated(consensus_key).sum()),
        "smoke_test_present": False, "nan_state_present": False,
    }
    atomic_json(schema, paths["output"] / "output_schema.json")
    report = {
        "status": "passed", "canonical_rows": {name: int(frame.shape[0]) for name, frame in outputs.items()},
        "output_checksums": {path.name: sha256(path) for path in sorted(paths["output"].iterdir()) if path.is_file()},
        "representation_sha256": sha256(paths["representation"]), "compact_sha256": sha256(paths["compact"]), "annotated_sha256": sha256(paths["annotated"]),
        "timestamp_utc": utc_now(),
    }
    atomic_json(report, paths["output"] / "metadata/c7_canonical_summary.json")
    return report


def main() -> int:
    """Small command-line surface used by clean-kernel notebook workers."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    contrast_parser = subparsers.add_parser("contrast")
    contrast_parser.add_argument("contrast_id")
    contrast_parser.add_argument("hierarchy", choices=["major", "fibroblast"])
    subparsers.add_parser("assemble")
    args = parser.parse_args()
    result = (
        run_contrast(args.contrast_id, args.hierarchy)
        if args.command == "contrast"
        else assemble_canonical_outputs()
    )
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
