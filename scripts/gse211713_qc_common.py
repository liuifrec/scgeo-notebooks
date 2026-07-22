#!/usr/bin/env python3
"""Sparse QC, per-mouse doublet assessment, and compact-object construction."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "configs/gse211713_dataset_c_v1.json").read_text())
DATA_DIR = Path(os.environ.get("SCGEO_GSE211713_DATA_DIR", "/home/liuyuchen/data/gse211713")).resolve()
RESULTS = ROOT / "results/public_validation/gse211713_dataset_c"
RAW = DATA_DIR / "gse211713_all20_raw_sparse.h5ad"
COMPACT = DATA_DIR / "gse211713_revision_qc_hvg.h5ad"
WARNING_BYTES = 20 * 1024**3
STOP_BYTES = 24 * 1024**3
SEED = 1729
QC_RULES = {
    "minimum_total_counts": 500,
    "minimum_detected_genes": 250,
    "maximum_mito_fraction": 0.25,
    "minimum_gene_count_ratio": 0.05,
    "extreme_mad_multiplier": 5.0,
    "haemoglobin_flag_fraction": 0.20,
    "ribosomal_flag_fraction": 0.60,
    "scrublet_expected_doublet_rate": 0.06,
    "scrublet_score_threshold": 0.25,
    "scrublet_very_high_score": 0.50,
    "hvg_count": 3500,
    "minimum_cells_per_retained_gene": 3,
    "minimum_cells_per_hvg_candidate": 20,
}


def rss_gib() -> float:
    return psutil.Process().memory_info().rss / 1024**3


def record(records: list[dict[str, Any]], stage: str, started: float) -> None:
    rss = rss_gib()
    records.append({"stage": stage, "rss_gib": rss, "elapsed_seconds": time.perf_counter() - started, "timestamp_utc": datetime.now(timezone.utc).isoformat()})
    if rss * 1024**3 >= STOP_BYTES:
        raise MemoryError(f"QC hard stop at {rss:.2f} GiB during {stage}")
    if rss * 1024**3 >= WARNING_BYTES:
        print(f"WARNING: QC RSS {rss:.2f} GiB during {stage}", file=sys.stderr)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    frame.to_csv(temporary, index=False, compression=compression)
    os.replace(temporary, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def mad_bounds(values: np.ndarray, multiplier: float) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    median = float(np.median(finite)); mad = float(np.median(np.abs(finite - median)))
    scaled = 1.4826 * mad
    return median - multiplier * scaled, median + multiplier * scaled


def compute_metrics() -> dict[str, Any]:
    started = time.perf_counter(); memory: list[dict[str, Any]] = []
    adata = ad.read_h5ad(RAW)
    if not sparse.isspmatrix_csr(adata.X):
        raise RuntimeError("Raw X must be CSR sparse")
    record(memory, "raw_loaded", started)
    symbols = adata.var["gene_symbol"].astype(str)
    mt = symbols.str.startswith("mt-").to_numpy()
    ribo = symbols.str.match(r"^Rp[sl]").to_numpy()
    hb = symbols.str.match(r"^Hb[ab](?:-|[0-9])").to_numpy()
    total = np.asarray(adata.X.sum(axis=1)).ravel().astype(np.float64)
    genes = adata.X.getnnz(axis=1).astype(np.int32)
    def masked_fraction(mask: np.ndarray) -> np.ndarray:
        numerator = np.asarray(adata.X[:, mask].sum(axis=1)).ravel() if mask.any() else np.zeros(adata.n_obs)
        return np.divide(numerator, total, out=np.zeros_like(total), where=total > 0)
    obs = adata.obs.copy()
    obs["total_counts"] = total
    obs["detected_genes"] = genes
    obs["mitochondrial_fraction"] = masked_fraction(mt)
    obs["ribosomal_fraction"] = masked_fraction(ribo)
    obs["haemoglobin_fraction"] = masked_fraction(hb)
    obs["gene_count_ratio"] = np.divide(genes, total, out=np.zeros_like(total), where=total > 0)
    record(memory, "sparse_qc_metrics", started)
    obs["extreme_high_library"] = False
    obs["extreme_low_library"] = False
    for gsm, index in obs.groupby("gsm", observed=False).groups.items():
        idx = obs.index.get_indexer(index)
        low_c, high_c = mad_bounds(np.log1p(total[idx]), QC_RULES["extreme_mad_multiplier"])
        low_g, high_g = mad_bounds(np.log1p(genes[idx]), QC_RULES["extreme_mad_multiplier"])
        obs.loc[index, "extreme_high_library"] = (np.log1p(total[idx]) > high_c) | (np.log1p(genes[idx]) > high_g)
        obs.loc[index, "extreme_low_library"] = (np.log1p(total[idx]) < low_c) | (np.log1p(genes[idx]) < low_g)
    obs["low_counts"] = obs["total_counts"] < QC_RULES["minimum_total_counts"]
    obs["low_genes"] = obs["detected_genes"] < QC_RULES["minimum_detected_genes"]
    obs["high_mito"] = obs["mitochondrial_fraction"] > QC_RULES["maximum_mito_fraction"]
    obs["low_complexity"] = obs["gene_count_ratio"] < QC_RULES["minimum_gene_count_ratio"]
    obs["high_haemoglobin"] = obs["haemoglobin_fraction"] > QC_RULES["haemoglobin_flag_fraction"]
    obs["high_ribosomal"] = obs["ribosomal_fraction"] > QC_RULES["ribosomal_flag_fraction"]
    obs["pre_doublet_retain"] = ~(obs["low_counts"] | obs["low_genes"] | obs["high_mito"] | obs["low_complexity"])
    path = RESULTS / "qc/per_cell_qc_pre_doublet.csv.gz"
    atomic_csv(obs.reset_index(), path, compression="gzip")
    summary = obs.groupby(["gsm", "irradiation_group"], observed=False).agg(
        cells=("gsm", "size"), median_counts=("total_counts", "median"), median_genes=("detected_genes", "median"),
        median_mito=("mitochondrial_fraction", "median"), pre_doublet_retained=("pre_doublet_retain", "sum"),
        low_counts=("low_counts", "sum"), low_genes=("low_genes", "sum"), high_mito=("high_mito", "sum"),
        low_complexity=("low_complexity", "sum"), extreme_high_library=("extreme_high_library", "sum"),
    ).reset_index()
    atomic_csv(summary, RESULTS / "qc/qc_summary_by_mouse_pre_doublet.csv")
    atomic_csv(pd.DataFrame([{"rule": key, "value": value, "source": "prespecified Dataset C revision rule; publication numeric threshold unavailable"} for key, value in QC_RULES.items()]), RESULTS / "qc/qc_filter_rules.csv")
    atomic_csv(pd.DataFrame(memory), RESULTS / "qc/memory_metrics_stage.csv")
    report = {"status": "passed", "cells": adata.n_obs, "genes": adata.n_vars, "peak_rss_gib": max(x["rss_gib"] for x in memory), "runtime_seconds": time.perf_counter() - started}
    atomic_json(report, RESULTS / "qc/metrics_report.json")
    print(json.dumps(report))
    return report


def doublet_sample(gsm: str) -> dict[str, Any]:
    local_packages = DATA_DIR / "python_packages"
    if local_packages.is_dir():
        sys.path.insert(0, str(local_packages))
    import importlib.metadata
    import scrublet as scr
    started = time.perf_counter(); memory: list[dict[str, Any]] = []
    path = DATA_DIR / "sample_h5ad" / f"{gsm}.h5ad"
    adata = ad.read_h5ad(path)
    record(memory, "sample_loaded", started)
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=QC_RULES["scrublet_expected_doublet_rate"], random_state=SEED)
    scores, _ = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30, verbose=False)
    threshold = QC_RULES["scrublet_score_threshold"]
    frame = pd.DataFrame({"cell_id": adata.obs_names, "gsm": gsm, "scrublet_score": scores, "scrublet_fixed_threshold_flag": scores >= threshold})
    output = RESULTS / "qc/doublet_scores" / f"{gsm}.csv.gz"
    atomic_csv(frame, output, compression="gzip")
    record(memory, "scrublet_complete", started)
    report = {
        "status": "passed", "gsm": gsm, "cells": adata.n_obs, "fixed_threshold": threshold,
        "automatic_threshold_reported_not_used": float(scrub.threshold_),
        "fixed_threshold_flags": int((scores >= threshold).sum()),
        "score_median": float(np.median(scores)), "score_max": float(np.max(scores)),
        "scrublet_version": importlib.metadata.version("scrublet"),
        "random_seed": SEED, "peak_rss_gib": max(x["rss_gib"] for x in memory),
        "runtime_seconds": time.perf_counter() - started,
    }
    atomic_json(report, RESULTS / "qc/doublet_scores" / f"{gsm}.json")
    print(json.dumps(report))
    return report


def run_doublets() -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(DATA_DIR / "python_packages") + os.pathsep + environment.get("PYTHONPATH", "")
    for values in CONFIG["samples"]:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--doublet-gsm", values[0]], check=True, env=environment)


def dispersion_hvgs(matrix: sparse.csr_matrix, var: pd.DataFrame, n_hvg: int) -> tuple[pd.DataFrame, np.ndarray]:
    n = matrix.shape[0]
    detected = matrix.getnnz(axis=0).astype(np.int32)
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    sumsq = np.asarray(matrix.power(2).sum(axis=0)).ravel().astype(np.float64)
    mean = totals / max(n, 1)
    variance = np.maximum(sumsq / max(n, 1) - mean**2, 0)
    dispersion = np.divide(variance, mean, out=np.zeros_like(variance), where=mean > 0)
    log_mean = np.log1p(mean); log_dispersion = np.log1p(dispersion)
    bins = pd.qcut(pd.Series(log_mean).rank(method="first"), q=20, labels=False, duplicates="drop").to_numpy()
    score = np.full(matrix.shape[1], -np.inf, dtype=np.float64)
    for value in np.unique(bins):
        idx = np.flatnonzero(bins == value); values = log_dispersion[idx]
        sd = values.std(ddof=1)
        score[idx] = (values - values.mean()) / sd if sd > 0 else 0
    symbols = var["gene_symbol"].astype(str)
    excluded = symbols.str.startswith("mt-") | symbols.str.match(r"^Rp[sl]") | symbols.str.match(r"^Hb[ab](?:-|[0-9])")
    candidates = (detected >= QC_RULES["minimum_cells_per_hvg_candidate"]) & ~excluded.to_numpy() & np.isfinite(score)
    eligible = np.flatnonzero(candidates)
    selected = eligible[np.argsort(score[eligible], kind="mergesort")[-n_hvg:]]
    hvg = np.zeros(matrix.shape[1], dtype=bool); hvg[selected] = True
    table = var.reset_index().copy()
    table["cells_detected"] = detected; table["total_counts"] = totals
    table["mean_counts"] = mean; table["variance"] = variance
    table["dispersion"] = dispersion; table["hvg_score"] = score
    table["hvg_candidate"] = candidates; table["highly_variable"] = hvg
    return table, hvg


def save_figures(obs: pd.DataFrame, summary: pd.DataFrame) -> None:
    figures = RESULTS / "qc/figures"; figures.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    groups = ["control", "10Gy", "17Gy"]; colors = ["#777777", "#E69F00", "#3B75AF"]
    for group, color in zip(groups, colors):
        part = obs[obs["irradiation_group"].eq(group)]
        axes[0,0].hist(np.log10(part.total_counts.clip(lower=1)), bins=80, histtype="step", density=True, color=color, label=group)
        axes[0,1].hist(part.detected_genes, bins=80, histtype="step", density=True, color=color)
        axes[1,0].hist(part.mitochondrial_fraction, bins=80, histtype="step", density=True, color=color)
    axes[0,0].set(xlabel="log10 total counts", ylabel="density", title="Library size"); axes[0,0].legend(frameon=False)
    axes[0,1].set(xlabel="detected genes", ylabel="density", title="Complexity")
    axes[1,0].set(xlabel="mitochondrial fraction", ylabel="density", title="Mitochondrial fraction")
    axes[1,1].barh(summary.gsm, summary.retained_fraction, color=summary.irradiation_group.map(dict(zip(groups, colors))))
    axes[1,1].set(xlabel="retained fraction", title="Retention by mouse", xlim=(0,1))
    fig.suptitle("GSE211713 sparse QC diagnostics; thresholds were prespecified, not fitted to the publication cell count")
    for extension in ("png", "svg"):
        fig.savefig(figures / f"gse211713_qc_diagnostics.{extension}", dpi=250, bbox_inches="tight", metadata={"Creator":"GSE211713 sparse QC","Date":None})
    plt.close(fig)
    alt = "Four-panel sparse-QC figure showing total-count, detected-gene and mitochondrial-fraction distributions for control, 10 Gy and 17 Gy cells, plus retained-cell fractions for each of 20 mouse libraries. Thresholds were prespecified and were not adjusted to reproduce the publication-retained count."
    (RESULTS / "qc/gse211713_qc_diagnostics_alt.txt").write_text(alt + "\n")


def build_compact() -> dict[str, Any]:
    started = time.perf_counter(); memory: list[dict[str, Any]] = []
    adata = ad.read_h5ad(RAW); record(memory, "raw_loaded", started)
    qc = pd.read_csv(RESULTS / "qc/per_cell_qc_pre_doublet.csv.gz").set_index("cell_id")
    doublets = pd.concat([pd.read_csv(RESULTS / "qc/doublet_scores" / f"{values[0]}.csv.gz") for values in CONFIG["samples"]], ignore_index=True).set_index("cell_id")
    qc = qc.join(doublets[["scrublet_score", "scrublet_fixed_threshold_flag"]], how="left")
    qc["doublet_high_confidence"] = (qc["scrublet_score"] >= QC_RULES["scrublet_very_high_score"]) | (qc["scrublet_fixed_threshold_flag"] & qc["extreme_high_library"])
    qc["retain"] = qc["pre_doublet_retain"] & ~qc["doublet_high_confidence"]
    reasons = []
    for row in qc.itertuples():
        item = []
        for column in ["low_counts", "low_genes", "high_mito", "low_complexity", "doublet_high_confidence"]:
            if getattr(row, column): item.append(column)
        reasons.append("retained" if not item else ";".join(item))
    qc["filter_reason"] = reasons
    retain = qc.reindex(adata.obs_names)["retain"].fillna(False).to_numpy(bool)
    retained = adata[retain].copy(); record(memory, "cells_filtered", started)
    detected = retained.X.getnnz(axis=0)
    gene_keep = detected >= QC_RULES["minimum_cells_per_retained_gene"]
    retained = retained[:, gene_keep].copy(); retained.X = retained.X.tocsr().astype(np.float32, copy=False)
    table, hvg = dispersion_hvgs(retained.X, retained.var, QC_RULES["hvg_count"])
    retained.var["cells_detected"] = table["cells_detected"].to_numpy()
    retained.var["total_counts"] = table["total_counts"].to_numpy()
    retained.var["hvg_score"] = table["hvg_score"].to_numpy()
    retained.var["highly_variable"] = hvg
    retained.obs = retained.obs.join(qc[["total_counts", "detected_genes", "mitochondrial_fraction", "ribosomal_fraction", "haemoglobin_fraction", "gene_count_ratio", "extreme_high_library", "extreme_low_library", "scrublet_score", "doublet_high_confidence", "retain", "filter_reason"]])
    retained.uns["qc"] = {
        "rules": QC_RULES, "publication_numeric_qc_thresholds_available": False,
        "publication_soupx_reproducible_from_deposit": False,
        "reason_soupx_not_reproduced": "empty-droplet/raw unfiltered matrices were not deposited among the processed MEX files",
        "hvg_method": "sparse raw-count mean-binned standardized log dispersion; condition-blind; mitochondrial/ribosomal/haemoglobin genes excluded from HVG candidacy",
        "raw_counts_location": "X", "normalized_matrix_stored": False,
        "biological_replicate_unit": "mouse_id/GSM",
        "source_raw_sha256": sha256_file(RAW),
    }
    record(memory, "genes_filtered_hvgs_selected", started)
    temporary = COMPACT.with_name(f".{COMPACT.name}.tmp.{os.getpid()}.h5ad")
    retained.write_h5ad(temporary, compression="gzip", compression_opts=4)
    backed = ad.read_h5ad(temporary, backed="r")
    try:
        if backed.shape != retained.shape or "CSRDataset" not in type(backed.X).__name__:
            raise RuntimeError("Compact backed validation failed")
    finally:
        backed.file.close()
    os.replace(temporary, COMPACT); record(memory, "compact_written_validated", started)
    atomic_csv(qc.reset_index(), RESULTS / "qc/retained_excluded_cells.csv.gz", compression="gzip")
    atomic_csv(table, RESULTS / "qc/gene_qc_hvg_table.csv")
    summary = qc.groupby(["gsm", "irradiation_group"], observed=False).agg(
        input_cells=("gsm", "size"), retained_cells=("retain", "sum"),
        doublet_high_confidence=("doublet_high_confidence", "sum"),
        median_counts=("total_counts", "median"), median_genes=("detected_genes", "median"),
        median_mito=("mitochondrial_fraction", "median"),
    ).reset_index()
    summary["retained_fraction"] = summary.retained_cells / summary.input_cells
    atomic_csv(summary, RESULTS / "qc/qc_summary_by_mouse.csv")
    doublet_summary = pd.DataFrame([json.loads((RESULTS / "qc/doublet_scores" / f"{values[0]}.json").read_text()) for values in CONFIG["samples"]])
    atomic_csv(doublet_summary, RESULTS / "qc/doublet_summary.csv")
    publication = int(CONFIG["dataset"]["publication_retained_cell_count"])
    discrepancy = {
        "geo_filtered_input_cells": adata.n_obs, "revision_retained_cells": retained.n_obs,
        "publication_retained_cells": publication, "revision_minus_publication": retained.n_obs - publication,
        "explanation": [
            "The publication gives no numeric cell-QC or doublet thresholds.",
            "The publication applied SoupX, but deposited processed files lack the empty-droplet/raw matrices required to reproduce it.",
            "Revision thresholds were prespecified independently and were not tuned to reproduce 102,869 cells.",
        ],
    }
    atomic_json(discrepancy, RESULTS / "qc/publication_count_discrepancy.json")
    atomic_csv(pd.DataFrame(memory), RESULTS / "qc/memory_compact_stage.csv")
    save_figures(qc, summary)
    report = {
        "status": "passed", "input_shape": list(adata.shape), "compact_shape": list(retained.shape),
        "retained_cells": retained.n_obs, "retained_genes": retained.n_vars,
        "hvg_count": int(hvg.sum()), "nnz": int(retained.X.nnz), "dtype": str(retained.X.dtype),
        "sparse_format": "csr", "sha256": sha256_file(COMPACT), "path": str(COMPACT),
        "peak_rss_gib": max(x["rss_gib"] for x in memory), "runtime_seconds": time.perf_counter() - started,
    }
    atomic_json(report, RESULTS / "qc/compact_h5ad_validation.json")
    print(json.dumps(report))
    return report


def run_all() -> None:
    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--compute-metrics"], check=True)
    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--run-doublets"], check=True)
    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--build-compact"], check=True)


def main() -> int:
    parser = argparse.ArgumentParser(); group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--compute-metrics", action="store_true"); group.add_argument("--doublet-gsm")
    group.add_argument("--run-doublets", action="store_true"); group.add_argument("--build-compact", action="store_true")
    group.add_argument("--all", action="store_true"); args = parser.parse_args()
    if args.compute_metrics: compute_metrics()
    elif args.doublet_gsm: doublet_sample(args.doublet_gsm)
    elif args.run_doublets: run_doublets()
    elif args.build_compact: build_compact()
    else: run_all()
    gc.collect(); return 0


if __name__ == "__main__":
    raise SystemExit(main())
