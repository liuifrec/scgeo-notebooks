"""Sparse, descriptive-only Phase 2 helpers for GSE249479.

The functions in this module never densify a complete cells-by-genes matrix.
They operate on CSR counts, materialize only one-dimensional summaries or small
signature blocks, and enforce 20/24 GiB warning/stop RSS thresholds.
"""

from __future__ import annotations

import gc
import hashlib
import importlib.metadata as im
import json
import os
import platform
import re
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path("/home/liuyuchen/hsc_memory_nature_2026/results/zeng2026_xenograft_cd34_rna_raw.h5ad")
DEFAULT_COMPACT = Path("/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_qc_hvg.h5ad")
DEFAULT_OUTPUT = ROOT / "results/public_validation/gse249479_dataset_b"
WARNING_GIB = 20.0
HARD_STOP_GIB = 24.0
SOURCE_DOI = "10.1038/s41586-026-10522-7"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def rss_gib() -> float:
    return psutil.Process(os.getpid()).memory_info().rss / 1024**3


def close_and_collect() -> None:
    try:
        import matplotlib.pyplot as plt

        plt.close("all")
    except Exception:
        pass
    gc.collect()


class MemoryLog:
    def __init__(self, stage: str):
        self.stage = stage
        self.rows: list[dict[str, Any]] = []
        self.peak_rss_gib = rss_gib()

    def check(self, operation: str, point: str) -> float:
        value = rss_gib()
        self.peak_rss_gib = max(self.peak_rss_gib, value)
        status = "hard_stop" if value >= HARD_STOP_GIB else "warning" if value >= WARNING_GIB else "ok"
        self.rows.append(
            {
                "stage": self.stage,
                "operation": operation,
                "point": point,
                "rss_gib": value,
                "warning_threshold_gib": WARNING_GIB,
                "hard_stop_gib": HARD_STOP_GIB,
                "status": status,
                "timestamp_utc": utc_now(),
            }
        )
        if value >= HARD_STOP_GIB:
            raise MemoryError(f"RSS {value:.2f} GiB reached the {HARD_STOP_GIB:.0f} GiB hard stop")
        return value

    @contextmanager
    def operation(self, name: str):
        close_and_collect()
        self.check(name, "before")
        try:
            yield
        finally:
            close_and_collect()
            self.check(name, "after")

    def write(self, path: Path) -> None:
        frame = pd.DataFrame(self.rows)
        if path.exists():
            old = pd.read_csv(path)
            frame = pd.concat([old, frame], ignore_index=True)
        atomic_csv(frame, path)


def paths() -> dict[str, Path]:
    return {
        "input": Path(os.environ.get("SCGEO_GSE249479_H5AD", DEFAULT_INPUT)).resolve(),
        "compact": Path(os.environ.get("SCGEO_GSE249479_COMPACT_H5AD", DEFAULT_COMPACT)).resolve(),
        "output": Path(os.environ.get("SCGEO_GSE249479_OUTPUT_DIR", DEFAULT_OUTPUT)).resolve(),
        "source_repo": Path(os.environ.get("SCGEO_SOURCE_REPO", ROOT.parent / "scgeo")).resolve(),
    }


def ensure_output_tree(output: Path) -> None:
    for name in ["audit", "figures", "alt_text", "metadata", "version_records", "execution", "executed_notebooks", "signatures"]:
        (output / name).mkdir(parents=True, exist_ok=True)


def atomic_json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(frame: pd.DataFrame, path: Path, **kwargs: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".csv.gz" if path.name.endswith(".csv.gz") else path.suffix
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}{suffix}")
    frame.to_csv(tmp, index=False, **kwargs)
    os.replace(tmp, path)


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(path: Path) -> str | None:
    import subprocess

    try:
        return subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"], check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        ).stdout.strip()
    except Exception:
        return None


def vector_sum(matrix: sparse.spmatrix, axis: int) -> np.ndarray:
    """Return a dense one-dimensional summary, never a dense 2-D expression matrix."""
    return np.array(matrix.sum(axis=axis), dtype=np.float64).ravel()


def classify_identifiers(names: pd.Index) -> dict[str, Any]:
    values = names.astype(str)
    ensembl = np.fromiter((bool(re.fullmatch(r"ENS[A-Z]*G\d+(?:\.\d+)?", x)) for x in values), bool)
    symbol = np.fromiter((bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", x)) and not x.startswith("ENS") for x in values), bool)
    if ensembl.mean() >= 0.95:
        kind = "ensembl_ids"
    elif symbol.mean() >= 0.95:
        kind = "gene_symbols_with_gencode_locus_names"
    else:
        kind = "mixed_or_unclassified"
    return {
        "identifier_type": kind,
        "n_identifiers": len(values),
        "n_ensembl_like": int(ensembl.sum()),
        "n_symbol_like": int(symbol.sum()),
        "unique": bool(names.is_unique),
        "examples": values[:20].tolist(),
        "remapping_performed": False,
    }


def gene_masks(names: pd.Index) -> dict[str, np.ndarray]:
    upper = names.astype(str).str.upper()
    return {
        "mitochondrial": np.array(upper.str.startswith("MT-"), dtype=bool),
        "ribosomal": np.array(upper.str.match(r"^RP[SL]\d+[A-Z0-9-]*$"), dtype=bool),
        "haemoglobin": np.array(upper.isin(["HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBQ1"]), dtype=bool),
    }


def fraction_for_mask(X: sparse.csr_matrix, mask: np.ndarray, totals: np.ndarray) -> np.ndarray:
    numerator = vector_sum(X[:, mask], axis=1) if mask.any() else np.zeros(X.shape[0], dtype=np.float64)
    return np.divide(numerator, totals, out=np.zeros_like(numerator), where=totals > 0) * 100.0


def add_library_outlier_flags(obs: pd.DataFrame) -> None:
    obs["extreme_library_size_flag"] = False
    for _, idx in obs.groupby("condition", observed=False).groups.items():
        vals = np.log1p(obs.loc[idx, "total_counts"].to_numpy(dtype=float))
        med = np.median(vals)
        mad = np.median(np.abs(vals - med))
        if mad > 0:
            robust_z = 0.67448975 * (vals - med) / mad
            obs.loc[idx, "extreme_library_size_flag"] = np.abs(robust_z) > 5.0


def exclusion_reasons(obs: pd.DataFrame) -> pd.Series:
    out: list[str] = []
    for zero, genes, mito in zip(obs["total_counts"].eq(0), obs["detected_genes"], obs["pct_counts_mito"]):
        reasons = []
        if zero:
            reasons.append("zero_total_counts")
        if genes <= 1000:
            reasons.append("detected_genes_le_1000")
        if mito >= 18.0:
            reasons.append("mitochondrial_fraction_ge_18pct")
        out.append(";".join(reasons))
    return pd.Series(out, index=obs.index, dtype="object")


def distribution_table(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["total_counts", "detected_genes", "pct_counts_mito", "pct_counts_ribo", "pct_counts_haemoglobin"]
    for condition, frame in [("ALL", obs), *list(obs.groupby("condition", observed=False))]:
        for metric in metrics:
            values = frame[metric].to_numpy(dtype=float)
            rows.append({
                "condition": str(condition), "metric": metric, "n": len(values),
                "min": float(np.min(values)), "q01": float(np.quantile(values, 0.01)),
                "q05": float(np.quantile(values, 0.05)), "median": float(np.median(values)),
                "q95": float(np.quantile(values, 0.95)), "q99": float(np.quantile(values, 0.99)),
                "max": float(np.max(values)), "mean": float(np.mean(values)),
            })
    return pd.DataFrame(rows)


def attach_souporcell(obs: pd.DataFrame, mapping_path: Path) -> dict[str, Any]:
    for col in ["souporcell_assignment", "souporcell_clade", "souporcell_group"]:
        obs[col] = pd.Series(pd.NA, index=obs.index, dtype="string")
    if not mapping_path.exists():
        return {"mapping_available": False, "n_mapped": 0, "mapping_path": str(mapping_path)}
    mapping = pd.read_csv(mapping_path, dtype=str).rename(columns={
        "Original_barcode": "obs_name", "SoupAssignment": "souporcell_assignment",
        "SoupClade": "souporcell_clade", "Group": "souporcell_group",
    }).set_index("obs_name")
    overlap = obs.index.intersection(mapping.index)
    obs.loc[overlap, ["souporcell_assignment", "souporcell_clade", "souporcell_group"]] = mapping.loc[
        overlap, ["souporcell_assignment", "souporcell_clade", "souporcell_group"]
    ].astype("string")
    return {
        "mapping_available": True, "n_mapping_rows": int(len(mapping)), "n_mapped": int(len(overlap)),
        "mapping_fraction": float(len(overlap) / len(obs)), "mapping_path": str(mapping_path),
        "interpretation": "Genetic lineage descriptor only; not a biological replicate unit.",
    }


def diagnostic_qc_figure(obs: pd.DataFrame, output: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    specs = [
        ("total_counts", "log10 total counts"), ("detected_genes", "Detected genes"),
        ("pct_counts_mito", "Mitochondrial fraction (%)"), ("pct_counts_ribo", "Ribosomal fraction (%)"),
    ]
    colors = {"PBS": "#4c78a8", "TNF": "#e45756", "LPS": "#72b7b2"}
    for ax, (metric, label) in zip(axes.ravel(), specs):
        for condition, frame in obs.groupby("condition", observed=False):
            vals = frame[metric].to_numpy(dtype=float)
            if metric == "total_counts":
                vals = np.log10(np.maximum(vals, 1))
            ax.hist(vals, bins=60, histtype="step", density=True, linewidth=1.2, color=colors.get(str(condition)), label=str(condition))
        ax.set_xlabel(label)
        ax.set_ylabel("Density")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("GSE249479 sparse QC distributions (diagnostic; cells are not replicates)")
    fig.tight_layout()
    for ext in ["png", "svg"]:
        fig.savefig(output / "figures" / f"01_sparse_qc_distributions.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    (output / "alt_text" / "01_sparse_qc_distributions.txt").write_text(
        "Four overlaid density panels compare PBS, TNF and LPS cell-level distributions for log10 total counts, detected genes, mitochondrial percentage and ribosomal percentage. The plot is diagnostic only; cells and libraries are not biological replicates.\n",
        encoding="utf-8",
    )


def compact_storage_info(path: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        return {
            "shape": [int(v) for v in x.attrs["shape"]],
            "encoding_type": x.attrs.get("encoding-type", "").decode() if isinstance(x.attrs.get("encoding-type"), bytes) else x.attrs.get("encoding-type"),
            "dtype": str(x["data"].dtype), "nnz": int(x["data"].size),
            "raw_present": "raw" in handle, "layers": sorted(handle.get("layers", {}).keys()),
            "file_size_bytes": path.stat().st_size,
        }


def write_phase_metadata(stage: str, output: Path, p: dict[str, Path], memory: MemoryLog, extra: dict[str, Any]) -> None:
    packages = {}
    for name in ["anndata", "h5py", "numpy", "pandas", "scipy", "psutil", "scanpy", "matplotlib", "nbformat", "nbclient"]:
        try:
            packages[name] = im.version(name)
        except im.PackageNotFoundError:
            packages[name] = None
    payload = {
        "stage": stage, "timestamp_utc": utc_now(), "inference_status": "descriptive_only",
        "python": sys.version, "platform": platform.platform(), "packages": packages,
        "notebook_repository_commit": git_commit(ROOT), "frozen_scgeo_commit": git_commit(p["source_repo"]),
        "scgeo_modified": False, "frozen_thresholds_tuned": False,
        "peak_rss_gib": memory.peak_rss_gib, "warning_threshold_gib": WARNING_GIB, "hard_stop_gib": HARD_STOP_GIB,
        "forbidden_steps_run": [], "full_matrix_densified": False,
    }
    payload.update(extra)
    atomic_json(payload, output / "metadata" / f"{stage}_metadata.json")
    atomic_json({"packages": packages, "python": platform.python_version()}, output / "version_records" / f"{stage}_versions.json")


def run_sparse_qc_and_compact() -> dict[str, Any]:
    p = paths()
    output, source, compact_path = p["output"], p["input"], p["compact"]
    ensure_output_tree(output)
    memory = MemoryLog("01_sparse_qc_and_compact_object")
    mapping_path = output / "audit" / "souporcell_cell_mapping.csv.gz"
    source_sha_before = sha256(source)

    with memory.operation("load_sparse_source"):
        adata = ad.read_h5ad(source)
        if not sparse.issparse(adata.X):
            raise RuntimeError("Source X is not sparse; refusing Phase 2")
        X = adata.X.tocsr(copy=False)
        if X.dtype != np.float32:
            X = X.astype(np.float32)
        adata.X = X

    with memory.operation("cell_qc_vectors"):
        totals = vector_sum(X, axis=1)
        detected = np.diff(X.indptr).astype(np.int32, copy=False)
        masks = gene_masks(adata.var_names)
        obs = adata.obs.copy()
        obs["original_obs_name"] = adata.obs_names.astype(str)
        obs["total_counts"] = totals
        obs["detected_genes"] = detected
        obs["pct_counts_mito"] = fraction_for_mask(X, masks["mitochondrial"], totals)
        obs["pct_counts_ribo"] = fraction_for_mask(X, masks["ribosomal"], totals)
        obs["pct_counts_haemoglobin"] = fraction_for_mask(X, masks["haemoglobin"], totals)
        add_library_outlier_flags(obs)
        obs["low_complexity_flag"] = obs["detected_genes"] <= 1000
        obs["exclusion_reason"] = exclusion_reasons(obs)
        obs["retained_for_compact"] = obs["exclusion_reason"].eq("")
        obs["inference_status"] = "descriptive_only"
        mapping_summary = attach_souporcell(obs, mapping_path)

    with memory.operation("gene_qc_vectors"):
        gene_totals = vector_sum(X, axis=0)
        gene_detected = np.array(X.getnnz(axis=0), dtype=np.int64).ravel()
        identifier = classify_identifiers(adata.var_names)
        gene_table = pd.DataFrame(index=adata.var_names)
        gene_table["gene_identifier"] = adata.var_names.astype(str)
        gene_table["cells_detected"] = gene_detected
        gene_table["total_counts"] = gene_totals
        gene_table["mitochondrial"] = masks["mitochondrial"]
        gene_table["ribosomal"] = masks["ribosomal"]
        gene_table["haemoglobin"] = masks["haemoglobin"]
        gene_table["retained_feature"] = (gene_detected >= 10) & (gene_totals > 0)
        gene_table["feature_filter_reason"] = np.where(gene_table["retained_feature"], "", "detected_in_fewer_than_10_cells_or_zero_total")

    atomic_csv(distribution_table(obs), output / "audit" / "01_qc_distribution_summary.csv")
    cell_columns = ["original_obs_name", "condition", "sample", "condition_batch", "total_counts", "detected_genes", "pct_counts_mito", "pct_counts_ribo", "pct_counts_haemoglobin", "extreme_library_size_flag", "low_complexity_flag", "retained_for_compact", "exclusion_reason", "souporcell_assignment", "souporcell_clade", "souporcell_group"]
    atomic_csv(obs[cell_columns], output / "audit" / "01_retained_excluded_cells.csv.gz", compression="gzip")
    atomic_csv(obs[cell_columns], output / "audit" / "01_cell_qc.csv.gz", compression="gzip")
    atomic_csv(gene_table.reset_index(drop=True), output / "audit" / "01_gene_filter_table.csv.gz", compression="gzip")
    filter_decisions = pd.DataFrame([
        {"criterion": "total_counts > 0", "action": "exclude if failed", "basis": "empty droplets/cells cannot be analysed"},
        {"criterion": "detected_genes > 1000", "action": "exclude if failed", "basis": "primary-study CB xenograft scMultiome RNA threshold"},
        {"criterion": "mitochondrial fraction < 18%", "action": "exclude if failed", "basis": "primary-study CB xenograft scMultiome RNA threshold"},
        {"criterion": "within-condition |robust z(log1p counts)| > 5", "action": "flag only; do not exclude", "basis": "diagnostic extreme library-size flag"},
        {"criterion": "gene detected in >=10 retained-or-source cells", "action": "retain feature", "basis": "conservative sparse feature filter"},
    ])
    atomic_csv(filter_decisions, output / "audit" / "01_qc_filter_decisions.csv")
    atomic_json(identifier, output / "audit" / "01_identifier_audit.json")
    diagnostic_qc_figure(obs, output)

    retained_cells = obs["retained_for_compact"].to_numpy(dtype=bool)
    retained_genes = gene_table["retained_feature"].to_numpy(dtype=bool)
    with memory.operation("subset_compact_sparse"):
        compact = adata[retained_cells, retained_genes].copy()
        compact.obs = obs.loc[compact.obs_names].copy()
        compact.var = gene_table.loc[compact.var_names].copy()
        compact.X = compact.X.tocsr().astype(np.float32, copy=False)

    with memory.operation("sparse_raw_count_hvg"):
        os.environ.setdefault("NUMBA_CACHE_DIR", str(output / "_numba_cache"))
        import scanpy as sc

        n_top = min(3000, compact.n_vars)
        sc.pp.highly_variable_genes(compact, flavor="seurat_v3", n_top_genes=n_top, subset=False, check_values=True)
        compact.var["hvg_method"] = "seurat_v3_raw_counts_sparse"
        compact.uns["hvg"] = {"method": "seurat_v3_raw_counts_sparse", "n_top_genes": int(n_top)}

    compact.uns["dataset_b_inference_status"] = "descriptive_only"
    compact.uns["replication_note"] = "No valid biological replicate identifier. Cells, conditions, sample, condition_batch, libraries and SouporCell clades are not replicates."
    compact.uns["source_h5ad"] = str(source)
    compact.uns["source_sha256"] = source_sha_before
    compact.uns["source_doi"] = SOURCE_DOI
    compact.uns["X_representation"] = "raw_sparse_counts"
    compact.uns["full_matrix_densified"] = False
    compact.uns["souporcell_mapping"] = mapping_summary
    compact.uns["filter_policy"] = {"detected_genes": ">1000", "pct_counts_mito": "<18", "feature_min_cells": 10, "extreme_library_size": "flag_only"}

    with memory.operation("atomic_compact_write"):
        compact_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = compact_path.with_name(f".{compact_path.name}.tmp.{os.getpid()}.h5ad")
        compact.write_h5ad(tmp, compression="gzip")
        backed = ad.read_h5ad(tmp, backed="r")
        try:
            backed_sparse = sparse.issparse(backed.X) or getattr(backed.X, "format", None) in {"csr", "csc"}
            if backed.shape != compact.shape or not backed_sparse:
                raise RuntimeError("Backed validation failed for temporary compact H5AD")
        finally:
            backed.file.close()
        os.replace(tmp, compact_path)

    with memory.operation("compact_validation_and_checksum"):
        validation = compact_storage_info(compact_path)
        validation.update({
            "path": str(compact_path), "sha256": sha256(compact_path),
            "source_path": str(source), "source_sha256_before": source_sha_before,
            "source_sha256_after": sha256(source), "source_immutable": True,
            "inference_status": "descriptive_only", "validation_status": "passed",
            "n_hvg": int(compact.var["highly_variable"].sum()),
        })
        atomic_json(validation, output / "audit" / "compact_object_validation_report.json")

    summary = {
        "status": "completed", "input_shape": [int(adata.n_obs), int(adata.n_vars)],
        "retained_cells": int(retained_cells.sum()), "excluded_cells": int((~retained_cells).sum()),
        "retained_genes": int(retained_genes.sum()), "excluded_genes": int((~retained_genes).sum()),
        "identifier_type": identifier["identifier_type"], "compact_h5ad": str(compact_path),
        "compact_shape": validation["shape"], "compact_sha256": validation["sha256"],
        "mapping_summary": mapping_summary, "peak_rss_gib": memory.peak_rss_gib,
        "inference_status": "descriptive_only",
    }
    atomic_json(summary, output / "audit" / "01_sparse_qc_summary.json")
    memory.write(output / "audit" / "phase2_memory_log.csv")
    write_phase_metadata("01_sparse_qc_and_compact_object", output, p, memory, summary)
    del compact, adata, X
    close_and_collect()
    return summary


def signature_score(X: sparse.csr_matrix, totals: np.ndarray, indices: np.ndarray, chunk: int = 32) -> np.ndarray:
    if len(indices) == 0:
        return np.full(X.shape[0], np.nan, dtype=np.float32)
    accum = np.zeros(X.shape[0], dtype=np.float64)
    scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=np.float64), where=totals > 0)
    for start in range(0, len(indices), chunk):
        block = X[:, indices[start : start + chunk]].tocoo(copy=True)
        block.data = np.log1p(block.data.astype(np.float64) * scale[block.row])
        accum += vector_sum(block, axis=1)
    return (accum / len(indices)).astype(np.float32)


def run_annotation_and_signatures() -> dict[str, Any]:
    p = paths()
    output, compact_path = p["output"], p["compact"]
    ensure_output_tree(output)
    memory = MemoryLog("02_annotation_and_signatures")
    with memory.operation("load_compact_sparse"):
        adata = ad.read_h5ad(compact_path)
        if not sparse.issparse(adata.X):
            raise RuntimeError("Compact X is not sparse")
        X = adata.X.tocsr(copy=False)
        totals = vector_sum(X, axis=1)

    spec = json.loads((ROOT / "configs/gse249479_signatures_v1.json").read_text(encoding="utf-8"))
    name_to_idx = {str(g).upper(): i for i, g in enumerate(adata.var_names)}
    definitions = []
    with memory.operation("sparse_signature_scoring"):
        for name, item in spec["signatures"].items():
            requested = [str(g) for g in item["genes"]]
            available = [g for g in requested if g.upper() in name_to_idx]
            idx = np.array([name_to_idx[g.upper()] for g in available], dtype=np.int64)
            adata.obs[f"score_{name}"] = signature_score(X, totals, idx)
            definitions.append({
                "signature": name, "provenance_class": item["provenance_class"],
                "source_table": item.get("source_table", item.get("source_signature", "")),
                "n_requested": len(requested), "n_available": len(available),
                "available_genes": ";".join(available), "missing_genes": ";".join(g for g in requested if g not in available),
                "scoring_method": "mean log1p(counts_per_10000) over available genes; sparse signature blocks only",
            })
    definitions_frame = pd.DataFrame(definitions)
    atomic_csv(definitions_frame, output / "signatures" / "signature_definitions.csv")

    lineage = {
        "HSC_quiescent": "score_hsc_quiescent_markers", "activated_HSC": "score_activated_hsc_markers",
        "MPP_progenitor": "score_mpp_progenitor_markers", "myeloid": "score_myeloid_markers",
        "megakaryocyte_erythroid": "score_megakaryocyte_erythroid_markers", "lymphoid": "score_lymphoid_markers",
    }
    with memory.operation("conservative_marker_annotation"):
        score_matrix = np.column_stack([adata.obs[col].to_numpy(dtype=float) for col in lineage.values()])
        z = (score_matrix - np.nanmedian(score_matrix, axis=0)) / np.maximum(np.nanstd(score_matrix, axis=0), 1e-8)
        order = np.argsort(z, axis=1)
        top = order[:, -1]
        margin = z[np.arange(len(z)), top] - z[np.arange(len(z)), order[:, -2]]
        peak = z[np.arange(len(z)), top]
        labels = np.array(list(lineage), dtype=object)[top]
        confident = (peak >= 0.5) & (margin >= 0.25)
        adata.obs["marker_inferred_label"] = np.where(confident, labels, "ambiguous_HSPC")
        adata.obs["annotation_confidence"] = np.where(confident & (margin >= 0.75), "moderate", np.where(confident, "low", "ambiguous"))
        adata.obs["annotation_score_margin"] = margin.astype(np.float32)
        adata.obs["annotation_note"] = "Marker-based descriptive label; not a primary-study label and not independently validated."

    score_cols = [c for c in adata.obs if c.startswith("score_")]
    summary_rows = []
    for condition, frame in adata.obs.groupby("condition", observed=False):
        for col in score_cols:
            values = frame[col].to_numpy(dtype=float)
            summary_rows.append({"condition": str(condition), "signature": col.removeprefix("score_"), "n_cells": len(values), "median": float(np.nanmedian(values)), "q25": float(np.nanquantile(values, .25)), "q75": float(np.nanquantile(values, .75))})
    atomic_csv(pd.DataFrame(summary_rows), output / "audit" / "02_signature_condition_summary.csv")
    confidence = adata.obs.groupby(["marker_inferred_label", "annotation_confidence"], observed=False).size().reset_index(name="n_cells")
    atomic_csv(confidence, output / "audit" / "02_annotation_confidence_table.csv")

    cell_out = adata.obs[["condition", "marker_inferred_label", "annotation_confidence", "annotation_score_margin", "souporcell_clade", *score_cols]].copy()
    cell_out.insert(0, "original_obs_name", adata.obs_names.astype(str))
    atomic_csv(cell_out, output / "audit" / "02_cell_signature_scores.csv.gz", compression="gzip")

    representation_plan = pd.DataFrame([
        {"representation": "PCA20", "execute_now": False, "feasibility": "feasible_with_sparse_centering-free_or_incremental_method", "batch_covariate": "none", "resource_plan": "HVG-only float32 input; fresh process; checkpoint and monitor RSS"},
        {"representation": "PCA30", "execute_now": False, "feasibility": "feasible_with_sparse_centering-free_or_incremental_method", "batch_covariate": "none", "resource_plan": "HVG-only float32 input; fresh process; checkpoint and monitor RSS"},
        {"representation": "PCA50", "execute_now": False, "feasibility": "feasible_with_sparse_centering-free_or_incremental_method", "batch_covariate": "none", "resource_plan": "HVG-only float32 input; fresh process; checkpoint and monitor RSS"},
        {"representation": "diffusion_map", "execute_now": False, "feasibility": "defer_until_sparse_neighbour_graph_memory_estimate", "batch_covariate": "none", "resource_plan": "Build from a selected PCA checkpoint only; sparse kNN graph; fresh process"},
        {"representation": "scVI_optional", "execute_now": False, "feasibility": "possible_separate_monitored_pass", "batch_covariate": "none", "resource_plan": "Raw counts; no condition covariate; small minibatches and workers; hard stop at 24 GiB"},
        {"representation": "UMAP_display_only", "execute_now": False, "feasibility": "feasible_after_representation_selection", "batch_covariate": "none", "resource_plan": "Display only; never infer replication from visual separation"},
        {"representation": "Harmony", "execute_now": False, "feasibility": "not_planned_no_genuine_noncondition_batch", "batch_covariate": "prohibited", "resource_plan": "Do not use condition, sample or condition_batch"},
        {"representation": "Scanorama", "execute_now": False, "feasibility": "not_planned_no_genuine_noncondition_batch", "batch_covariate": "prohibited", "resource_plan": "Do not use condition, sample or condition_batch"},
    ])
    atomic_csv(representation_plan, output / "audit" / "02_representation_plan.csv")

    with memory.operation("write_annotation_metadata_only"):
        tmp = compact_path.with_name(f".{compact_path.name}.annotation.tmp.{os.getpid()}.h5ad")
        adata.uns["signature_scoring"] = {"status": "descriptive_only", "config": "configs/gse249479_signatures_v1.json", "no_inferential_p_values": True}
        adata.uns["representation_plan_status"] = "prepared_not_executed"
        adata.write_h5ad(tmp, compression="gzip")
        check = ad.read_h5ad(tmp, backed="r")
        try:
            backed_sparse = sparse.issparse(check.X) or getattr(check.X, "format", None) in {"csr", "csc"}
            if check.shape != adata.shape or not backed_sparse:
                raise RuntimeError("Annotated compact-object validation failed")
        finally:
            check.file.close()
        os.replace(tmp, compact_path)

    validation = compact_storage_info(compact_path)
    validation.update({"path": str(compact_path), "sha256": sha256(compact_path), "validation_status": "passed_after_annotation", "inference_status": "descriptive_only", "n_hvg": int(adata.var["highly_variable"].sum())})
    atomic_json(validation, output / "audit" / "compact_object_validation_report.json")
    summary = {
        "status": "completed", "shape": [int(adata.n_obs), int(adata.n_vars)],
        "label_counts": adata.obs["marker_inferred_label"].value_counts().to_dict(),
        "confidence_counts": adata.obs["annotation_confidence"].value_counts().to_dict(),
        "signatures_scored": list(spec["signatures"]), "inferential_p_values": False,
        "representations_executed": [], "peak_rss_gib": memory.peak_rss_gib,
        "compact_sha256": validation["sha256"], "inference_status": "descriptive_only",
    }
    atomic_json(summary, output / "audit" / "02_annotation_summary.json")
    memory.write(output / "audit" / "phase2_memory_log.csv")
    write_phase_metadata("02_annotation_and_signatures", output, p, memory, summary)
    del adata, X
    close_and_collect()
    return summary
