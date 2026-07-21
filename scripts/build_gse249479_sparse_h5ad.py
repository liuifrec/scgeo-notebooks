#!/usr/bin/env python3
"""Build a sparse raw-count AnnData object from public GSE249479 RNA MEX files."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from gse249479_memory_safe import (
    MemoryAudit,
    bytes_to_gb,
    ensure_output_tree,
    inspect_h5ad_storage,
    load_config,
    memory_threshold_bytes,
    configure_temp_environment,
    relative_or_absolute,
    require_repo_local_dataset_paths,
    repo_root,
    require_active_branch,
    resolve_path,
    sha256_file,
    write_json,
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_barcodes(path: Path) -> pd.Index:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if frame.shape[1] < 1:
        raise RuntimeError(f"Barcode file has no columns: {path}")
    return pd.Index(frame.iloc[:, 0].astype(str), name="barcode_original")


def read_features(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if frame.shape[1] < 2:
        raise RuntimeError(f"Feature file must contain at least gene ID and symbol columns: {path}")
    columns = ["gene_id", "gene_symbol", "feature_type", "genome"]
    frame = frame.iloc[:, : min(frame.shape[1], len(columns))].copy()
    frame.columns = columns[: frame.shape[1]]
    frame["gene_id"] = frame["gene_id"].astype(str)
    frame["gene_symbol"] = frame["gene_symbol"].astype(str)
    return frame


def make_unique(values: pd.Series, *, join: str = "-") -> tuple[pd.Index, dict[str, Any]]:
    seen: dict[str, int] = {}
    output: list[str] = []
    duplicate_count = 0
    empty_count = 0
    for raw in values.astype(str).tolist():
        base = raw.strip()
        if not base or base.lower() == "nan":
            empty_count += 1
            base = "missing_feature_name"
        if base not in seen:
            seen[base] = 0
            output.append(base)
            continue
        duplicate_count += 1
        seen[base] += 1
        candidate = f"{base}{join}{seen[base]}"
        while candidate in seen:
            seen[base] += 1
            candidate = f"{base}{join}{seen[base]}"
        seen[candidate] = 0
        output.append(candidate)
    return pd.Index(output), {
        "input_count": int(len(output)),
        "duplicate_count": int(duplicate_count),
        "empty_count": int(empty_count),
        "join": join,
        "method": "stable_suffix_only_for_duplicates",
        "var_names_modified": bool(duplicate_count or empty_count),
    }


def feature_signature(features: pd.DataFrame) -> pd.DataFrame:
    cols = [col for col in ["gene_id", "gene_symbol", "feature_type"] if col in features.columns]
    return features[cols].reset_index(drop=True)


def read_mex_matrix(matrix_path: Path, n_features: int, n_barcodes: int) -> sp.csr_matrix:
    with gzip.open(matrix_path, "rb") as handle:
        matrix = scipy.io.mmread(handle)
    if not sp.issparse(matrix):
        raise RuntimeError(f"Matrix Market reader did not return a sparse matrix for {matrix_path}")
    matrix = matrix.tocsr()
    if matrix.shape != (n_features, n_barcodes):
        raise RuntimeError(
            f"Matrix shape mismatch for {matrix_path.name}: observed {matrix.shape}, "
            f"expected {(n_features, n_barcodes)} from features/barcodes"
        )
    cell_by_gene = matrix.T.tocsr()
    cell_by_gene.sort_indices()
    return cell_by_gene


def source_path_from_manifest(row: pd.Series) -> Path:
    return Path(str(row["local_path"])).expanduser().resolve()


def load_manifest(output_dir: Path) -> pd.DataFrame:
    path = output_dir / "audit" / "download_manifest.csv"
    if not path.exists():
        raise FileNotFoundError(f"Download manifest is required before reconstruction: {path}")
    return pd.read_csv(path)


def sample_file_map(manifest: pd.DataFrame) -> dict[str, dict[str, Any]]:
    samples: dict[str, dict[str, Any]] = {}
    for row in manifest.itertuples(index=False):
        sample = str(row.sample_accession)
        samples.setdefault(sample, {"condition": str(row.condition), "files": {}, "checksums": {}})
        samples[sample]["files"][str(row.file_role)] = Path(str(row.local_path)).expanduser().resolve()
        samples[sample]["checksums"][str(row.file_role)] = {
            "filename": str(row.filename),
            "sha256": str(row.sha256),
            "size_bytes": int(row.observed_size_bytes),
            "official_url": str(row.official_url),
        }
    return samples


def serializable_source_files(samples: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    payload: dict[str, dict[str, Any]] = {}
    for accession, sample in samples.items():
        payload[accession] = {
            "condition": sample["condition"],
            "files": {
                role: {
                    **sample["checksums"][role],
                    "local_path": str(sample["files"][role]),
                }
                for role in sorted(sample["files"])
            },
        }
    return payload


def build_var(features: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    var = features.copy()
    preferred = "gene_symbol" if "gene_symbol" in var.columns else "gene_id"
    var_names, handling = make_unique(var[preferred])
    handling["preferred_var_name_source"] = preferred
    handling["gene_ids_preserved"] = "gene_id" in var.columns
    handling["gene_symbols_preserved"] = "gene_symbol" in var.columns
    var.index = var_names
    return var, handling


def append_memory_records(output_dir: Path, audit: MemoryAudit) -> Path:
    path = output_dir / "audit" / "reconstruction_memory.csv"
    new_records = audit.dataframe()
    new_records.insert(0, "script", Path(__file__).name)
    if path.exists():
        existing = pd.read_csv(path)
        combined = pd.concat([existing, new_records], ignore_index=True)
    else:
        combined = new_records
    combined.to_csv(path, index=False)
    return path


def run(
    *,
    root: Path,
    data_dir: Path,
    output_dir: Path,
    output_h5ad: Path,
    memory_threshold_gb: float | None = None,
) -> dict[str, Any]:
    config = load_config(root)
    configure_temp_environment(root, config, create=True)
    require_repo_local_dataset_paths(root, data_dir, output_h5ad)
    require_active_branch(root, config["required_git_branch"])
    ensure_output_tree(output_dir)
    threshold_bytes = int((memory_threshold_gb or float(config["default_memory_threshold_gb"])) * 1024**3)
    audit = MemoryAudit(threshold_bytes)
    manifest = load_manifest(output_dir)
    samples = sample_file_map(manifest)

    sample_matrices: list[sp.csr_matrix] = []
    obs_frames: list[pd.DataFrame] = []
    first_features: pd.DataFrame | None = None
    sample_summaries: list[dict[str, Any]] = []

    with audit.section("read_sparse_mex_by_sample"):
        for sample in config["geo"]["samples"]:
            accession = sample["accession"]
            if accession not in samples:
                raise RuntimeError(f"Downloaded files are missing from manifest for {accession}")
            files = samples[accession]["files"]
            for role in ("barcodes", "features", "matrix"):
                if role not in files:
                    raise RuntimeError(f"Downloaded {role} file is missing for {accession}")
            barcodes = read_barcodes(files["barcodes"])
            features = read_features(files["features"])
            if first_features is None:
                first_features = features
            elif not feature_signature(features).equals(feature_signature(first_features)):
                raise RuntimeError(
                    f"Feature axis differs for {accession}; Phase 1A refuses to realign genes without review."
                )
            matrix = read_mex_matrix(files["matrix"], len(features), len(barcodes))
            cell_ids = pd.Index([f"{accession}_{barcode}" for barcode in barcodes], name="cell_id")
            obs = pd.DataFrame(
                {
                    "condition": sample["condition"],
                    "library_id": accession,
                    "source_dataset": config["dataset"]["accession"],
                    "barcode_original": barcodes.astype(str).to_numpy(copy=False),
                    "cell_id": cell_ids.astype(str).to_numpy(copy=False),
                },
                index=cell_ids,
            )
            sample_matrices.append(matrix)
            obs_frames.append(obs)
            sample_summaries.append(
                {
                    "sample_accession": accession,
                    "condition": sample["condition"],
                    "n_cells": int(matrix.shape[0]),
                    "n_genes": int(matrix.shape[1]),
                    "nnz": int(matrix.nnz),
                    "matrix_dtype": str(matrix.dtype),
                    "sparse_format": matrix.getformat(),
                }
            )

    if first_features is None:
        raise RuntimeError("No features were read from the public MEX files.")

    with audit.section("assemble_sparse_anndata"):
        X = sp.vstack(sample_matrices, format="csr")
        X.sort_indices()
        obs = pd.concat(obs_frames, axis=0)
        if not obs.index.is_unique:
            raise RuntimeError("Globally unique cell identifiers are not unique after sample concatenation.")
        var, var_handling = build_var(first_features)
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["gse249479_public_reconstruction"] = {
            "timestamp_utc": now_utc(),
            "source_dataset": config["dataset"]["accession"],
            "source": "official GEO supplementary processed RNA MEX files",
            "no_normalization_performed": True,
            "no_downstream_analysis_performed": True,
            "biological_replicate_assumption": "not_made",
            "condition_and_library_id_are_not_biological_replicates": True,
            "download_manifest": str(output_dir / "audit" / "download_manifest.csv"),
            "source_files": serializable_source_files(samples),
            "var_name_handling": var_handling,
            "sample_summaries": sample_summaries,
        }

    with audit.section("write_sparse_h5ad"):
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        if output_h5ad.exists():
            raise RuntimeError(f"Refusing to overwrite existing reconstructed H5AD: {output_h5ad}")
        adata.write_h5ad(output_h5ad)

    with audit.section("verify_written_h5ad_without_materializing_X"):
        storage = inspect_h5ad_storage(output_h5ad)

    output_sha256 = sha256_file(output_h5ad)
    condition_counts = adata.obs["condition"].value_counts().rename_axis("condition").reset_index(name="n_cells")
    condition_counts.to_csv(output_dir / "audit" / "reconstruction_condition_counts.csv", index=False)
    memory_path = append_memory_records(output_dir, audit)

    summary = {
        "timestamp_utc": now_utc(),
        "output_h5ad": str(output_h5ad),
        "output_h5ad_display": relative_or_absolute(output_h5ad, root),
        "output_h5ad_size_bytes": int(output_h5ad.stat().st_size),
        "output_h5ad_sha256": output_sha256,
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "nnz": int(adata.X.nnz),
        "sparse_format": adata.X.getformat(),
        "matrix_dtype": str(adata.X.dtype),
        "raw_integer_like_counts_retained": bool(np.issubdtype(adata.X.dtype, np.integer)),
        "raw_present": False,
        "layers": [],
        "obs_fields": list(adata.obs.columns),
        "var_fields": list(adata.var.columns),
        "obs_names_unique": bool(adata.obs_names.is_unique),
        "var_names_unique": bool(adata.var_names.is_unique),
        "condition_counts": condition_counts.to_dict(orient="records"),
        "sample_summaries": sample_summaries,
        "var_name_handling": var_handling,
        "storage_verification": storage,
        "peak_rss_gb": bytes_to_gb(audit.peak_rss_bytes),
        "memory_log": str(memory_path),
    }
    write_json(output_dir / "audit" / "reconstruction_summary.json", summary)

    del adata
    del X
    del sample_matrices
    return summary


def parse_args() -> argparse.Namespace:
    root = repo_root()
    config = load_config(root)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default=os.environ.get(config["data_dir_env"], config["default_data_dir"]),
        help="Directory containing downloaded public GEO files.",
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get(config["output_dir_env"], config["default_output_dir"]),
        help="Directory for audit artifacts.",
    )
    parser.add_argument(
        "--output-h5ad",
        default=os.environ.get(config["input_h5ad_env"], config["default_input_h5ad"]),
        help="Sparse reconstructed H5AD path.",
    )
    parser.add_argument(
        "--memory-threshold-gb",
        type=float,
        default=float(os.environ.get(config["memory_threshold_gb_env"], config["default_memory_threshold_gb"])),
        help="RSS threshold in GB.",
    )
    return parser.parse_args()


def main() -> int:
    root = repo_root()
    args = parse_args()
    summary = run(
        root=root,
        data_dir=resolve_path(root, args.data_dir),
        output_dir=resolve_path(root, args.output_dir),
        output_h5ad=resolve_path(root, args.output_h5ad),
        memory_threshold_gb=float(args.memory_threshold_gb),
    )
    json.dump(summary, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
