#!/usr/bin/env python3
"""Build per-GSM and combined sparse H5AD files for GSE211713."""

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
import numpy as np
import pandas as pd
import psutil
from scipy import io, sparse


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "configs/gse211713_dataset_c_v1.json").read_text())
DATA_DIR = Path(os.environ.get("SCGEO_GSE211713_DATA_DIR", "/home/liuyuchen/data/gse211713")).resolve()
WARNING_BYTES = 20 * 1024**3
STOP_BYTES = 24 * 1024**3


def rss_bytes() -> int:
    return psutil.Process().memory_info().rss


def guard(stage: str, records: list[dict[str, Any]]) -> None:
    rss = rss_bytes()
    records.append({"stage": stage, "rss_bytes": rss, "rss_gib": rss / 1024**3, "timestamp_utc": datetime.now(timezone.utc).isoformat()})
    if rss >= STOP_BYTES:
        raise MemoryError(f"Hard stop at {rss / 1024**3:.2f} GiB during {stage}")
    if rss >= WARNING_BYTES:
        print(f"WARNING: RSS {rss / 1024**3:.2f} GiB during {stage}", file=sys.stderr)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_record(gsm: str) -> dict[str, Any]:
    for values in CONFIG["samples"]:
        record = dict(zip(CONFIG["sample_columns"], values))
        if record["geo_accession"] == gsm:
            return record
    raise KeyError(gsm)


def paths_for(record: dict[str, Any]) -> dict[str, Path]:
    gsm, title = record["geo_accession"], record["sample_title"]
    directory = DATA_DIR / "processed_mex" / gsm
    stem = f"{gsm}_{title}"
    return {
        "barcodes": directory / f"{stem}_barcodes.tsv.gz",
        "features": directory / f"{stem}_genes.tsv.gz",
        "matrix": directory / f"{stem}_count_matrix.mtx.gz",
        "h5ad": DATA_DIR / "sample_h5ad" / f"{gsm}.h5ad",
    }


def primary_group(dose: int, month: int | None) -> str:
    if dose == 0:
        return "control"
    if dose == 17 and month in (1, 2):
        return "17Gy_early"
    if dose == 17 and month in (4, 5):
        return "17Gy_late"
    return "not_primary"


def write_atomic(adata: ad.AnnData, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.h5ad")
    adata.write_h5ad(temporary, compression="gzip", compression_opts=4)
    backed = ad.read_h5ad(temporary, backed="r")
    try:
        if backed.shape != adata.shape:
            raise RuntimeError(f"Backed shape mismatch for {path}")
        if getattr(backed.X, "format", None) not in {"csr", "csc"} and "CSRDataset" not in type(backed.X).__name__ and "CSCDataset" not in type(backed.X).__name__:
            raise RuntimeError(f"Backed X is not sparse for {path}")
    finally:
        backed.file.close()
    os.replace(temporary, path)


def build_sample(gsm: str) -> dict[str, Any]:
    started = time.perf_counter(); memory: list[dict[str, Any]] = []
    record = sample_record(gsm); paths = paths_for(record)
    report_path = DATA_DIR / "manifests/samples" / f"{gsm}_reconstruction.json"
    if paths["h5ad"].is_file() and report_path.is_file():
        previous = json.loads(report_path.read_text())
        if previous.get("status") == "passed" and previous.get("sha256") == sha256_file(paths["h5ad"]):
            backed = ad.read_h5ad(paths["h5ad"], backed="r")
            try:
                if list(backed.shape) == previous.get("shape") and "CSRDataset" in type(backed.X).__name__:
                    print(json.dumps({"gsm": gsm, "shape": previous["shape"], "status": "reused_validated"}))
                    return previous
            finally:
                backed.file.close()
    guard("start", memory)
    for key in ("barcodes", "features", "matrix"):
        if not paths[key].is_file():
            raise FileNotFoundError(paths[key])
    barcodes = pd.read_csv(paths["barcodes"], sep="\t", header=None, names=["barcode"], dtype=str)
    features = pd.read_csv(paths["features"], sep="\t", header=None, names=["ensembl_id", "gene_symbol", "feature_type"], dtype=str)
    guard("metadata_loaded", memory)
    matrix = io.mmread(paths["matrix"]).T.tocsr().astype(np.float32, copy=False)
    matrix.sort_indices()
    guard("matrix_loaded_csr", memory)
    if matrix.shape != (len(barcodes), len(features)):
        raise RuntimeError(f"MEX dimension mismatch for {gsm}: {matrix.shape}")
    dose = int(record["dose_gy"]); month = record["time_month"]
    obs_names = pd.Index(gsm + "_" + barcodes["barcode"], name="cell_id")
    obs = pd.DataFrame(index=obs_names)
    obs["barcode"] = barcodes["barcode"].to_numpy()
    obs["mouse_id"] = gsm
    obs["gsm"] = gsm
    obs["dose_gy"] = dose
    obs["month_post_irradiation"] = np.nan if month is None else int(month)
    obs["irradiation_group"] = {0: "control", 10: "10Gy", 17: "17Gy"}[dose]
    obs["primary_contrast_group"] = primary_group(dose, month)
    obs["biological_replicate_eligible"] = True
    var = features.set_index("ensembl_id")
    var.index.name = "ensembl_id"
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.uns["source"] = {
        "geo_accession": gsm, "series": "GSE211713", "assembly": "mm10",
        "matrix_sha256": sha256_file(paths["matrix"]),
        "barcodes_sha256": sha256_file(paths["barcodes"]),
        "features_sha256": sha256_file(paths["features"]),
        "inference_unit": "mouse/GSM library",
    }
    write_atomic(adata, paths["h5ad"])
    guard("written_and_backed_validated", memory)
    report = {
        "status": "passed", "gsm": gsm, "shape": list(adata.shape), "nnz": int(matrix.nnz),
        "dtype": str(matrix.dtype), "sparse_format": "csr", "output": str(paths["h5ad"]),
        "sha256": sha256_file(paths["h5ad"]), "runtime_seconds": time.perf_counter() - started,
        "peak_rss_gib": max(x["rss_gib"] for x in memory), "memory": memory,
    }
    temporary = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, report_path)
    print(json.dumps({k: report[k] for k in ["gsm", "shape", "nnz", "peak_rss_gib", "status"]}))
    return report


def build_all_samples() -> None:
    for values in CONFIG["samples"]:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--gsm", values[0]], check=True)


def concatenate() -> dict[str, Any]:
    started = time.perf_counter(); memory: list[dict[str, Any]] = []
    output = DATA_DIR / "gse211713_all20_raw_sparse.h5ad"
    path = DATA_DIR / "manifests/gse211713_raw_sparse_validation.json"
    if output.is_file() and path.is_file():
        previous = json.loads(path.read_text())
        if previous.get("status") == "passed" and previous.get("sha256") == sha256_file(output):
            backed = ad.read_h5ad(output, backed="r")
            try:
                if list(backed.shape) == previous.get("shape") and "CSRDataset" in type(backed.X).__name__:
                    print(json.dumps({"shape": previous["shape"], "sha256": previous["sha256"], "status": "reused_validated"}))
                    return previous
            finally:
                backed.file.close()
    objects: list[ad.AnnData] = []
    guard("concat_start", memory)
    for values in CONFIG["samples"]:
        gsm = values[0]
        obj = ad.read_h5ad(paths_for(sample_record(gsm))["h5ad"])
        if not sparse.isspmatrix_csr(obj.X):
            obj.X = sparse.csr_matrix(obj.X, dtype=np.float32)
        objects.append(obj)
        guard(f"loaded_{gsm}", memory)
    combined = ad.concat(objects, axis=0, join="inner", merge="same", uns_merge=None, index_unique=None)
    combined.X = combined.X.tocsr().astype(np.float32, copy=False)
    combined.var = objects[0].var.copy()
    combined.uns["dataset"] = {
        "accession": "GSE211713", "mouse_libraries": 20,
        "inference_unit": "mouse/GSM library", "cross_sectional": True,
        "source_sample_h5ad_sha256": {values[0]: sha256_file(paths_for(sample_record(values[0]))["h5ad"]) for values in CONFIG["samples"]},
    }
    guard("concatenated", memory)
    write_atomic(combined, output)
    guard("combined_written_validated", memory)
    backed = ad.read_h5ad(output, backed="r")
    try:
        validation = {"shape": list(backed.shape), "x_class": type(backed.X).__name__, "obs_fields": list(backed.obs.columns), "var_fields": list(backed.var.columns)}
    finally:
        backed.file.close()
    report = {
        "status": "passed", "output": str(output), "sha256": sha256_file(output),
        "shape": list(combined.shape), "nnz": int(combined.X.nnz), "dtype": str(combined.X.dtype),
        "sparse_format": "csr", "validation": validation,
        "runtime_seconds": time.perf_counter() - started, "peak_rss_gib": max(x["rss_gib"] for x in memory), "memory": memory,
    }
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)
    print(json.dumps({k: report[k] for k in ["shape", "nnz", "peak_rss_gib", "sha256", "status"]}))
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gsm")
    group.add_argument("--all-samples", action="store_true")
    group.add_argument("--concatenate", action="store_true")
    args = parser.parse_args()
    if args.gsm:
        build_sample(args.gsm)
    elif args.all_samples:
        build_all_samples()
    else:
        concatenate()
    gc.collect()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
