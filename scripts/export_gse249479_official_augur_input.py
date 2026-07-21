"""Export the immutable sparse 3,000-HVG input for official R Augur."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=ROOT / "configs/gse249479_official_augur_validation_v1.json", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    source = Path(config["compact_h5ad"])
    if sha256(source) != config["compact_sha256"]:
        raise RuntimeError("Compact H5AD checksum mismatch")
    output = ROOT / config["validation_output"] / "input"
    output.mkdir(parents=True, exist_ok=True)
    targets = [output / "expression_genes_by_cells.mtx.gz", output / "cells.csv", output / "genes.csv"]
    if any(path.exists() for path in targets):
        raise FileExistsError("Official Augur export already exists; refusing to overwrite")

    backed = ad.read_h5ad(source, backed="r")
    hvg = backed.var["highly_variable"].to_numpy(bool)
    if int(hvg.sum()) != int(config["n_hvg"]):
        raise RuntimeError("Expected exactly 3,000 HVGs")
    block = backed[:, hvg].X
    if hasattr(block, "to_memory"):
        block = block.to_memory()
    matrix = block.tocsr().astype(np.float32)
    obs = backed.obs.copy()
    genes = backed.var_names[hvg].astype(str).to_numpy()
    backed.file.close()

    totals = obs["total_counts"].to_numpy(np.float32)
    factors = np.divide(np.float32(1e4), totals, out=np.zeros_like(totals), where=totals > 0)
    matrix = matrix.multiply(factors[:, None]).tocsr().astype(np.float32)
    np.log1p(matrix.data, out=matrix.data)
    if not sparse.isspmatrix_csr(matrix) or matrix.dtype != np.float32:
        raise RuntimeError("Export matrix is not CSR float32")

    mapped = obs["marker_inferred_label"].astype(str).map(config["detailed_state_mapping"])
    if mapped.isna().any():
        raise RuntimeError("Unmapped detailed state label")
    cells = pd.DataFrame({
        "cell_id": obs.index.astype(str),
        "condition": obs["condition"].astype(str).to_numpy(),
        "state_detailed": mapped.to_numpy(),
        "marker_inferred_label": obs["marker_inferred_label"].astype(str).to_numpy(),
        "inference_status": "descriptive_only",
    })
    gene_frame = pd.DataFrame({"gene_id": genes, "highly_variable": True})

    mtx_tmp = output / "expression_genes_by_cells.mtx.gz.tmp"
    with gzip.open(mtx_tmp, "wb", compresslevel=6) as handle:
        mmwrite(handle, matrix.transpose().tocoo(), field="real", precision=8)
    os.replace(mtx_tmp, targets[0])
    cells.to_csv(targets[1], index=False)
    gene_frame.to_csv(targets[2], index=False)

    manifest = {
        "inference_status": "descriptive_only",
        "source_h5ad": str(source),
        "source_sha256": config["compact_sha256"],
        "shape_cells_by_genes": [int(matrix.shape[0]), int(matrix.shape[1])],
        "sparse_source_encoding": "csr_float32",
        "export_encoding": "MatrixMarket_coordinate_genes_by_cells_gzip",
        "nnz": int(matrix.nnz),
        "normalization": config["normalization"],
        "files": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size} for path in targets},
    }
    atomic_text(output / "input_manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
