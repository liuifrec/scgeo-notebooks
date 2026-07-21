#!/usr/bin/env python3
"""Run GSE249479 Phase 3A notebooks in fresh worker and kernel processes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/public_validation/gse249479_dataset_b"
NOTEBOOKS = [
    ("notebooks/public_validation/gse249479/03a_pca_representations.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse249479/03b_diffusion_representation.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse249479/03c_scvi_representation.ipynb", "sc_atac"),
    ("notebooks/public_validation/gse249479/03d_representation_quality.ipynb", "scgeo_pre"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def output_bytes(nb: nbformat.NotebookNode) -> int:
    total = 0
    for cell in nb.cells:
        for output in cell.get("outputs", []):
            if "text" in output:
                value = output["text"]
                total += len("".join(value) if isinstance(value, list) else str(value))
            for value in output.get("data", {}).values():
                total += len("".join(value) if isinstance(value, list) else str(value))
    return total


def execute_worker(notebook_rel: str, kernel: str, timeout: int) -> int:
    source = ROOT / notebook_rel
    executed = OUTPUT / "executed_notebooks" / notebook_rel
    report_path = OUTPUT / "execution" / f"{source.stem}_execution_report.json"
    before = sha256(source)
    started = time.perf_counter()
    row: dict[str, Any]
    try:
        nb = nbformat.read(source, as_version=4)
        embedded = output_bytes(nb)
        if embedded:
            raise RuntimeError(f"Source notebook contains {embedded} output bytes")
        client = NotebookClient(nb, timeout=timeout, kernel_name=kernel, resources={"metadata": {"path": str(ROOT)}})
        client.execute()
        after = sha256(source)
        if after != before:
            raise RuntimeError("Source notebook changed during execution")
        executed.parent.mkdir(parents=True, exist_ok=True)
        tmp = executed.with_name(f".{executed.name}.tmp.{os.getpid()}")
        nbformat.write(nb, tmp)
        os.replace(tmp, executed)
        row = {
            "notebook": notebook_rel, "kernel": kernel, "status": "passed",
            "runtime_seconds": time.perf_counter() - started,
            "source_sha256_before": before, "source_sha256_after": after,
            "source_output_bytes": embedded, "executed_notebook": str(executed.relative_to(ROOT)),
            "executed_notebook_sha256": sha256(executed), "executed_output_bytes": output_bytes(nb),
            "worker_pid": os.getpid(),
        }
        code = 0
    except Exception as exc:
        row = {
            "notebook": notebook_rel, "kernel": kernel, "status": "failed",
            "runtime_seconds": time.perf_counter() - started, "source_sha256_before": before,
            "source_sha256_after": sha256(source), "error": repr(exc), "traceback": traceback.format_exc(),
            "worker_pid": os.getpid(),
        }
        code = 1
    report_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_report = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}")
    tmp_report.write_text(json.dumps(row, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp_report, report_path)
    print(json.dumps(row, indent=2), flush=True)
    return code


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker-notebook")
    parser.add_argument("--kernel")
    parser.add_argument("--timeout", type=int, default=7200)
    args = parser.parse_args()
    if args.worker_notebook:
        return execute_worker(args.worker_notebook, args.kernel, args.timeout)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    executed_root = OUTPUT / "executed_notebooks"
    ignored = subprocess.run(["git", "-C", str(ROOT), "check-ignore", "-q", str(executed_root.relative_to(ROOT))]).returncode == 0
    if not ignored:
        raise RuntimeError("Executed-notebook directory must be ignored by Git")

    env = os.environ.copy()
    env.update({
        "SCGEO_GSE249479_COMPACT_H5AD": "/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_qc_hvg.h5ad",
        "SCGEO_GSE249479_REPRESENTATION_H5AD": "/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_representations.h5ad",
        "SCGEO_GSE249479_OUTPUT_DIR": str(OUTPUT),
        "SCGEO_SOURCE_REPO": "/home/liuyuchen/Github/scgeo",
        "NUMBA_CACHE_DIR": str(OUTPUT / "_numba_cache_representations"),
        "MPLCONFIGDIR": str(OUTPUT / "_matplotlib_cache_representations"),
        "PYTHONHASHSEED": "0",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    })
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    rows = []
    overall = 0
    for notebook, kernel in NOTEBOOKS:
        print(f"[gse249479-phase3] fresh worker: {notebook} ({kernel})", flush=True)
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker-notebook", notebook, "--kernel", kernel, "--timeout", str(args.timeout)],
            cwd=ROOT, env=env, check=False,
        )
        report_path = OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json"
        rows.append(json.loads(report_path.read_text(encoding="utf-8")))
        if result.returncode != 0:
            overall = result.returncode
            break
    report = {
        "workflow": "gse249479_phase3a_representations", "status": "passed" if overall == 0 else "failed",
        "fresh_worker_per_notebook": True, "notebooks": rows,
    }
    path = OUTPUT / "execution" / "gse249479_phase3a_execution_report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return overall


if __name__ == "__main__":
    raise SystemExit(main())
