#!/usr/bin/env python3
"""Execute the GSE211713 Phase C0 audit notebook in a clean kernel."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path
from typing import Any

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/gse211713_dataset_c_v1.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def notebook_output_bytes(nb: nbformat.NotebookNode) -> int:
    total = 0
    for cell in nb.cells:
        for output in cell.get("outputs", []):
            total += len(str(output))
    return total


def ensure_ignored(path: Path) -> None:
    relative = path.resolve().relative_to(ROOT)
    result = subprocess.run(["git", "-C", str(ROOT), "check-ignore", "-q", str(relative)], check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Executed notebook path must be ignored by Git: {relative}")


def execute(kernel_name: str, timeout: int) -> dict[str, Any]:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    notebook = ROOT / "notebooks/public_validation/gse211713/00_study_design_replication_and_file_audit.ipynb"
    output_root = (ROOT / config["output_dir"]).resolve()
    executed = (ROOT / config["executed_notebook_dir"] / notebook.relative_to(ROOT)).resolve()
    ensure_ignored(executed)

    source_hash_before = sha256_file(notebook)
    nb = nbformat.read(notebook, as_version=4)
    if notebook_output_bytes(nb):
        raise RuntimeError("Source notebook contains embedded outputs")
    if any(cell.get("execution_count") is not None for cell in nb.cells if cell.cell_type == "code"):
        raise RuntimeError("Source notebook contains execution counts")

    cache = output_root / "_jupyter_cache"
    cache.mkdir(parents=True, exist_ok=True)
    old_environment = {key: os.environ.get(key) for key in ["IPYTHONDIR", "JUPYTER_RUNTIME_DIR", "MPLCONFIGDIR", "PYTHONHASHSEED"]}
    os.environ.update({
        "IPYTHONDIR": str(cache / "ipython"),
        "JUPYTER_RUNTIME_DIR": str(cache / "runtime"),
        "MPLCONFIGDIR": str(cache / "matplotlib"),
        "PYTHONHASHSEED": "0",
    })
    for key in ["IPYTHONDIR", "JUPYTER_RUNTIME_DIR", "MPLCONFIGDIR"]:
        Path(os.environ[key]).mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    try:
        NotebookClient(
            nb, timeout=timeout, kernel_name=kernel_name,
            resources={"metadata": {"path": str(ROOT)}},
        ).execute()
    finally:
        for key, value in old_environment.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    runtime = time.perf_counter() - started

    if sha256_file(notebook) != source_hash_before:
        raise RuntimeError("Source notebook changed during clean-kernel execution")
    executed.parent.mkdir(parents=True, exist_ok=True)
    temporary = executed.with_name(f".{executed.name}.tmp.{os.getpid()}")
    nbformat.write(nb, temporary)
    os.replace(temporary, executed)

    report = {
        "status": "passed", "notebook": str(notebook.relative_to(ROOT)),
        "source_sha256": source_hash_before, "source_outputs": 0,
        "executed_notebook": str(executed.relative_to(ROOT)),
        "executed_sha256": sha256_file(executed), "kernel": kernel_name,
        "runtime_seconds": runtime,
    }
    report_path = output_root / "audit/execution_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_report = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}")
    temporary_report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary_report, report_path)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel-name", default="scgeo_pre")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    print(json.dumps(execute(args.kernel_name, args.timeout), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
