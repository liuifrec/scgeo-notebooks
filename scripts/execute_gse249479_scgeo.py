#!/usr/bin/env python3
"""Execute GSE249479 Phase 3B notebooks in fresh worker and kernel processes."""

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
    "notebooks/public_validation/gse249479/04a_scgeo_pbs_tnf.ipynb",
    "notebooks/public_validation/gse249479/04b_scgeo_pbs_lps.ipynb",
    "notebooks/public_validation/gse249479/04c_cross_condition_state_summary.ipynb",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def output_bytes(nb: nbformat.NotebookNode) -> int:
    total = 0
    for cell in nb.cells:
        for output in cell.get("outputs", []):
            if "text" in output:
                value = output["text"]; total += len("".join(value) if isinstance(value, list) else str(value))
            for value in output.get("data", {}).values():
                total += len("".join(value) if isinstance(value, list) else str(value))
    return total


def worker(notebook_rel: str, timeout: int) -> int:
    source = ROOT / notebook_rel
    executed = OUTPUT / "executed_notebooks" / notebook_rel
    report_path = OUTPUT / "execution" / f"{source.stem}_execution_report.json"
    before, started = sha256(source), time.perf_counter()
    try:
        nb = nbformat.read(source, as_version=4)
        embedded = output_bytes(nb)
        if embedded or any(c.get("execution_count") is not None for c in nb.cells):
            raise RuntimeError("Source notebook must be output-free and clean-kernel")
        NotebookClient(nb, timeout=timeout, kernel_name="scgeo_pre", resources={"metadata": {"path": str(ROOT)}}).execute()
        after = sha256(source)
        if after != before: raise RuntimeError("Source notebook changed during execution")
        executed.parent.mkdir(parents=True, exist_ok=True)
        tmp = executed.with_name(f".{executed.name}.tmp.{os.getpid()}"); nbformat.write(nb, tmp); os.replace(tmp, executed)
        row = {"notebook": notebook_rel, "kernel": "scgeo_pre", "status": "passed", "inference_status": "descriptive_only", "runtime_seconds": time.perf_counter()-started,
               "source_sha256_before": before, "source_sha256_after": after, "source_output_bytes": embedded,
               "executed_notebook": str(executed.relative_to(ROOT)), "executed_notebook_sha256": sha256(executed),
               "executed_output_bytes": output_bytes(nb), "worker_pid": os.getpid()}
        code = 0
    except Exception as exc:
        row = {"notebook": notebook_rel, "kernel": "scgeo_pre", "status": "failed", "inference_status": "descriptive_only", "runtime_seconds": time.perf_counter()-started,
               "source_sha256_before": before, "source_sha256_after": sha256(source), "error": repr(exc),
               "traceback": traceback.format_exc(), "worker_pid": os.getpid()}; code = 1
    report_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}"); tmp.write_text(json.dumps(row, indent=2, sort_keys=True)+"\n"); os.replace(tmp, report_path)
    print(json.dumps(row, indent=2), flush=True); return code


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--worker-notebook"); parser.add_argument("--timeout",type=int,default=14400); args=parser.parse_args()
    if args.worker_notebook: return worker(args.worker_notebook,args.timeout)
    OUTPUT.mkdir(parents=True,exist_ok=True)
    if subprocess.run(["git","-C",str(ROOT),"check-ignore","-q",str((OUTPUT/"executed_notebooks").relative_to(ROOT))]).returncode:
        raise RuntimeError("Executed-notebook directory must be ignored")
    env=os.environ.copy(); env.update({
        "SCGEO_GSE249479_REPRESENTATION_H5AD":"/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_representations.h5ad",
        "SCGEO_GSE249479_COMPACT_H5AD":"/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_qc_hvg.h5ad",
        "SCGEO_GSE249479_OUTPUT_DIR":str(OUTPUT), "SCGEO_SOURCE_REPO":"/home/liuyuchen/Github/scgeo",
        "PYTHONPATH":"/home/liuyuchen/Github/scgeo:"+str(ROOT/"scripts"), "NUMBA_CACHE_DIR":str(OUTPUT/"_numba_cache_scgeo"),
        "MPLCONFIGDIR":str(OUTPUT/"_matplotlib_cache_scgeo"), "PYTHONHASHSEED":"0", "OMP_NUM_THREADS":"1", "MKL_NUM_THREADS":"1", "OPENBLAS_NUM_THREADS":"1"})
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True,exist_ok=True); Path(env["MPLCONFIGDIR"]).mkdir(parents=True,exist_ok=True)
    rows=[]; overall=0
    for notebook in NOTEBOOKS:
        print(f"[gse249479-phase3b] fresh worker: {notebook}",flush=True)
        result=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--worker-notebook",notebook,"--timeout",str(args.timeout)],cwd=ROOT,env=env,check=False)
        report=OUTPUT/"execution"/f"{Path(notebook).stem}_execution_report.json"; rows.append(json.loads(report.read_text()))
        if result.returncode: overall=result.returncode; break
    final={"workflow":"gse249479_phase3b_descriptive_scgeo","status":"passed" if overall==0 else "failed","inference_status":"descriptive_only","fresh_worker_per_notebook":True,"notebooks":rows}
    (OUTPUT/"execution"/"gse249479_phase3b_execution_report.json").write_text(json.dumps(final,indent=2,sort_keys=True)+"\n")
    return overall


if __name__ == "__main__": raise SystemExit(main())
