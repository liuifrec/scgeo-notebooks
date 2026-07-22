#!/usr/bin/env python3
"""Execute repaired GSE211713 C7 v2 notebooks in isolated processes."""

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

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "results/public_validation/gse211713_dataset_c"
OUTPUT = RESULT_ROOT / "scgeo_v2"
EXECUTED = RESULT_ROOT / "executed_notebooks_v2"
CONFIG = ROOT / "configs/gse211713_scgeo_v1.json"
NOTEBOOKS = [
    ("notebooks/public_validation/gse211713/05a_scgeo_control_vs_17gy_early.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse211713/05b_scgeo_control_vs_17gy_late.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse211713/05c_scgeo_17gy_early_vs_late.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse211713/05d_secondary_and_fibroblast_analysis.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse211713/05e_cross_contrast_summary.ipynb", "scgeo_pre"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def output_bytes(notebook: nbformat.NotebookNode) -> int:
    return sum(len(str(output)) for cell in notebook.cells for output in cell.get("outputs", []))


def execute_one(notebook_rel: str, kernel: str, timeout: int) -> int:
    source = ROOT / notebook_rel
    executed = EXECUTED / notebook_rel
    report_path = OUTPUT / "execution" / f"{source.stem}_execution_report.json"
    before = sha256(source)
    started = time.perf_counter()
    try:
        notebook = nbformat.read(source, as_version=4)
        if output_bytes(notebook) or any(cell.get("execution_count") is not None for cell in notebook.cells):
            raise RuntimeError("Source notebook is not output-free and clean-kernel")
        NotebookClient(notebook, timeout=timeout, kernel_name=kernel, resources={"metadata": {"path": str(ROOT)}}).execute()
        if sha256(source) != before:
            raise RuntimeError("Source notebook changed during execution")
        executed.parent.mkdir(parents=True, exist_ok=True)
        temporary = executed.with_name(f".{executed.name}.tmp.{os.getpid()}")
        nbformat.write(notebook, temporary)
        os.replace(temporary, executed)
        report = {
            "status": "passed", "notebook": notebook_rel, "kernel": kernel,
            "runtime_seconds": time.perf_counter() - started, "source_sha256": before,
            "source_output_bytes": 0, "executed_notebook": str(executed.relative_to(ROOT)),
            "executed_notebook_sha256": sha256(executed), "worker_pid": os.getpid(),
        }
        code = 0
    except Exception as exc:
        report = {
            "status": "failed", "notebook": notebook_rel, "kernel": kernel,
            "runtime_seconds": time.perf_counter() - started, "source_sha256": before,
            "error": repr(exc), "traceback": traceback.format_exc(), "worker_pid": os.getpid(),
        }
        code = 1
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, report_path)
    print(json.dumps(report, indent=2), flush=True)
    return code


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker-notebook")
    parser.add_argument("--kernel")
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--resume-secondary-and-summary-after-module-invocation-fix", action="store_true")
    parser.add_argument("--finalize-report-only", action="store_true")
    args = parser.parse_args()
    if args.worker_notebook:
        return execute_one(args.worker_notebook, args.kernel, args.timeout)

    config = json.loads(CONFIG.read_text())
    resume = args.resume_secondary_and_summary_after_module_invocation_fix
    finalize_only = args.finalize_report_only
    if finalize_only:
        notebooks_to_run = []
        for notebook, _ in NOTEBOOKS:
            report_path = OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json"
            if not report_path.is_file() or json.loads(report_path.read_text()).get("status") != "passed":
                raise RuntimeError("Cannot finalize C7 report: a successful notebook report is absent")
    elif resume:
        required_primary = [
            OUTPUT / "per_contrast/control_vs_17gy_early__major/report.json",
            OUTPUT / "per_contrast/control_vs_17gy_late__major/report.json",
            OUTPUT / "per_contrast/17gy_early_vs_late__major/report.json",
        ]
        if any(not path.is_file() or json.loads(path.read_text()).get("status") != "passed" for path in required_primary):
            raise RuntimeError("Cannot resume C7: successful primary reports are absent")
        notebooks_to_run = NOTEBOOKS[-2:]
    else:
        if OUTPUT.exists() and any(OUTPUT.iterdir()):
            raise RuntimeError(f"Canonical C7 v2 output directory must start empty: {OUTPUT}")
        OUTPUT.mkdir(parents=True, exist_ok=True)
        notebooks_to_run = NOTEBOOKS
    if subprocess.run(["git", "-C", str(ROOT), "check-ignore", "-q", str(EXECUTED.relative_to(ROOT))]).returncode:
        raise RuntimeError("Executed-notebook v2 directory is not ignored")
    env = os.environ.copy()
    env.update({
        "SCGEO_GSE211713_REPRESENTATION_V2_H5AD": config["default_representation_h5ad"],
        "SCGEO_GSE211713_C7_V2_OUTPUT_DIR": str(OUTPUT),
        "SCGEO_SOURCE_REPO": "/home/liuyuchen/Github/scgeo",
        "NUMBA_CACHE_DIR": str(OUTPUT / "_numba_cache"),
        "MPLCONFIGDIR": str(OUTPUT / "_matplotlib_cache"),
        "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    })
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    rows = []
    started = time.perf_counter()
    status = "passed"
    for notebook, kernel in notebooks_to_run:
        print(f"[gse211713-c7-v2] fresh worker: {notebook}", flush=True)
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker-notebook", notebook, "--kernel", kernel, "--timeout", str(args.timeout)],
            cwd=ROOT, env=env, check=False,
        )
        report_path = OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json"
        rows.append(json.loads(report_path.read_text()))
        if result.returncode:
            status = "failed"
            break
    if resume or finalize_only:
        rows = [
            json.loads((OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json").read_text())
            for notebook, _ in NOTEBOOKS
        ]
        status = "passed" if all(row.get("status") == "passed" for row in rows) else "failed"
    per_contrast_reports = {}
    for path in sorted((OUTPUT / "per_contrast").glob("*/report.json")):
        per_contrast_reports[str(path.relative_to(OUTPUT))] = json.loads(path.read_text())
    aggregate = {
        "workflow": "gse211713_phase_c7_v2_replicate_aware_scgeo", "status": status,
        "fresh_process_per_heavy_stage": True,
        "runtime_seconds": float(sum(row["runtime_seconds"] for row in rows)),
        "report_finalization_runtime_seconds": time.perf_counter() - started,
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "repository_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "repository_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True).splitlines(),
        "frozen_scgeo_commit": subprocess.check_output(["git", "-C", "/home/liuyuchen/Github/scgeo", "rev-parse", "HEAD"], text=True).strip(),
        "frozen_scgeo_clean": not bool(subprocess.check_output(["git", "-C", "/home/liuyuchen/Github/scgeo", "status", "--short"], text=True).strip()),
        "input_checksums": {
            "compact": sha256(Path(config["default_compact_h5ad"])),
            "annotated": sha256(Path(config["default_annotated_h5ad"])),
            "representation_v2": sha256(Path(config["default_representation_h5ad"])),
        },
        "configuration_sha256": sha256(CONFIG), "notebooks": rows,
        "per_contrast_reports": per_contrast_reports,
        "peak_rss_gib": max([float(x.get("peak_rss_gib", 0)) for x in per_contrast_reports.values()] or [0]),
        "peak_gpu_allocated_gib": 0.0, "peak_gpu_reserved_gib": 0.0,
        "warnings": [
            "Exact mouse-label fractions have discrete small-sample resolution; related contrasts reuse control mice.",
            *(["Secondary/fibroblast jobs and summary resumed after a module-invocation-only wrapper defect; three primary analyses were not rerun."] if (resume or finalize_only) else []),
        ],
        "deviations_from_approved_protocol": [],
        "invalidated_run_manifest": config["invalidation_manifest"],
        "repair_attempt_manifests": [
            "results/public_validation/gse211713_dataset_c/invalidated_runs/c6_c7_repair_attempt_1_notebook_escape/invalidation_manifest.json",
            "results/public_validation/gse211713_dataset_c/invalidated_runs/c6_c7_repair_attempt_2_obs_names_object/invalidation_manifest.json",
            "results/public_validation/gse211713_dataset_c/invalidated_runs/c6_c7_repair_attempt_3_c7_module_invocation/invalidation_manifest.json"
        ],
        "resumed_secondary_and_summary_only": bool(resume or finalize_only),
        "report_finalized_without_computation": finalize_only,
    }
    if status == "passed":
        aggregate["output_checksums"] = {
            str(path.relative_to(OUTPUT)): sha256(path)
            for path in sorted(OUTPUT.rglob("*"))
            if path.is_file() and path.name != "gse211713_phase_c7_v2_execution_report.json"
        }
    report_path = OUTPUT / "execution/gse211713_phase_c7_v2_execution_report.json"
    temporary = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, report_path)
    print(json.dumps({"status": status, "aggregate_report": str(report_path)}, indent=2), flush=True)
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
