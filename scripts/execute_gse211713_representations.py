#!/usr/bin/env python3
"""Execute repaired GSE211713 C6 v2 notebooks in isolated fresh processes."""

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
OUTPUT = RESULT_ROOT / "c6_v2"
EXECUTED = RESULT_ROOT / "executed_notebooks_v2"
CONFIG = ROOT / "configs/gse211713_representations_v1.json"
NOTEBOOKS = [
    ("notebooks/public_validation/gse211713/04a_pca_representations.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse211713/04b_diffusion_representation.ipynb", "scgeo_pre"),
    ("notebooks/public_validation/gse211713/04c_scvi_representation.ipynb", "sc_atac"),
    ("notebooks/public_validation/gse211713/04d_representation_quality.ipynb", "scgeo_pre"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def output_bytes(nb: nbformat.NotebookNode) -> int:
    return sum(len(str(output)) for cell in nb.cells for output in cell.get("outputs", []))


def execute_one(notebook_rel: str, kernel: str, timeout: int) -> int:
    source = ROOT / notebook_rel
    executed = EXECUTED / notebook_rel
    report_path = OUTPUT / "execution" / f"{source.stem}_execution_report.json"
    before = sha256(source)
    started = time.perf_counter()
    try:
        notebook = nbformat.read(source, as_version=4)
        if output_bytes(notebook) or any(cell.get("execution_count") is not None for cell in notebook.cells):
            raise RuntimeError("Source notebook is not output-free")
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-notebook")
    parser.add_argument("--kernel")
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--resume-assembly-after-identifier-export-fix", action="store_true")
    parser.add_argument("--finalize-report-only", action="store_true")
    args = parser.parse_args()
    if args.worker_notebook:
        return execute_one(args.worker_notebook, args.kernel, args.timeout)

    config = json.loads(CONFIG.read_text())
    representation = Path(config["default_representation_h5ad"])
    resume_assembly = args.resume_assembly_after_identifier_export_fix
    finalize_only = args.finalize_report_only
    if finalize_only:
        notebooks_to_run = []
        for notebook, _ in NOTEBOOKS:
            report_path = OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json"
            if not report_path.is_file() or json.loads(report_path.read_text()).get("status") != "passed":
                raise RuntimeError("Cannot finalize C6 report: a successful notebook report is absent")
    elif resume_assembly:
        if representation.exists():
            raise RuntimeError(f"Representation v2 target already exists: {representation}")
        for name in ["pca", "diffusion", "scvi"]:
            stage = OUTPUT / f"metadata/{name}_stage.json"
            if not stage.is_file() or json.loads(stage.read_text()).get("status") != "passed":
                raise RuntimeError(f"Cannot resume: successful {name} stage is absent")
        notebooks_to_run = [NOTEBOOKS[-1]]
    else:
        if OUTPUT.exists() and any(OUTPUT.iterdir()):
            raise RuntimeError(f"C6 v2 output directory must start empty: {OUTPUT}")
        if representation.exists():
            raise RuntimeError(f"Representation v2 target already exists: {representation}")
        OUTPUT.mkdir(parents=True, exist_ok=True)
        notebooks_to_run = NOTEBOOKS
    if subprocess.run(["git", "-C", str(ROOT), "check-ignore", "-q", str(EXECUTED.relative_to(ROOT))]).returncode:
        raise RuntimeError("Executed-notebook v2 directory is not ignored")
    env = os.environ.copy()
    env.update({
        "SCGEO_GSE211713_COMPACT_H5AD": config["default_compact_h5ad"],
        "SCGEO_GSE211713_ANNOTATED_H5AD": config["default_annotated_h5ad"],
        "SCGEO_GSE211713_REPRESENTATION_V2_H5AD": config["default_representation_h5ad"],
        "SCGEO_GSE211713_C6_V2_OUTPUT_DIR": str(OUTPUT),
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
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker-notebook", notebook, "--kernel", kernel, "--timeout", str(args.timeout)],
            cwd=ROOT, env=env, check=False,
        )
        report_path = OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json"
        rows.append(json.loads(report_path.read_text()))
        if result.returncode:
            status = "failed"
            break
    if resume_assembly or finalize_only:
        rows = [
            json.loads((OUTPUT / "execution" / f"{Path(notebook).stem}_execution_report.json").read_text())
            for notebook, _ in NOTEBOOKS
        ]
        status = "passed" if all(row.get("status") == "passed" for row in rows) else "failed"
    output_checksums = {}
    if status == "passed":
        for path in sorted(OUTPUT.rglob("*")):
            if path.is_file() and path.name != "gse211713_phase_c6_v2_execution_report.json":
                output_checksums[str(path.relative_to(OUTPUT))] = sha256(path)
        output_checksums[str(representation)] = sha256(representation)
    stage_reports = {}
    for name in ["pca", "diffusion", "scvi", "assembly_quality"]:
        path = OUTPUT / f"metadata/{name}_stage.json"
        if path.exists():
            stage_reports[name] = json.loads(path.read_text())
    aggregate = {
        "workflow": "gse211713_phase_c6_v2_representations", "status": status,
        "fresh_process_per_heavy_stage": True,
        "runtime_seconds": float(sum(row["runtime_seconds"] for row in rows)),
        "report_finalization_runtime_seconds": time.perf_counter() - started,
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "repository_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "repository_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True).splitlines(),
        "frozen_scgeo_commit": subprocess.check_output(["git", "-C", "/home/liuyuchen/Github/scgeo", "rev-parse", "HEAD"], text=True).strip(),
        "frozen_scgeo_clean": not bool(subprocess.check_output(["git", "-C", "/home/liuyuchen/Github/scgeo", "status", "--short"], text=True).strip()),
        "input_checksums": {"compact": sha256(Path(config["default_compact_h5ad"])), "annotated": sha256(Path(config["default_annotated_h5ad"]))},
        "configuration_sha256": sha256(CONFIG), "output_checksums": output_checksums,
        "notebooks": rows, "stage_reports": stage_reports,
        "peak_rss_gib": max([float(x.get("peak_rss_gib", 0)) for x in stage_reports.values()] or [0]),
        "peak_gpu_allocated_gib": max([float(x.get("peak_gpu_allocated_gib", 0)) for x in stage_reports.values()] or [0]),
        "peak_gpu_reserved_gib": max([float(x.get("peak_gpu_reserved_gib", 0)) for x in stage_reports.values()] or [0]),
        "warnings": (["Assembly resumed after a non-numerical observation-name dtype export defect; PCA, diffusion, and scVI were not rerun."] if (resume_assembly or finalize_only) else []),
        "deviations_from_approved_protocol": [],
        "invalidated_run_manifest": config["invalidation_manifest"],
        "repair_attempt_manifests": [
            "results/public_validation/gse211713_dataset_c/invalidated_runs/c6_c7_repair_attempt_1_notebook_escape/invalidation_manifest.json",
            "results/public_validation/gse211713_dataset_c/invalidated_runs/c6_c7_repair_attempt_2_obs_names_object/invalidation_manifest.json",
        ],
        "resumed_assembly_only": bool(resume_assembly or finalize_only),
        "report_finalized_without_computation": finalize_only,
    }
    path = OUTPUT / "execution/gse211713_phase_c6_v2_execution_report.json"
    path.write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": status, "aggregate_report": str(path)}, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
