#!/usr/bin/env python3
"""Execute Dataset B GSE249479 Phase 1 validation notebooks in clean kernels."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import nbformat
from nbclient import NotebookClient


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(root: Path) -> dict[str, Any]:
    with (root / "configs" / "gse249479_dataset_b_v1.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_from_repo(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def repo_dataset_data_dir(root: Path) -> Path:
    return (root.resolve() / "data" / "gse249479").resolve()


def workflow_temp_dir(root: Path, config: dict[str, Any]) -> Path:
    env_name = config.get("temp_dir_env", "SCGEO_GSE249479_TMPDIR")
    default_value = config.get("default_temp_dir", ".tmp")
    return resolve_from_repo(root, os.environ.get(env_name, default_value))


def require_repo_local_dataset_paths(root: Path, data_dir: Path, input_h5ad: Path | None = None) -> None:
    expected_data_dir = repo_dataset_data_dir(root)
    if data_dir.resolve() != expected_data_dir:
        raise RuntimeError(f"GSE249479 data_dir must be repo-local {expected_data_dir}; observed {data_dir.resolve()}")
    if input_h5ad is not None and not is_relative_to(input_h5ad.resolve(), expected_data_dir):
        raise RuntimeError(f"GSE249479 input H5AD must live under {expected_data_dir}; observed {input_h5ad.resolve()}")


def format_report_path(root: Path, path: Path | None) -> str | None:
    if path is None:
        return None
    path = path.resolve()
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        return False
    return True


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def active_git_branch(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "branch", "--show-current"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def notebook_output_bytes(nb: nbformat.NotebookNode) -> int:
    total = 0
    for cell in nb.get("cells", []):
        for output in cell.get("outputs", []):
            if "text" in output:
                text = output["text"]
                total += len("".join(text) if isinstance(text, list) else str(text))
            data = output.get("data", {})
            if isinstance(data, dict):
                for value in data.values():
                    total += len("".join(value) if isinstance(value, list) else str(value))
    return total


def ensure_executed_dir_is_safe(root: Path, executed_dir: Path) -> None:
    executed_dir = executed_dir.resolve()
    if executed_dir == root:
        raise RuntimeError("--executed-dir cannot be the repository root")
    if is_relative_to(executed_dir, root):
        rel = executed_dir.relative_to(root)
        result = subprocess.run(["git", "-C", str(root), "check-ignore", "-q", str(rel)], check=False)
        if result.returncode != 0:
            raise RuntimeError(f"--executed-dir must be ignored by Git when inside the repository: {rel}")


@contextmanager
def patched_env(updates: dict[str, str]):
    old_values = {key: os.environ.get(key) for key in updates}
    try:
        os.environ.update(updates)
        yield
    finally:
        for key, value in old_values.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def collect_artifacts(output_dir: Path) -> dict[str, list[str]]:
    groups = {
        "audit_csv": "audit/*.csv",
        "metadata": "metadata/*.json",
        "version_records": "version_records/*.json",
    }
    artifacts: dict[str, list[str]] = {}
    for label, pattern in groups.items():
        artifacts[label] = sorted(str(path.relative_to(output_dir)) for path in output_dir.glob(pattern))
    return artifacts


def execute_one(
    notebook_path: Path,
    *,
    root: Path,
    data_dir: Path,
    input_h5ad: Path,
    output_dir: Path,
    source_repo: Path,
    memory_threshold_gb: float,
    executed_dir: Path,
    kernel_name: str,
    timeout: int,
    allow_source_outputs: bool,
    keep_executed: bool,
) -> dict[str, Any]:
    notebook_rel = notebook_path.relative_to(root)
    source_sha_before = sha256_file(notebook_path)
    nb = nbformat.read(notebook_path, as_version=4)
    embedded_bytes = notebook_output_bytes(nb)
    if embedded_bytes and not allow_source_outputs:
        raise RuntimeError(f"{notebook_rel} contains {embedded_bytes} bytes of embedded outputs")

    temp_dir = workflow_temp_dir(root, load_config(root))
    temp_dir.mkdir(parents=True, exist_ok=True)
    env = {
        "SCGEO_GSE249479_DATA_DIR": str(data_dir),
        "SCGEO_GSE249479_H5AD": str(input_h5ad),
        "SCGEO_GSE249479_OUTPUT_DIR": str(output_dir),
        "SCGEO_SOURCE_REPO": str(source_repo),
        "SCGEO_GSE249479_MEMORY_THRESHOLD_GB": str(memory_threshold_gb),
        "NUMBA_CACHE_DIR": str(output_dir / "_numba_cache"),
        "MPLCONFIGDIR": str(output_dir / "_matplotlib_cache"),
        "TMPDIR": str(temp_dir),
        "TEMP": str(temp_dir),
        "TMP": str(temp_dir),
        "PYTHONHASHSEED": "0",
    }
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["TMPDIR"]).mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    with patched_env(env):
        client = NotebookClient(
            nb,
            timeout=timeout,
            kernel_name=kernel_name,
            resources={"metadata": {"path": str(root)}},
        )
        client.execute()
    runtime = time.perf_counter() - started

    source_sha_after = sha256_file(notebook_path)
    if source_sha_after != source_sha_before:
        raise RuntimeError(f"Source notebook changed during execution: {notebook_rel}")

    executed_path = None
    executed_sha256 = None
    executed_output_bytes = None
    if keep_executed:
        executed_path = executed_dir / notebook_rel
        if executed_path.resolve() == notebook_path.resolve():
            raise RuntimeError(f"Refusing to overwrite source notebook: {notebook_rel}")
        executed_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_executed_path = executed_path.with_name(f".{executed_path.name}.tmp.{os.getpid()}")
        nbformat.write(nb, tmp_executed_path)
        os.replace(tmp_executed_path, executed_path)
        executed_sha256 = sha256_file(executed_path)
        executed_output_bytes = notebook_output_bytes(nb)

    row: dict[str, Any] = {
        "notebook": str(notebook_rel),
        "status": "passed",
        "runtime_seconds": round(runtime, 3),
        "source_sha256_before": source_sha_before,
        "source_sha256_after": source_sha_after,
        "source_output_bytes": embedded_bytes,
        "executed_notebook": format_report_path(root, executed_path),
    }
    if keep_executed:
        row["executed_notebook_sha256"] = executed_sha256
        row["executed_notebook_output_bytes"] = executed_output_bytes
    return row


def parse_args() -> argparse.Namespace:
    root = repo_root()
    config = load_config(root)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default=os.environ.get(config["data_dir_env"], config["default_data_dir"]),
        help="Directory for immutable public GEO downloads and reconstructed H5AD.",
    )
    parser.add_argument(
        "--input-h5ad",
        default=os.environ.get(config["input_h5ad_env"], config["default_input_h5ad"]),
        help="GSE249479 raw H5AD object to audit.",
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get(config["output_dir_env"], config["default_output_dir"]),
        help="Directory for audit tables, metadata, and execution report.",
    )
    parser.add_argument(
        "--source-repo",
        default=os.environ.get(config["source_repository_env"], config["default_source_repository"]),
        help="Read-only ScGeo source repository path used only for metadata.",
    )
    parser.add_argument(
        "--memory-threshold-gb",
        type=float,
        default=float(os.environ.get(config["memory_threshold_gb_env"], config["default_memory_threshold_gb"])),
        help="RSS threshold in GB. The notebook stops cleanly above this limit.",
    )
    parser.add_argument(
        "--notebook",
        action="append",
        dest="notebooks",
        help="Notebook path relative to repository root. Defaults to Phase 1 notebook 00 only.",
    )
    parser.add_argument("--kernel", default="python3", help="Jupyter kernel name.")
    parser.add_argument("--timeout", type=int, default=1200, help="Per-notebook timeout in seconds.")
    parser.add_argument(
        "--keep-executed",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Save executed review notebook copies with embedded outputs. Enabled by default.",
    )
    parser.add_argument(
        "--executed-dir",
        default=None,
        help="Directory for executed notebook copies. Defaults to <output-dir>/executed_notebooks.",
    )
    parser.add_argument("--allow-source-outputs", action="store_true", help="Allow embedded outputs in source notebooks.")
    parser.add_argument("--keep-going", action="store_true", help="Continue after notebook execution failure.")
    return parser.parse_args()


def main() -> int:
    root = repo_root()
    config = load_config(root)
    args = parse_args()
    branch = active_git_branch(root)
    required_branch = config["required_git_branch"]
    if branch != required_branch:
        print(
            f"[gse249479-validation] refusing to run: active branch is {branch!r}, "
            f"expected {required_branch!r}",
            file=sys.stderr,
            flush=True,
        )
        return 2

    data_dir = resolve_from_repo(root, args.data_dir)
    input_h5ad = resolve_from_repo(root, args.input_h5ad)
    output_dir = resolve_from_repo(root, args.output_dir)
    source_repo = resolve_from_repo(root, args.source_repo)
    require_repo_local_dataset_paths(root, data_dir, input_h5ad)
    executed_dir = (
        resolve_from_repo(root, args.executed_dir)
        if args.executed_dir
        else (output_dir / "executed_notebooks").resolve()
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "execution").mkdir(parents=True, exist_ok=True)
    if args.keep_executed:
        ensure_executed_dir_is_safe(root, executed_dir)
        executed_dir.mkdir(parents=True, exist_ok=True)

    notebook_list = args.notebooks or config["phase1_notebooks"]
    notebook_paths = [resolve_from_repo(root, item) for item in notebook_list]

    report: dict[str, Any] = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "workflow": config["workflow_name"],
        "schema_version": config["schema_version"],
        "active_branch": branch,
        "required_branch": required_branch,
        "data_dir": str(data_dir),
        "input_h5ad": str(input_h5ad),
        "input_h5ad_exists": input_h5ad.exists(),
        "output_dir": str(output_dir),
        "source_repo": str(source_repo),
        "memory_threshold_gb": float(args.memory_threshold_gb),
        "keep_executed": bool(args.keep_executed),
        "executed_dir": format_report_path(root, executed_dir) if args.keep_executed else None,
        "notebooks": [],
    }

    overall_status = 0
    for notebook_path in notebook_paths:
        if not notebook_path.exists():
            raise FileNotFoundError(notebook_path)
        rel = notebook_path.relative_to(root)
        print(f"[gse249479-validation] executing {rel}", flush=True)
        try:
            row = execute_one(
                notebook_path,
                root=root,
                data_dir=data_dir,
                input_h5ad=input_h5ad,
                output_dir=output_dir,
                source_repo=source_repo,
                memory_threshold_gb=float(args.memory_threshold_gb),
                executed_dir=executed_dir,
                kernel_name=args.kernel,
                timeout=args.timeout,
                allow_source_outputs=args.allow_source_outputs,
                keep_executed=bool(args.keep_executed),
            )
            print(f"[gse249479-validation] passed {rel} in {row['runtime_seconds']}s", flush=True)
        except Exception as exc:  # noqa: BLE001 - write clean notebook failure report.
            overall_status = 1
            row = {
                "notebook": str(rel),
                "status": "failed",
                "error": repr(exc),
                "traceback": traceback.format_exc(),
            }
            print(f"[gse249479-validation] failed {rel}: {exc}", file=sys.stderr, flush=True)
            if not args.keep_going:
                report["notebooks"].append(row)
                break
        report["notebooks"].append(row)

    report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    report["artifacts"] = collect_artifacts(output_dir)
    report_path = output_dir / "execution" / "gse249479_validation_execution_report.json"
    tmp_report_path = report_path.with_name(f".{report_path.name}.tmp.{os.getpid()}")
    with tmp_report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    os.replace(tmp_report_path, report_path)
    print(f"[gse249479-validation] wrote {report_path}", flush=True)
    return overall_status


if __name__ == "__main__":
    raise SystemExit(main())
