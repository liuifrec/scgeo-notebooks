#!/usr/bin/env python3
"""Execute the synthetic revision notebooks in clean kernels.

The runner does not rerun the manuscript simulation suite. It passes the frozen
benchmark directory through SCGEO_BENCHMARK_DIR and lets the notebooks load the
committed manifest/checksum files.
"""

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
    config_path = root / "configs" / "manuscript_benchmark_v1.json"
    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_from_repo(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def format_report_path(root: Path, path: Path) -> str:
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
        raise RuntimeError("--executed-dir cannot be the repository root because it could overwrite source notebooks")

    if is_relative_to(executed_dir, root):
        rel = executed_dir.relative_to(root)
        result = subprocess.run(
            ["git", "-C", str(root), "check-ignore", "-q", str(rel)],
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"--executed-dir must be ignored by Git when it is inside the repository: {rel}"
            )


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
        "figures_png": "figures/*.png",
        "figures_svg": "figures/*.svg",
        "figure_sources": "figure_sources/*.csv",
        "alt_text": "alt_text/*.txt",
        "metadata": "metadata/*.json",
    }
    artifacts: dict[str, list[str]] = {}
    for label, pattern in groups.items():
        artifacts[label] = sorted(str(path.relative_to(output_dir)) for path in output_dir.glob(pattern))
    return artifacts


def execute_one(
    notebook_path: Path,
    root: Path,
    benchmark_dir: Path,
    output_dir: Path,
    executed_dir: Path,
    kernel_name: str,
    timeout: int,
    allow_source_outputs: bool,
    keep_executed: bool,
) -> dict[str, Any]:
    notebook_rel = notebook_path.relative_to(root)
    nb = nbformat.read(notebook_path, as_version=4)
    embedded_bytes = notebook_output_bytes(nb)
    if embedded_bytes and not allow_source_outputs:
        raise RuntimeError(
            f"{notebook_rel} contains {embedded_bytes} bytes of embedded outputs"
        )

    env = {
        "SCGEO_BENCHMARK_DIR": str(benchmark_dir),
        "SCGEO_REVISION_OUTPUT_DIR": str(output_dir),
        "MPLCONFIGDIR": str(output_dir / "_matplotlib_cache"),
    }

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

    executed_path = None
    executed_sha256 = None
    executed_output_bytes = None
    if keep_executed:
        executed_path = executed_dir / notebook_rel
        if executed_path.resolve() == notebook_path.resolve():
            raise RuntimeError(f"Refusing to overwrite source notebook: {notebook_rel}")
        executed_path.parent.mkdir(parents=True, exist_ok=True)
        nbformat.write(nb, executed_path)
        executed_sha256 = sha256_file(executed_path)
        executed_output_bytes = notebook_output_bytes(nb)

    result = {
        "notebook": str(notebook_rel),
        "runtime_seconds": round(runtime, 3),
        "source_output_bytes": embedded_bytes,
        "executed_notebook": format_report_path(root, executed_path) if executed_path else None,
        "status": "passed",
    }
    if keep_executed:
        result["executed_notebook_sha256"] = executed_sha256
        result["executed_notebook_output_bytes"] = executed_output_bytes
    return result


def parse_args() -> argparse.Namespace:
    root = repo_root()
    config = load_config(root)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-dir",
        default=os.environ.get(config["benchmark_dir_env"], config["default_benchmark_dir"]),
        help="Frozen synthetic benchmark directory. Defaults to SCGEO_BENCHMARK_DIR or the config default.",
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get(config["notebook_output_env"], config["default_output_dir"]),
        help="Directory for generated figures, source tables, alt text, metadata, and execution reports.",
    )
    parser.add_argument(
        "--notebook",
        action="append",
        dest="notebooks",
        help="Notebook path relative to the repository root. Can be repeated. Defaults to all revision notebooks.",
    )
    parser.add_argument("--kernel", default="python3", help="Jupyter kernel name.")
    parser.add_argument("--timeout", type=int, default=1200, help="Per-notebook timeout in seconds.")
    parser.add_argument("--tutorial-timeout", type=int, default=300, help="Timeout in seconds for tutorial notebooks.")
    parser.add_argument(
        "--keep-executed",
        action="store_true",
        help="Save executed notebook copies with embedded outputs under --executed-dir.",
    )
    parser.add_argument(
        "--save-executed",
        action="store_true",
        dest="keep_executed",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--executed-dir",
        default=None,
        help=(
            "Directory for executed notebook copies. Defaults to "
            "<output-dir>/executed_notebooks and is used only with --keep-executed."
        ),
    )
    parser.add_argument(
        "--allow-source-outputs",
        action="store_true",
        help="Allow embedded outputs in source notebooks. Off by default to keep notebook JSON small.",
    )
    parser.add_argument("--keep-going", action="store_true", help="Continue after a notebook execution failure.")
    return parser.parse_args()


def main() -> int:
    root = repo_root()
    config = load_config(root)
    args = parse_args()

    benchmark_dir = resolve_from_repo(root, args.benchmark_dir)
    output_dir = resolve_from_repo(root, args.output_dir)
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

    notebook_list = args.notebooks or config["revision_notebooks"]
    notebook_paths = [resolve_from_repo(root, item) for item in notebook_list]

    report: dict[str, Any] = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "benchmark_dir": str(benchmark_dir),
        "output_dir": str(output_dir),
        "protocol_version": config["protocol_version"],
        "expected_source_commit": config["expected_source_commit"],
        "notebooks": [],
    }
    if args.keep_executed:
        report["keep_executed"] = True
        report["executed_dir"] = format_report_path(root, executed_dir)

    overall_status = 0
    for notebook_path in notebook_paths:
        if not notebook_path.exists():
            raise FileNotFoundError(notebook_path)
        rel = notebook_path.relative_to(root)
        print(f"[revision-notebooks] executing {rel}", flush=True)
        timeout = args.tutorial_timeout if "/tutorials/" in str(rel) else args.timeout
        try:
            result = execute_one(
                notebook_path=notebook_path,
                root=root,
                benchmark_dir=benchmark_dir,
                output_dir=output_dir,
                executed_dir=executed_dir,
                kernel_name=args.kernel,
                timeout=timeout,
                allow_source_outputs=args.allow_source_outputs,
                keep_executed=args.keep_executed,
            )
            print(f"[revision-notebooks] passed {rel} in {result['runtime_seconds']}s", flush=True)
        except Exception as exc:  # noqa: BLE001 - report notebook failures with context.
            overall_status = 1
            result = {
                "notebook": str(rel),
                "status": "failed",
                "error": repr(exc),
                "traceback": traceback.format_exc(),
            }
            print(f"[revision-notebooks] failed {rel}: {exc}", file=sys.stderr, flush=True)
            if not args.keep_going:
                report["notebooks"].append(result)
                break
        report["notebooks"].append(result)

    report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    report["artifacts"] = collect_artifacts(output_dir)
    report_path = output_dir / "execution" / "revision_notebook_execution_report.json"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    print(f"[revision-notebooks] wrote {report_path}", flush=True)
    return overall_status


if __name__ == "__main__":
    raise SystemExit(main())
