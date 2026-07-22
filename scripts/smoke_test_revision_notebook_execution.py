#!/usr/bin/env python3
"""Smoke test executed-copy behavior for the revision notebook runner.

The test uses the lightweight tutorial notebook so it does not rerun the frozen
benchmark suite. It compares generated figure and CSV hashes across a baseline
run and a run with --keep-executed enabled. SVG hashes are normalized for
Matplotlib-generated element ids; PNG and CSV hashes are byte-for-byte.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import nbformat


DEFAULT_NOTEBOOK = "notebooks/tutorials/01_quickstart_perturbation_report.ipynb"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def notebook_output_bytes(path: Path) -> int:
    nb = nbformat.read(path, as_version=4)
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


def normalized_svg_bytes(path: Path) -> bytes:
    text = path.read_text(encoding="utf-8")
    volatile_ids = re.findall(r"\b[mp][0-9a-f]{10,}\b", text)
    replacements: dict[str, str] = {}
    for volatile_id in volatile_ids:
        if volatile_id not in replacements:
            replacements[volatile_id] = f"{volatile_id[0]}_stable_{len(replacements)}"
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text.encode("utf-8")


def artifact_hash(path: Path) -> str:
    if path.suffix == ".svg":
        return sha256_bytes(normalized_svg_bytes(path))
    return sha256_file(path)


def generated_hashes(output_dir: Path) -> dict[str, str]:
    patterns = ["figures/*.png", "figures/*.svg", "figure_sources/*.csv"]
    hashes: dict[str, str] = {}
    for pattern in patterns:
        for path in sorted(output_dir.glob(pattern)):
            hashes[str(path.relative_to(output_dir))] = artifact_hash(path)
    return hashes


def run_runner(root: Path, args: list[str]) -> None:
    command = [sys.executable, str(root / "scripts" / "execute_revision_notebooks.py"), *args]
    subprocess.run(command, cwd=root, check=True)


def load_report(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--notebook",
        default=DEFAULT_NOTEBOOK,
        help=f"Notebook to smoke-test, relative to repo root. Defaults to {DEFAULT_NOTEBOOK}.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the temporary output directory for debugging.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = repo_root()
    notebook_rel = Path(args.notebook)
    notebook_path = root / notebook_rel
    if not notebook_path.exists():
        raise FileNotFoundError(notebook_path)

    source_hash_before = sha256_file(notebook_path)
    if notebook_output_bytes(notebook_path) != 0:
        raise AssertionError(f"Source notebook already has embedded outputs: {notebook_rel}")

    temp_root = Path(tempfile.mkdtemp(prefix="scgeo-revision-runner-smoke-"))
    try:
        output_dir = temp_root / "outputs"
        executed_dir = temp_root / "executed_notebooks"

        run_runner(
            root,
            [
                "--notebook",
                str(notebook_rel),
                "--output-dir",
                str(output_dir),
            ],
        )
        baseline_generated = generated_hashes(output_dir)
        if not baseline_generated:
            raise AssertionError("Baseline run did not generate any figures or figure-source CSVs")
        baseline_report = load_report(output_dir / "execution" / "revision_notebook_execution_report.json")
        baseline_rows = baseline_report.get("notebooks", [])
        if len(baseline_rows) != 1:
            raise AssertionError(f"Expected one baseline notebook report row, observed {len(baseline_rows)}")
        if baseline_rows[0].get("executed_notebook") is not None:
            raise AssertionError("Baseline run unexpectedly recorded an executed notebook path")
        if "executed_notebook_sha256" in baseline_rows[0]:
            raise AssertionError("Baseline run unexpectedly recorded an executed notebook checksum")
        if (output_dir / "executed_notebooks").exists():
            raise AssertionError("Baseline run unexpectedly created an executed_notebooks directory")
        if sha256_file(notebook_path) != source_hash_before:
            raise AssertionError("Source notebook changed during baseline execution")
        if notebook_output_bytes(notebook_path) != 0:
            raise AssertionError("Source notebook gained outputs during baseline execution")

        run_runner(
            root,
            [
                "--notebook",
                str(notebook_rel),
                "--output-dir",
                str(output_dir),
                "--keep-executed",
                "--executed-dir",
                str(executed_dir),
            ],
        )
        keep_generated = generated_hashes(output_dir)
        if keep_generated != baseline_generated:
            raise AssertionError(
                "Generated figure or figure-source CSV hashes changed with --keep-executed"
            )
        if sha256_file(notebook_path) != source_hash_before:
            raise AssertionError("Source notebook changed during --keep-executed execution")
        if notebook_output_bytes(notebook_path) != 0:
            raise AssertionError("Source notebook gained outputs during --keep-executed execution")

        executed_path = executed_dir / notebook_rel
        if not executed_path.exists():
            raise AssertionError(f"Executed notebook copy was not created at {executed_path}")
        executed_output_bytes = notebook_output_bytes(executed_path)
        if executed_output_bytes <= 0:
            raise AssertionError("Executed notebook copy does not contain embedded outputs")

        report_path = output_dir / "execution" / "revision_notebook_execution_report.json"
        report = load_report(report_path)
        notebook_reports = report.get("notebooks", [])
        if len(notebook_reports) != 1:
            raise AssertionError(f"Expected one notebook report row, observed {len(notebook_reports)}")
        row = notebook_reports[0]
        if row.get("executed_notebook_sha256") != sha256_file(executed_path):
            raise AssertionError("Execution report SHA-256 does not match the executed notebook copy")
        if row.get("executed_notebook_output_bytes", 0) <= 0:
            raise AssertionError("Execution report did not record executed notebook outputs")

        print(
            "smoke test passed: source notebooks stayed output-free, executed copy has outputs, "
            "source hash is unchanged, and generated figure/CSV hashes are stable"
        )
    finally:
        if args.keep_temp:
            print(f"kept temporary directory: {temp_root}")
        else:
            shutil.rmtree(temp_root, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
