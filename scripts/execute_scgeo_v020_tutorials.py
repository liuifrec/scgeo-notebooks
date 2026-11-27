"""Execute independent ScGeo 0.2.0 tutorials; never edit source notebooks."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_DIR = ROOT / "notebooks" / "tutorials" / "v0.2.0"
EXAMPLES = ROOT / "scripts" / "scgeo_v020_examples"
sys.path.insert(0, str(EXAMPLES))
from tutorial_support import VERSION, RELEASE_COMMIT, verify_examples


def validate_sources():
    records = verify_examples(ROOT)
    for record in records:
        path = ROOT / record["path"]
        compile(path.read_bytes(), str(path), "exec")
    compile((EXAMPLES / "tutorial_support.py").read_bytes(), "tutorial_support.py", "exec")
    notebooks = sorted(NOTEBOOK_DIR.glob("*.ipynb"))
    if [p.name[:2] for p in notebooks] != [f"{i:02d}" for i in range(12)]:
        raise RuntimeError("Expected all 12 numbered tutorials exactly once.")
    for path in notebooks:
        book = nbformat.read(path, as_version=4)
        nbformat.validate(book)
        for i, cell in enumerate(book.cells):
            if cell.cell_type == "code":
                if cell.execution_count is not None or cell.outputs:
                    raise ValueError(f"Tracked notebook contains outputs: {path}:{i}")
                compile(cell.source, f"{path.name}:{i}", "exec")
    print(f"PASS: {len(notebooks)} clean notebook sources and {len(records)} exact release examples.")
    return notebooks


def audit_exports(out):
    """Check every actual PlotResult manifest and each file it references."""
    manifests, figures, tables = 0, 0, 0
    for path in sorted(out.rglob("*.json")):
        # Other gallery summaries can be large or contain unrelated schemas.
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, dict) or value.get("schema_version") != "scgeo.plot_export.v1":
            continue
        manifests += 1
        if value["software"]["scgeo"] != VERSION:
            raise RuntimeError(f"Incorrect plotting version: {path}")
        if set(value["figures"]) != {"png", "svg", "pdf"}:
            raise RuntimeError(f"Missing figure format: {path}")
        for kind in ("figures", "tables"):
            for record in value[kind].values():
                name = Path(record["file"])
                if name.name != str(name):
                    raise RuntimeError(f"Unexpected export path: {path}")
                export = path.parent / name
                if not export.is_file() or hashlib.sha256(export.read_bytes()).hexdigest() != record["sha256"]:
                    raise RuntimeError(f"Missing or changed exported file: {export}")
                if kind == "figures":
                    figures += 1
                else:
                    tables += 1
    return {"plot_manifests": manifests, "figures": figures, "source_tables": tables}


def execute_one(path, destination, kernel, timeout):
    source_before = path.read_bytes()
    lesson_root = destination / "lessons" / path.stem
    lesson_root.mkdir(parents=True)
    book = nbformat.read(path, as_version=4)
    environment = dict(os.environ, SCGEO_TUTORIAL_OUTPUT_ROOT=str(lesson_root),
                       SCGEO_TUTORIAL_EXPECTED_PYTHON=sys.executable, MPLBACKEND="Agg")
    started = time.monotonic()
    executed_path = destination / "executed_notebooks" / path.name
    executed_path.parent.mkdir(parents=True, exist_ok=True)
    client = NotebookClient(book, timeout=timeout, kernel_name=kernel,
                            allow_errors=False, resources={"metadata": {"path": str(ROOT)}})
    try:
        client.execute(cwd=str(ROOT), env=environment)
    finally:
        # Keep an executed/partial copy for debugging, even when a cell raises.
        nbformat.write(book, executed_path)
        if path.read_bytes() != source_before:
            raise RuntimeError(f"Source notebook changed: {path}")
    outputs = list(lesson_root.glob("*/tutorial_provenance.json"))
    if len(outputs) != 1:
        raise RuntimeError(f"Missing or ambiguous output directory for {path.name}")
    out = outputs[0].parent
    provenance = json.loads(outputs[0].read_text(encoding="utf-8"))
    if provenance["package_release"] != VERSION or provenance["release_source_commit"] != RELEASE_COMMIT:
        raise RuntimeError(f"Incorrect package provenance for {path.name}")
    stats = audit_exports(out)
    status = "executed"
    if path.name.startswith("11_"):
        public = json.loads((out / "public_data_status.json").read_text(encoding="utf-8"))
        if public["status"] == "not_run":
            status = "public_data_not_run"
        elif public["status"] != "executed":
            raise RuntimeError("Unrecognized public-data status")
    if status == "executed" and not stats["plot_manifests"]:
        raise RuntimeError(f"No actual figure manifests: {path.name}")
    record = {"notebook": path.name, "status": status, "seconds": round(time.monotonic()-started, 2),
              **stats, "output": str(out.relative_to(destination))}
    print(json.dumps(record), flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lesson", nargs="+", choices=[f"{i:02d}" for i in range(12)])
    parser.add_argument("--kernel", default="python3")
    parser.add_argument("--timeout", type=int, default=900, help="Per-cell timeout in seconds")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "scgeo_v020_tutorials")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    notebooks = validate_sources()
    if args.validate_only:
        return
    if args.lesson:
        notebooks = [path for path in notebooks if path.name[:2] in args.lesson]
    import scgeo as sg
    if sg.__version__ != VERSION:
        raise RuntimeError(f"Install ScGeo {VERSION} into the tutorial environment first.")
    base = args.output_dir.resolve()
    base.mkdir(parents=True, exist_ok=True)
    destination = Path(tempfile.mkdtemp(prefix="run-", dir=base))
    report = {"package_release": VERSION, "release_source_commit": RELEASE_COMMIT,
              "status": "running", "notebooks": []}
    try:
        for path in notebooks:
            print(f"Executing {path.name}", flush=True)
            report["notebooks"].append(execute_one(path, destination, args.kernel, args.timeout))
        report["status"] = "success"
    except BaseException:
        report["status"] = "failed"
        raise
    finally:
        (destination / "execution_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"Execution record: {destination}", flush=True)


if __name__ == "__main__":
    main()
