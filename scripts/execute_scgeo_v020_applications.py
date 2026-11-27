"""Execute public application notebooks in isolated kernels, without changing sources."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = {
    "B": "04_gse249479_inflammation.ipynb",
    "C": "05_gse211713_radiation_response.ipynb",
}


MOLECULAR_NOTEBOOKS = {
    "B": "06_gse249479_shift_explanations.ipynb",
    "C": "07_gse211713_shift_explanations.ipynb",
}


def validate():
    for name in [*NOTEBOOKS.values(), *MOLECULAR_NOTEBOOKS.values()]:
        path = ROOT / "notebooks/applications/v0.2.0" / name
        book = nbformat.read(path, as_version=4)
        nbformat.validate(book)
        for i, cell in enumerate(book.cells):
            if cell.cell_type == "code":
                if cell.outputs or cell.execution_count is not None:
                    raise ValueError(f"Tracked notebook contains outputs: {name}:{i}")
                compile(cell.source, f"{name}:{i}", "exec")
    compile((ROOT / "scripts/scgeo_v020_applications.py").read_bytes(),
            "scgeo_v020_applications.py", "exec")
    compile((ROOT / "scripts/scgeo_public_molecular.py").read_bytes(),
            "scgeo_public_molecular.py", "exec")
    print("PASS: four clean application notebook sources and helper syntax")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", choices=["B", "C", "all"], default="all")
    p.add_argument("--kernel", default="scgeo-v020-applications")
    p.add_argument("--timeout", type=int, default=1800, help="Per-cell timeout in seconds")
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--validate-only", action="store_true")
    p.add_argument("--molecular", action="store_true", help="Run only the new molecular notebooks")
    p.add_argument("--version", choices=["0.2.0", "0.2.1"],
                   default=os.environ.get("SCGEO_APPLICATION_VERSION", "0.2.0"))
    args = p.parse_args()
    validate()
    if args.validate_only:
        return
    import scgeo as sg
    if args.molecular and args.version != "0.2.1":
        p.error("The molecular track requires --version 0.2.1")
    notebooks = MOLECULAR_NOTEBOOKS if args.molecular else NOTEBOOKS
    if sg.__version__ != args.version:
        raise RuntimeError(f"Install the selected ScGeo {args.version}.")
    if args.output_dir is None:
        p.error("--output-dir is required for real applications")
    base = args.output_dir.expanduser().resolve()
    # Reject a source or any of its descendants before creating directories.
    datasets = list(notebooks) if args.dataset == "all" else [args.dataset]
    for key in ["SCGEO_APPLICATION_FROZEN_RESULTS", *[f"SCGEO_APPLICATION_{d}_ROOT" for d in datasets]]:
        if not os.environ.get(key):
            p.error(f"Set {key}")
        source = Path(os.environ[key]).expanduser().resolve()
        if not source.is_dir() or base == source or source in base.parents:
            p.error(f"Invalid or overlapping source/output directory: {key}")
    base.mkdir(parents=True, exist_ok=True)
    dest = Path(tempfile.mkdtemp(prefix="application-run-", dir=base))
    summary = {"status": "running", "package": args.version,
               "track": "molecular" if args.molecular else "geometry", "datasets": []}
    try:
        for dataset in datasets:
            source = ROOT / "notebooks/applications/v0.2.0" / notebooks[dataset]
            before = source.read_bytes()
            folder = dest / dataset
            folder.mkdir()
            book = nbformat.read(source, as_version=4)
            env = dict(os.environ, MPLBACKEND="Agg", PYTHONNOUSERSITE="1",
                       SCGEO_APPLICATION_OUTPUT_ROOT=str(folder),
                       SCGEO_APPLICATION_EXPECTED_PYTHON=sys.executable,
                       SCGEO_APPLICATION_VERSION=args.version)
            env.pop("PYTHONPATH", None)
            env.pop("PYTHONHOME", None)
            client = NotebookClient(book, timeout=args.timeout, kernel_name=args.kernel,
                                    allow_errors=False, resources={"metadata": {"path": str(ROOT)}})
            try:
                print("Executing", source.name, flush=True)
                client.execute(cwd=str(ROOT), env=env)
            finally:
                nbformat.write(book, dest / source.name)
                if source.read_bytes() != before:
                    raise RuntimeError("Tracked notebook was modified")
            manifests = list(folder.glob("*/application_manifest.json"))
            if len(manifests) != 1:
                raise RuntimeError("Expected exactly one completed dataset manifest")
            manifest = json.loads(manifests[0].read_text())
            if (manifest["status"] != "success" or manifest["dataset"] != dataset
                    or manifest["package"] != args.version
                    or Path(manifest["python_executable"]).resolve() != Path(sys.executable).resolve()):
                raise RuntimeError("Application failed or used the wrong interpreter")
            record = {"dataset": dataset, "manifest": str(manifests[0].relative_to(dest)),
                      "notebook": source.name,
                      **{k: manifest[k] for k in ("plot_manifests", "figures", "source_tables")}}
            if args.molecular:
                from hashlib import sha256
                molecular = json.loads((manifests[0].parent / "molecular_manifest.json").read_text())
                for name, digest in molecular["artifact_sha256"].items():
                    if Path(name).name != name:
                        raise RuntimeError("Unexpected molecular artifact path")
                    with (manifests[0].parent / name).open("rb") as f:
                        if sha256(f.read()).hexdigest() != digest:
                            raise RuntimeError("Molecular output checksum mismatch")
                record.update(molecular_states=molecular["state_rows"],
                              molecular_artifacts=len(molecular["artifact_sha256"]),
                              basis_scope=molecular["basis"]["basis_scope"],
                              mapping_residual=molecular["basis"]["max_abs_residual"])
            summary["datasets"].append(record)
            print(json.dumps(record), flush=True)
        summary["status"] = "success"
    except BaseException:
        summary["status"] = "failed"
        raise
    finally:
        (dest / "execution_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print("Execution record:", dest, flush=True)


if __name__ == "__main__":
    main()
