"""Shared execution/export helpers for the independent 0.2.0 tutorial track."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

VERSION = "0.2.0"
RELEASE_COMMIT = "46ba9017ceef17255292edd0140602ed13ac5305"
FORMATS = ("png", "svg", "pdf")


def verify_examples(root: Path):
    """Check that the copied example sources are exactly the published tag's blobs."""
    example_dir = root / "scripts" / "scgeo_v020_examples"
    records = []
    for part in ("SOURCE_PART1.json", "SOURCE_PART2.json"):
        source = json.loads((example_dir / part).read_text(encoding="utf-8"))
        if source["repository"] != "liuifrec/scgeo" or source["commit"] != RELEASE_COMMIT:
            raise RuntimeError("Unexpected example provenance")
        for record in source["files"]:
            path = root / record["path"]
            if path.parent != example_dir or path.suffix != ".py":
                raise RuntimeError("Unexpected example path")
            data = path.read_bytes()
            blob_sha = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
            if blob_sha != record["git_blob_sha"]:
                raise RuntimeError(f"Example source differs from release: {path.name}")
            records.append(record)
    if len(records) != 13 or len({r["path"] for r in records}) != 13:
        raise RuntimeError("Incomplete release example inventory")
    return records


def setup(root: Path, lesson: str) -> Path:
    import scgeo as sg
    if sg.__version__ != VERSION or importlib.metadata.version("scgeo") != VERSION:
        raise RuntimeError("Use the separate ScGeo 0.2.0 tutorial environment and restart the kernel.")
    expected_python = os.environ.get("SCGEO_TUTORIAL_EXPECTED_PYTHON")
    if expected_python and Path(sys.executable).resolve() != Path(expected_python).resolve():
        raise RuntimeError("Notebook is running in the wrong Python interpreter.")
    records = verify_examples(root)
    explicit = os.environ.get("SCGEO_TUTORIAL_OUTPUT_ROOT")
    base = Path(explicit).expanduser().resolve() if explicit else root / "results" / "scgeo_v020_tutorials"
    base.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix=f"{lesson}-", dir=base))
    try:
        commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"],
                                         text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    versions = {name: importlib.metadata.version(name) for name in
                ("scgeo", "anndata", "numpy", "pandas", "scipy", "scikit-learn",
                 "matplotlib", "plotly", "nbformat", "nbclient", "ipykernel")}
    metadata = {
        "lesson": lesson, "package_release": VERSION, "release_source_commit": RELEASE_COMMIT,
        "source_examples": records, "companion_commit": commit,
        "python": sys.version, "python_executable": sys.executable,
        "scgeo_import": sg.__file__, "platform": platform.platform(),
        "packages": versions,
        "scope": "tutorial; synthetic unless explicitly supplied public inputs; not manuscript reproduction",
    }
    (out / "tutorial_provenance.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"ScGeo {sg.__version__} from {sg.__file__}\nOutputs: {out}")
    return out


def save_preview(result, out, name, *, scope="fully synthetic tutorial; no biological claims"):
    """Export the actual plotted sources and show the rendered PNG in the notebook."""
    import matplotlib.pyplot as plt
    from IPython.display import Image, display
    result.metadata["tutorial_scope"] = scope
    result.metadata["package_release_source_commit"] = RELEASE_COMMIT
    files = result.save(out / name, formats=FORMATS, dpi=145)
    display(Image(filename=str(out / f"{name}.png"), width=850))
    plt.close(result.fig)
    return files


def show_gallery(directory, *, limit=2):
    """Show actual exported images, not a hand-recreated preview."""
    from IPython.display import Image, display
    directory = Path(directory)
    previews = sorted(directory.rglob("*.png"))
    if not previews:
        raise RuntimeError(f"No gallery PNG exports found in {directory}")
    print(f"{len(previews)} PNG figures; complete source tables and metadata are under {directory}")
    for path in previews[:limit]:
        print(path.name)
        display(Image(filename=str(path), width=850))
