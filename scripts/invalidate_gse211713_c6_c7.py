#!/usr/bin/env python3
"""Archive the invalidated GSE211713 C6/C7 run without modifying C0-C5."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "results/public_validation/gse211713_dataset_c"
DEST = RESULT_ROOT / "invalidated_runs/c6_c7_model_switch_run"
OLD_REPRESENTATION = Path("/home/liuyuchen/data/gse211713/gse211713_revision_representations.h5ad")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_files() -> list[Path]:
    files: set[Path] = set()
    for directory in ["representations", "scgeo", "controls", "figures", "figure_sources", "alt_text"]:
        source = RESULT_ROOT / directory
        if source.exists():
            files.update(path for path in source.rglob("*") if path.is_file())
    model = RESULT_ROOT / "models/scvi"
    if model.exists():
        files.update(path for path in model.rglob("*") if path.is_file())
    execution = RESULT_ROOT / "execution"
    for pattern in ["04*.json", "05*.json", "gse211713_phase_c6*.json", "gse211713_phase_c7*.json"]:
        files.update(path for path in execution.glob(pattern) if path.is_file())
    metadata = RESULT_ROOT / "metadata"
    for pattern in [
        "gse211713_pca*", "gse211713_diffusion*", "gse211713_scvi*",
        "gse211713_representation*", "gse211713_scgeo*",
        "gse211713_secondary*", "control_vs_17gy_*", "17gy_early_vs_late*",
    ]:
        files.update(path for path in metadata.glob(pattern) if path.is_file())
    executed = RESULT_ROOT / "executed_notebooks/notebooks/public_validation/gse211713"
    for pattern in ["04*.ipynb", "05*.ipynb"]:
        files.update(path for path in executed.glob(pattern) if path.is_file())
    return sorted(files)


def main() -> int:
    if DEST.exists():
        raise RuntimeError(f"Invalidation destination already exists: {DEST}")
    DEST.mkdir(parents=True)
    records: list[dict[str, object]] = []
    for source in selected_files():
        relative = source.relative_to(RESULT_ROOT)
        target = DEST / "workspace_artifacts" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        records.append(
            {
                "original_path": str(source),
                "archived_path": str(target.relative_to(DEST)),
                "size_bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        )
    if OLD_REPRESENTATION.is_file():
        target = DEST / "external_artifacts" / OLD_REPRESENTATION.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OLD_REPRESENTATION, target)
        records.append(
            {
                "original_path": str(OLD_REPRESENTATION),
                "archived_path": str(target.relative_to(DEST)),
                "size_bytes": OLD_REPRESENTATION.stat().st_size,
                "sha256": sha256(OLD_REPRESENTATION),
            }
        )
    manifest = {
        "schema_version": "gse211713_c6_c7_invalidation_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_status": "invalidated_noncanonical",
        "reason": "The model-switch C6/C7 run failed strict acceptance and is not canonical.",
        "known_defects": [
            "per-mouse coverage eligibility was not enforced before ScGeo",
            "biological-mouse resampling was mislabeled as cell-level resampling",
            "PCA20 and diffusion ScGeo sensitivity results were not persisted",
            "summary construction recursively duplicated evidence rows 31-fold",
            "smoke-test and unintended states entered aggregate outputs",
            "scVI maximum epochs changed from 100 to 40 after an export failure",
            "scVI history was malformed",
            "resource and GPU peak records were incomplete",
            "aggregate C6/C7 workflow reports retained failed status",
        ],
        "frozen_c0_c5_artifacts_moved": False,
        "files": records,
    }
    path = DEST / "invalidation_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "archived", "files": len(records), "manifest": str(path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
