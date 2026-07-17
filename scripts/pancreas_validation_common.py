"""Shared utilities for the public pancreas Dataset D validation notebooks.

These helpers handle paths, checksums, artifact writing, and the prespecified
cosine/consensus bookkeeping used by the public validation notebooks. They do
not modify ScGeo and do not change any synthetic benchmark protocol settings.
"""

from __future__ import annotations

import hashlib
import importlib.metadata as importlib_metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(root: Path | None = None) -> dict[str, Any]:
    root = repo_root() if root is None else Path(root)
    with (root / "configs" / "pancreas_dataset_d_v1.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_path(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def configured_paths(config: dict[str, Any], root: Path | None = None) -> dict[str, Path]:
    root = repo_root() if root is None else Path(root)
    data_dir = resolve_path(
        root,
        os.environ.get(config["data_dir_env"], config["default_data_dir"]),
    )
    output_dir = resolve_path(
        root,
        os.environ.get(config["output_dir_env"], config["default_output_dir"]),
    )
    source_repo = resolve_path(
        root,
        os.environ.get(config["source_repository_env"], config["default_source_repository"]),
    )
    return {"root": root, "data_dir": data_dir, "output_dir": output_dir, "source_repo": source_repo}


def ensure_runtime_env(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR", str(output_dir / "_numba_cache"))
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / "_matplotlib_cache"))
    os.environ.setdefault("PYTHONHASHSEED", "0")
    Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)


def ensure_output_tree(output_dir: Path) -> None:
    for subdir in [
        "figures",
        "figure_sources",
        "alt_text",
        "metadata",
        "version_records",
        "intermediates",
        "execution",
    ]:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)


def rel_display(path: Path, root: Path | None = None) -> str:
    root = repo_root() if root is None else Path(root)
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(path: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def package_version(name: str) -> str | None:
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return None


def package_versions(names: list[str]) -> dict[str, str | None]:
    return {name: package_version(name) for name in names}


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if pd.isna(value) if np.isscalar(value) else False:
        return None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def write_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n", encoding="utf-8")
    return path


def write_dataframe(output_dir: Path, stem: str, df: pd.DataFrame) -> Path:
    path = output_dir / "figure_sources" / f"{stem}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def write_alt_text(output_dir: Path, stem: str, text: str) -> Path:
    path = output_dir / "alt_text" / f"{stem}.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized = " ".join(str(text).split())
    path.write_text(normalized + "\n", encoding="utf-8")
    return path


def save_figure(output_dir: Path, stem: str, fig: Any) -> dict[str, str]:
    png_path = output_dir / "figures" / f"{stem}.png"
    svg_path = output_dir / "figures" / f"{stem}.svg"
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, format="svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(png_path, format="png", bbox_inches="tight", dpi=200, metadata={"Software": "matplotlib"})
    return {"png": str(png_path), "svg": str(svg_path)}


def write_metadata(
    output_dir: Path,
    stem: str,
    config: dict[str, Any],
    extra: dict[str, Any] | None = None,
) -> Path:
    paths = configured_paths(config)
    payload: dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "workflow": config["workflow_name"],
        "schema_version": config["schema_version"],
        "dataset": config["dataset"]["name"],
        "public_dataset": True,
        "biological_validation_scope": "public pancreas developmental-dynamics validation",
        "no_artificial_treatment_control": True,
        "cellrank_velocity_kernel_independent_of_scvelo": False,
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "notebook_repository_commit": git_commit(paths["root"]),
        "source_repository_commit": git_commit(paths["source_repo"]) if paths["source_repo"].exists() else None,
        "scgeo_package_modified": False,
        "frozen_synthetic_protocol_modified": False,
        "threshold_policy": config["threshold_policy"],
        "packages": package_versions(config["version_record_packages"]),
    }
    if extra:
        payload.update(extra)
    return write_json(output_dir / "metadata" / f"{stem}_metadata.json", payload)


def version_record(output_dir: Path, stem: str, config: dict[str, Any], extra: dict[str, Any] | None = None) -> Path:
    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": config["schema_version"],
        "packages": package_versions(config["version_record_packages"]),
    }
    if extra:
        payload.update(extra)
    return write_json(output_dir / "version_records" / f"{stem}_versions.json", payload)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if not np.isfinite(denom) or denom <= 1e-12:
        return float("nan")
    return float(np.dot(a, b) / denom)


def classify_cosine(value: float, config: dict[str, Any]) -> str:
    thresholds = config["scgeo_alignment_defaults"]
    if not np.isfinite(value):
        return "unavailable"
    if value >= float(thresholds["alignment_pos_thr"]):
        return "aligned"
    if value <= float(thresholds["alignment_neg_thr"]):
        return "discordant"
    return "neutral"


def consensus_from_rows(rows: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    principal = rows[rows["principal_representation"].astype(bool)].copy()
    usable = principal[principal["status"].eq("usable") & principal["class"].isin(["aligned", "discordant", "neutral"])]
    n_usable = int(usable.shape[0])
    min_usable = int(config["scgeo_alignment_defaults"]["min_usable_representations"])
    if n_usable < min_usable:
        return {
            "consensus_class": "unavailable",
            "n_usable_representations": n_usable,
            "agreement_fraction": np.nan,
            "status_reason": f"fewer than {min_usable} usable principal representations",
        }
    counts = usable["class"].value_counts()
    winner = str(counts.index[0])
    fraction = float(counts.iloc[0] / n_usable)
    threshold = float(config["scgeo_alignment_defaults"]["class_fraction_threshold"])
    if fraction >= threshold:
        return {
            "consensus_class": winner,
            "n_usable_representations": n_usable,
            "agreement_fraction": fraction,
            "status_reason": f"{winner} fraction {fraction:.2f} across principal representations",
        }
    return {
        "consensus_class": "unstable",
        "n_usable_representations": n_usable,
        "agreement_fraction": fraction,
        "status_reason": f"no class reaches prespecified fraction {threshold:.2f}",
    }


def clean_label(value: Any) -> str:
    return str(value).replace("/", "_").replace(" ", "_").replace("+", "plus").replace("-", "_")
