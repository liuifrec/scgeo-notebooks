"""Memory-safe helpers for Dataset B GSE249479 Phase 1 audits.

The helpers inspect H5AD structure without materializing the expression matrix.
They are intentionally limited to metadata and storage audits for Phase 1.
"""

from __future__ import annotations

import gc
import hashlib
import importlib.metadata as importlib_metadata
import json
import os
import platform
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import psutil


class MemoryLimitExceeded(RuntimeError):
    """Raised when RSS exceeds the configured safety threshold."""


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(root: Path | None = None) -> dict[str, Any]:
    root = repo_root() if root is None else Path(root)
    with (root / "configs" / "gse249479_dataset_b_v1.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_path(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def configured_paths(config: dict[str, Any], root: Path | None = None) -> dict[str, Path]:
    root = repo_root() if root is None else Path(root)
    data_dir = resolve_path(root, os.environ.get(config["data_dir_env"], config["default_data_dir"]))
    input_h5ad = resolve_path(root, os.environ.get(config["input_h5ad_env"], config["default_input_h5ad"]))
    output_dir = resolve_path(root, os.environ.get(config["output_dir_env"], config["default_output_dir"]))
    source_repo = resolve_path(root, os.environ.get(config["source_repository_env"], config["default_source_repository"]))
    return {
        "root": root,
        "data_dir": data_dir,
        "input_h5ad": input_h5ad,
        "output_dir": output_dir,
        "source_repo": source_repo,
    }


def memory_threshold_bytes(config: dict[str, Any]) -> int:
    gb = float(os.environ.get(config["memory_threshold_gb_env"], config["default_memory_threshold_gb"]))
    return int(gb * 1024**3)


def ensure_output_tree(output_dir: Path) -> None:
    for subdir in ["audit", "metadata", "version_records", "execution", "executed_notebooks"]:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)


def rss_bytes() -> int:
    return int(psutil.Process(os.getpid()).memory_info().rss)


def bytes_to_gb(value: int | float | None) -> float | None:
    if value is None:
        return None
    return float(value) / 1024**3


@dataclass
class SectionRecord:
    section: str
    rss_before_bytes: int
    rss_after_bytes: int
    rss_delta_bytes: int
    threshold_bytes: int
    status: str
    timestamp_utc: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "section": self.section,
            "rss_before_gb": bytes_to_gb(self.rss_before_bytes),
            "rss_after_gb": bytes_to_gb(self.rss_after_bytes),
            "rss_delta_gb": bytes_to_gb(self.rss_delta_bytes),
            "threshold_gb": bytes_to_gb(self.threshold_bytes),
            "status": self.status,
            "timestamp_utc": self.timestamp_utc,
        }


class MemoryAudit:
    def __init__(self, threshold_bytes: int):
        self.threshold_bytes = int(threshold_bytes)
        self.records: list[SectionRecord] = []
        self.peak_rss_bytes = rss_bytes()

    @contextmanager
    def section(self, name: str):
        self.close_figures_and_collect()
        before = rss_bytes()
        self.peak_rss_bytes = max(self.peak_rss_bytes, before)
        if before > self.threshold_bytes:
            self._record(name, before, before, "threshold_exceeded_before")
            raise MemoryLimitExceeded(
                f"RSS before section {name!r} is {bytes_to_gb(before):.2f} GB, "
                f"above threshold {bytes_to_gb(self.threshold_bytes):.2f} GB"
            )
        status = "passed"
        try:
            yield
        except Exception:
            status = "failed"
            raise
        finally:
            self.close_figures_and_collect()
            after = rss_bytes()
            self.peak_rss_bytes = max(self.peak_rss_bytes, after)
            self._record(name, before, after, status)
            if after > self.threshold_bytes:
                raise MemoryLimitExceeded(
                    f"RSS after section {name!r} is {bytes_to_gb(after):.2f} GB, "
                    f"above threshold {bytes_to_gb(self.threshold_bytes):.2f} GB"
                )

    def _record(self, name: str, before: int, after: int, status: str) -> None:
        self.records.append(
            SectionRecord(
                section=name,
                rss_before_bytes=before,
                rss_after_bytes=after,
                rss_delta_bytes=after - before,
                threshold_bytes=self.threshold_bytes,
                status=status,
                timestamp_utc=datetime.now(timezone.utc).isoformat(),
            )
        )

    @staticmethod
    def close_figures_and_collect() -> None:
        try:
            import matplotlib.pyplot as plt

            plt.close("all")
        except Exception:
            pass
        gc.collect()

    def dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([record.as_dict() for record in self.records])


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


def active_git_branch(path: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(path), "branch", "--show-current"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def require_active_branch(root: Path, expected_branch: str) -> str:
    branch = active_git_branch(root)
    if branch != expected_branch:
        raise RuntimeError(f"Active git branch must be {expected_branch!r}; observed {branch!r}.")
    return branch


def package_version(name: str) -> str | None:
    if name == "python":
        return platform.python_version()
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return None


def write_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n", encoding="utf-8")
    os.replace(tmp_path, path)
    return path


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if pd.isna(value) if np.isscalar(value) else False:
        return None
    raise TypeError(f"{type(value).__name__} is not JSON serializable")


def relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def source_file_provenance(path: Path, root: Path) -> dict[str, Any]:
    path = path.resolve()
    exists = path.exists()
    payload: dict[str, Any] = {
        "absolute_path": str(path),
        "display_path": relative_or_absolute(path, root),
        "exists": bool(exists),
        "immutable_source": True,
        "write_mode_used": False,
        "overwritten": False,
    }
    if exists:
        stat = path.stat()
        payload.update(
            {
                "file_size_bytes": int(stat.st_size),
                "file_size_gb": bytes_to_gb(stat.st_size),
                "modified_time_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                "sha256": sha256_file(path),
            }
        )
    else:
        payload.update(
            {
                "file_size_bytes": None,
                "file_size_gb": None,
                "modified_time_utc": None,
                "sha256": None,
            }
        )
    return payload


def _decode_attr(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        return [_decode_attr(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def h5_shape(value: h5py.Dataset | h5py.Group) -> tuple[int, ...] | None:
    if isinstance(value, h5py.Dataset):
        return tuple(int(x) for x in value.shape)
    shape = value.attrs.get("shape")
    if shape is None:
        return None
    return tuple(int(x) for x in shape)


def dense_size_bytes(shape: tuple[int, ...] | None, dtype: str) -> int | None:
    if not shape:
        return None
    return int(np.prod(shape, dtype=np.int64) * np.dtype(dtype).itemsize)


def h5_matrix_info(handle: h5py.File, key: str) -> dict[str, Any]:
    if key not in handle:
        return {"key": key, "present": False}
    value = handle[key]
    if isinstance(value, h5py.Dataset):
        shape = h5_shape(value)
        dtype = str(value.dtype)
        itemsize = int(value.dtype.itemsize)
        dense_bytes = int(np.prod(shape, dtype=np.int64) * itemsize) if shape else None
        return {
            "key": key,
            "present": True,
            "storage": "dense_dataset",
            "sparse_format": None,
            "encoding_type": _decode_attr(value.attrs.get("encoding-type", "dense")),
            "shape": shape,
            "dtype": dtype,
            "nnz": None,
            "data_bytes": int(value.size * itemsize),
            "indices_bytes": None,
            "indptr_bytes": None,
            "estimated_sparse_bytes": None,
            "estimated_dense_bytes": dense_bytes,
            "estimated_dense_float32_bytes": dense_size_bytes(shape, "float32"),
            "estimated_dense_float64_bytes": dense_size_bytes(shape, "float64"),
            "property_status": "available_without_materializing_X",
        }
    encoding = _decode_attr(value.attrs.get("encoding-type", "group"))
    sparse_format = None
    if encoding in {"csr_matrix", "csc_matrix"}:
        sparse_format = encoding.replace("_matrix", "")
    shape = h5_shape(value)
    data = value.get("data")
    indices = value.get("indices")
    indptr = value.get("indptr")
    data_bytes = int(data.size * data.dtype.itemsize) if data is not None else None
    indices_bytes = int(indices.size * indices.dtype.itemsize) if indices is not None else None
    indptr_bytes = int(indptr.size * indptr.dtype.itemsize) if indptr is not None else None
    sparse_bytes = sum(x for x in [data_bytes, indices_bytes, indptr_bytes] if x is not None)
    dense_bytes = None
    if shape and data is not None:
        dense_bytes = int(np.prod(shape, dtype=np.int64) * data.dtype.itemsize)
    return {
        "key": key,
        "present": True,
        "storage": "sparse_group" if data is not None and indices is not None and indptr is not None else "group",
        "encoding_type": encoding,
        "sparse_format": sparse_format,
        "shape": shape,
        "dtype": str(data.dtype) if data is not None else None,
        "index_dtype": str(indices.dtype) if indices is not None else None,
        "indptr_dtype": str(indptr.dtype) if indptr is not None else None,
        "nnz": int(data.size) if data is not None else None,
        "data_bytes": data_bytes,
        "indices_bytes": indices_bytes,
        "indptr_bytes": indptr_bytes,
        "estimated_sparse_bytes": sparse_bytes if sparse_bytes else None,
        "estimated_dense_bytes": dense_bytes,
        "estimated_dense_float32_bytes": dense_size_bytes(shape, "float32"),
        "estimated_dense_float64_bytes": dense_size_bytes(shape, "float64"),
        "property_status": "available_without_materializing_X" if data is not None else "unavailable_without_loading_matrix",
    }


def inspect_h5ad_storage(path: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        x_info = h5_matrix_info(handle, "X")
        layers = []
        if "layers" in handle:
            for name in sorted(handle["layers"].keys()):
                info = h5_matrix_info(handle, f"layers/{name}")
                info["layer"] = name
                layers.append(info)
        raw_info = {"present": "raw" in handle}
        if "raw" in handle:
            raw_info["X"] = h5_matrix_info(handle, "raw/X")
            raw_info["keys"] = sorted(handle["raw"].keys())
        return {
            "file": str(path),
            "file_size_bytes": path.stat().st_size,
            "X": x_info,
            "layers": layers,
            "raw": raw_info,
            "top_level_keys": sorted(handle.keys()),
        }


def inspect_obs_metadata(path: Path, config: dict[str, Any]) -> dict[str, Any]:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        obs = adata.obs.copy()
        var = adata.var.copy()
        var_columns = list(adata.var.columns)
        obs_columns = list(obs.columns)
        obs_names = pd.Index(adata.obs_names)
        var_names = pd.Index(adata.var_names)
        lower_map = {col: col.lower() for col in obs_columns}
        qc_columns = [
            col
            for col, low in lower_map.items()
            if any(keyword.lower() in low for keyword in config["qc_field_keywords"])
        ]
        replicate_columns = [
            col
            for col, low in lower_map.items()
            if any(keyword.lower() in low for keyword in config["replicate_metadata_keywords"])
        ]
        clade_columns = [
            col
            for col, low in lower_map.items()
            if any(keyword.lower() in low for keyword in config["clade_metadata_keywords"])
        ]
        condition_columns = [col for col in config["condition_column_candidates"] if col in obs.columns]
        condition_counts: list[dict[str, Any]] = []
        for col in condition_columns:
            counts = obs[col].astype(str).value_counts(dropna=False)
            for value, count in counts.items():
                condition_counts.append({"column": col, "value": value, "n_cells": int(count)})
        metadata_summary = []
        for col in sorted(set(replicate_columns + clade_columns + condition_columns + qc_columns)):
            series = obs[col]
            metadata_summary.append(
                {
                    "column": col,
                    "dtype": str(series.dtype),
                    "n_missing": int(series.isna().sum()),
                    "n_unique": int(series.astype(str).nunique(dropna=False)),
                    "example_values": "; ".join(series.astype(str).drop_duplicates().head(8).tolist()),
                }
            )
        return {
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "obs_names_unique": bool(obs_names.is_unique),
            "var_names_unique": bool(var_names.is_unique),
            "n_duplicate_obs_names": int(obs_names.duplicated().sum()),
            "n_duplicate_var_names": int(var_names.duplicated().sum()),
            "obs_columns": obs_columns,
            "var_columns": var_columns,
            "qc_columns": qc_columns,
            "replicate_candidate_columns": replicate_columns,
            "clade_candidate_columns": clade_columns,
            "condition_candidate_columns": condition_columns,
            "condition_counts": condition_counts,
            "metadata_summary": metadata_summary,
            "obs": obs,
            "var": var,
        }
    finally:
        adata.file.close()
        del adata
        gc.collect()


def _field_category(field_name: str, config: dict[str, Any], condition_columns: list[str]) -> str:
    low = field_name.lower()
    if field_name in condition_columns:
        return "condition"
    for category, keywords in config["replication_field_categories"].items():
        if any(keyword.lower() in low for keyword in keywords):
            return category
    return "other_metadata"


def _candidate_field_names(obs: pd.DataFrame, config: dict[str, Any], condition_columns: list[str]) -> list[str]:
    candidates: list[str] = []
    for col in obs.columns:
        category = _field_category(str(col), config, condition_columns)
        if category != "other_metadata":
            candidates.append(str(col))
    return candidates


def _values_per_condition(obs: pd.DataFrame, field: str, condition_columns: list[str]) -> tuple[str, str]:
    if not condition_columns:
        return "unavailable_no_condition_field", "unavailable"
    summaries: dict[str, dict[str, list[str]]] = {}
    nested_flags: list[bool] = []
    field_values = obs[field].astype(str)
    for condition in condition_columns:
        cond_values = obs[condition].astype(str)
        condition_summary: dict[str, list[str]] = {}
        for cond_value, idx in cond_values.groupby(cond_values, observed=False).groups.items():
            vals = sorted(field_values.loc[idx].dropna().unique().tolist())
            condition_summary[str(cond_value)] = vals[:25]
        summaries[condition] = condition_summary

        cross = pd.crosstab(field_values, cond_values)
        if cross.shape[0] == 0 or cross.shape[1] == 0:
            nested_flags.append(False)
        else:
            each_field_in_one_condition = bool((cross.gt(0).sum(axis=1) <= 1).all())
            each_condition_in_one_field = bool((cross.gt(0).sum(axis=0) <= 1).all())
            nested_flags.append(each_field_in_one_condition or each_condition_in_one_field)
    nested_status = "yes" if any(nested_flags) else "no"
    return json.dumps(summaries, sort_keys=True), nested_status


def _eligibility_decision(category: str, nested_status: str) -> tuple[str, str]:
    if category in {"biological_donor", "xenograft_recipient"}:
        if nested_status == "yes":
            return (
                "not_eligible_confounded_with_condition",
                "Biological replicate-like field is nested or confounded with condition; do not use as replicate unit without external confirmation.",
            )
        return (
            "candidate_requires_external_confirmation",
            "Field name is compatible with a biological donor or xenograft recipient unit, but Phase 1 does not infer validity beyond metadata.",
        )
    if category == "source_pool":
        return ("not_eligible_source_pool", "Source pool is not automatically a biological replicate.")
    if category == "genetic_clade":
        return ("not_eligible_genetic_clade", "Genetic clade is lineage structure, not a biological replicate unit.")
    if category == "experimental_library":
        return ("not_eligible_library", "Experimental library is a technical unit and is not automatically a biological replicate.")
    if category == "condition":
        return ("not_eligible_condition", "Condition labels are not biological replicate units.")
    if category == "inferred_cell_label":
        return ("not_eligible_inferred_cell_label", "Inferred cell labels are not biological replicate units.")
    return ("not_eligible_unknown", "Field is not a prespecified biological replicate category.")


def candidate_replicate_audit(metadata: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    obs = metadata.get("obs")
    if obs is None:
        return pd.DataFrame(
            columns=[
                "field_name",
                "field_category",
                "n_unique_values",
                "values_per_condition",
                "missing_fraction",
                "nested_or_confounded_with_condition",
                "final_eligibility_decision",
                "eligibility_reason",
            ]
        )
    condition_columns = metadata.get("condition_candidate_columns", [])
    rows: list[dict[str, Any]] = []
    for field in _candidate_field_names(obs, config, condition_columns):
        series = obs[field]
        category = _field_category(field, config, condition_columns)
        values_per_condition, nested_status = _values_per_condition(obs, field, condition_columns)
        decision, reason = _eligibility_decision(category, nested_status)
        rows.append(
            {
                "field_name": field,
                "field_category": category,
                "n_unique_values": int(series.astype(str).nunique(dropna=False)),
                "values_per_condition": values_per_condition,
                "missing_fraction": float(series.isna().mean()),
                "nested_or_confounded_with_condition": nested_status,
                "final_eligibility_decision": decision,
                "eligibility_reason": reason,
            }
        )
    return pd.DataFrame(rows)


def replication_decision(candidate_fields: pd.DataFrame) -> dict[str, Any]:
    if candidate_fields.empty:
        return {
            "valid_biological_replicate_units": False,
            "decision": "no_explicit_replicate_metadata_found",
            "reason": "No candidate metadata fields were found. Cells are not treated as replicates.",
        }
    eligible = candidate_fields[
        candidate_fields["final_eligibility_decision"].eq("candidate_requires_external_confirmation")
    ]
    if not eligible.empty:
        return {
            "valid_biological_replicate_units": "candidate_requires_confirmation",
            "decision": "explicit_biological_replicate_candidates_present",
            "candidate_columns": eligible["field_name"].tolist(),
            "reason": "Biological donor or xenograft recipient candidate metadata exists, but Phase 1 requires external confirmation before use as replicate units.",
        }
    return {
        "valid_biological_replicate_units": False,
        "decision": "no_eligible_biological_replicate_field",
        "candidate_columns": candidate_fields["field_name"].tolist(),
        "reason": "Candidate fields are condition, clade, library, source-pool, inferred-label, or confounded biological-replicate-like fields. Do not treat cells, libraries, or condition labels as biological replicates.",
    }


def feasibility_assessment(matrix_info: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    x_info = matrix_info.get("X", {})
    shape = x_info.get("shape")
    dense_bytes = x_info.get("estimated_dense_bytes")
    dense_gb = bytes_to_gb(dense_bytes) if dense_bytes is not None else None
    dense_float32_gb = bytes_to_gb(x_info.get("estimated_dense_float32_bytes"))
    dense_float64_gb = bytes_to_gb(x_info.get("estimated_dense_float64_bytes"))
    sparse_gb = bytes_to_gb(x_info.get("estimated_sparse_bytes"))
    stop_gb = float(config["stop_before_estimated_rss_gb"])
    rows = []
    rows.append(
        {
            "step": "PCA",
            "phase1_status": "not_run",
            "estimated_feasibility": "requires_incremental_or_sparse_aware_plan" if dense_gb and dense_gb >= stop_gb else "likely_feasible_with_sparse_workflow",
            "memory_basis": (
                f"dense current dtype {dense_gb:.2f} GB; dense float32 {dense_float32_gb:.2f} GB; "
                f"dense float64 {dense_float64_gb:.2f} GB; sparse storage {sparse_gb:.2f} GB"
                if dense_gb is not None and sparse_gb is not None and dense_float32_gb is not None and dense_float64_gb is not None
                else "input unavailable"
            ),
            "required_next_action": "Use backed/sparse-aware preprocessing, avoid full dense scaling, and checkpoint before PCA.",
        }
    )
    rows.append(
        {
            "step": "Scanorama",
            "phase1_status": "not_run",
            "estimated_feasibility": "high_risk_until_batch_sizes_known",
            "memory_basis": "Integration can allocate dense corrected embeddings and pairwise workspace; no execution in Phase 1.",
            "required_next_action": "Run only after sample/batch metadata is validated and a memory estimate is below 10 GB RSS.",
        }
    )
    rows.append(
        {
            "step": "Harmony",
            "phase1_status": "not_run",
            "estimated_feasibility": "possibly_feasible_after_pca_only",
            "memory_basis": "Harmony operates on PCA embeddings, not full counts, but valid batch/replicate metadata is required.",
            "required_next_action": "Defer until PCA checkpoint and biological replicate/batch columns are confirmed.",
        }
    )
    rows.append(
        {
            "step": "scVI",
            "phase1_status": "not_run",
            "estimated_feasibility": "deferred_first_pass",
            "memory_basis": "User requested no scVI in first pass; GPU/CPU memory and DataLoader worker plan must be assessed later.",
            "required_next_action": "Prepare a separate scVI pass with small DataLoader worker counts and RSS monitoring.",
        }
    )
    if shape:
        for row in rows:
            row["n_cells"] = int(shape[0])
            row["n_genes"] = int(shape[1])
    return pd.DataFrame(rows)


def write_metadata(output_dir: Path, config: dict[str, Any], paths: dict[str, Path], extra: dict[str, Any]) -> Path:
    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "workflow": config["workflow_name"],
        "schema_version": config["schema_version"],
        "phase": config["phase"],
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "notebook_repository_commit": git_commit(paths["root"]),
        "source_repository_commit": git_commit(paths["source_repo"]) if paths["source_repo"].exists() else None,
        "scgeo_package_modified": False,
        "frozen_scgeo_thresholds_tuned": False,
        "rna_velocity_forced": False,
        "full_expression_matrix_densified": False,
        "packages": {name: package_version(name) for name in config["version_record_packages"]},
    }
    payload.update(extra)
    return write_json(output_dir / "metadata" / "00_metadata_replication_and_memory_audit_metadata.json", payload)
