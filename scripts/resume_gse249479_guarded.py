#!/usr/bin/env python3
"""Cautiously resume GSE249479 public MEX reconstruction.

This runner intentionally avoids the original all-in-one 00a notebook path.
It validates/resumes downloads, reconstructs one sample per fresh subprocess,
and concatenates guarded sparse H5AD outputs only after compatibility checks.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import resource
import shutil
import signal
import subprocess
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


USER_AGENT = "scgeo-gse249479-guarded-resume/1.0"
CHUNK_BYTES = 1024 * 1024
RNA_FILE_ROLES = ("barcodes", "features", "matrix")
REPORT_NAME = "gse249479_guarded_reconstruction_report.json"


@dataclass(frozen=True)
class DownloadTarget:
    sample_accession: str
    condition: str
    role: str
    filename: str
    official_url: str
    expected_size_bytes: int
    file_type: str | None
    local_path: Path


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return str(value)


def atomic_write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with tmp_path.open("w", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_path, path)
    return path


def atomic_write_json(path: Path, payload: dict[str, Any]) -> Path:
    return atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n")


def atomic_write_csv(path: Path, rows: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    tmp_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_path, path)
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text_lines(lines: list[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def rss_bytes() -> int:
    try:
        import psutil

        return int(psutil.Process(os.getpid()).memory_info().rss)
    except Exception:
        with Path("/proc/self/statm").open("r", encoding="utf-8") as handle:
            pages = int(handle.read().split()[1])
        return pages * os.sysconf("SC_PAGE_SIZE")


def peak_rss_bytes() -> int:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    # Linux reports ru_maxrss in KiB.
    return int(usage.ru_maxrss) * 1024


def bytes_to_gb(value: int | float | None) -> float | None:
    if value is None:
        return None
    return float(value) / 1024**3


def run_capture(args: list[str]) -> dict[str, Any]:
    started = time.perf_counter()
    try:
        completed = subprocess.run(args, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return {
            "args": args,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "runtime_seconds": round(time.perf_counter() - started, 3),
        }
    except Exception as exc:
        return {
            "args": args,
            "returncode": None,
            "stdout": "",
            "stderr": repr(exc),
            "runtime_seconds": round(time.perf_counter() - started, 3),
        }


def load_config(root: Path) -> dict[str, Any]:
    with (root / "configs" / "gse249479_dataset_b_v1.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_path(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def sample_bucket(accession: str) -> str:
    prefix = "".join(ch for ch in accession if ch.isalpha())
    digits = "".join(ch for ch in accession if ch.isdigit())
    return f"{prefix}{digits[:-3]}nnn"


def url_open(url: str, *, headers: dict[str, str] | None = None, timeout: int = 120):
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)
    request = urllib.request.Request(url, headers=request_headers)
    return urllib.request.urlopen(request, timeout=timeout)


def fetch_text(url: str, *, timeout: int = 60) -> str:
    with url_open(url, timeout=timeout) as response:
        body = response.read()
    return body.decode("utf-8", errors="replace")


def parse_geo_filelist(text: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        name = (row.get("Name") or "").strip()
        if not name:
            continue
        size_text = (row.get("Size") or "").strip()
        rows.append(
            {
                "record_kind": (row.get("#Archive/File") or "").strip(),
                "filename": name,
                "modified_time_geo": (row.get("Time") or "").strip(),
                "size_bytes": int(size_text) if size_text.isdigit() else None,
                "file_type": (row.get("Type") or "").strip() or None,
            }
        )
    return rows


def expected_targets(config: dict[str, Any], filelist_rows: list[dict[str, Any]], data_dir: Path) -> list[DownloadTarget]:
    by_name = {row["filename"]: row for row in filelist_rows if str(row["record_kind"]).lower() == "file"}
    targets: list[DownloadTarget] = []
    template = config["geo"]["sample_url_template"]
    for sample in config["geo"]["samples"]:
        accession = sample["accession"]
        for role in RNA_FILE_ROLES:
            filename = sample["expected_files"][role]
            if filename not in by_name:
                raise RuntimeError(f"Expected GEO file is absent from filelist.txt: {filename}")
            record = by_name[filename]
            if record["size_bytes"] is None:
                raise RuntimeError(f"GEO filelist has no byte size for {filename}")
            targets.append(
                DownloadTarget(
                    sample_accession=accession,
                    condition=sample["condition"],
                    role=role,
                    filename=filename,
                    official_url=template.format(
                        sample_bucket=sample_bucket(accession),
                        accession=accession,
                        filename=urllib.parse.quote(filename),
                    ),
                    expected_size_bytes=int(record["size_bytes"]),
                    file_type=record["file_type"],
                    local_path=data_dir / "source" / accession / filename,
                )
            )
    return targets


def file_record(path: Path, *, root: Path | None = None, include_sha256: bool = False) -> dict[str, Any]:
    stat_result = path.stat()
    display = str(path)
    if root is not None:
        try:
            display = str(path.resolve().relative_to(root.resolve()))
        except ValueError:
            display = str(path)
    record = {
        "path": str(path),
        "display_path": display,
        "size_bytes": int(stat_result.st_size),
        "modified_time_utc": datetime.fromtimestamp(stat_result.st_mtime, timezone.utc).isoformat(),
    }
    if include_sha256:
        record["sha256"] = sha256_file(path)
    return record


def collect_existing_files(data_dir: Path) -> list[dict[str, Any]]:
    if not data_dir.exists():
        return []
    rows = []
    for path in sorted(p for p in data_dir.rglob("*") if p.is_file()):
        rows.append(file_record(path, root=data_dir, include_sha256=False))
    return rows


def collect_prior_artifacts(output_dir: Path, data_dir: Path) -> list[dict[str, Any]]:
    roots = [output_dir, data_dir]
    suffixes = {".json", ".csv", ".log", ".txt", ".md", ".sha256"}
    rows: list[dict[str, Any]] = []
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(p for p in root.rglob("*") if p.is_file()):
            if any(part.startswith("_") and "cache" in part for part in path.parts):
                continue
            name = path.name.lower()
            if path.suffix.lower() not in suffixes and not any(token in name for token in ("report", "reconstruct", "download", "memory", "audit")):
                continue
            record = file_record(path, root=root, include_sha256=False)
            if path.suffix.lower() == ".json" and path.stat().st_size <= 2_000_000 and any(
                token in name for token in ("report", "reconstruct", "download", "memory", "audit", "provenance", "metadata")
            ):
                try:
                    record["json_preview"] = json.loads(path.read_text(encoding="utf-8"))
                except Exception as exc:
                    record["json_preview_error"] = repr(exc)
            rows.append(record)
    return rows


def disk_usage_record(path: Path) -> dict[str, Any]:
    usage = shutil.disk_usage(path)
    return {
        "path": str(path),
        "total_bytes": int(usage.total),
        "used_bytes": int(usage.used),
        "free_bytes": int(usage.free),
        "free_gb": bytes_to_gb(usage.free),
    }


def collect_environment_inventory(data_dir: Path, output_dir: Path) -> dict[str, Any]:
    data_root = data_dir
    if not data_root.exists():
        data_root = data_dir.parent if data_dir.parent.exists() else Path.home()
    try:
        import psutil

        memory = psutil.virtual_memory()._asdict()
        swap = psutil.swap_memory()._asdict()
    except Exception:
        memory = {}
        swap = {}
    return {
        "timestamp_utc": now_utc(),
        "platform": platform.platform(),
        "python": sys.version,
        "python_executable": sys.executable,
        "current_rss_gb": bytes_to_gb(rss_bytes()),
        "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
        "memory": memory,
        "swap": swap,
        "commands": {
            "free_h": run_capture(["free", "-h"]),
            "swapon_show": run_capture(["swapon", "--show"]),
            "df_h": run_capture(["df", "-h", str(data_root), "/tmp"]),
        },
        "disk_usage": {
            "data_root": disk_usage_record(data_root),
            "tmp": disk_usage_record(Path("/tmp")),
        },
        "existing_files": collect_existing_files(data_dir),
        "prior_artifacts": collect_prior_artifacts(output_dir, data_dir),
    }


def load_known_sha256(output_dir: Path) -> dict[str, str]:
    known: dict[str, str] = {}
    manifest_path = output_dir / "audit" / "gse249479_guarded_download_manifest.json"
    if manifest_path.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            for row in payload.get("downloads", []):
                filename = str(row.get("filename") or "")
                digest = row.get("expected_sha256") or row.get("observed_sha256")
                if filename and digest:
                    known[filename] = str(digest)
        except Exception:
            pass
    return known


def gzip_magic_ok(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def gzip_stream_test(path: Path) -> None:
    if not gzip_magic_ok(path):
        raise RuntimeError(f"gzip magic bytes are invalid: {path}")
    with gzip.open(path, "rb") as handle:
        while handle.read(CHUNK_BYTES):
            pass


def validate_download_file(path: Path, expected_size: int, expected_sha256: str | None, *, gzip_test: bool) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    observed_size = path.stat().st_size
    observed_sha256 = sha256_file(path)
    if observed_size != expected_size:
        raise RuntimeError(f"size mismatch: observed {observed_size}, expected {expected_size}")
    if expected_sha256 is not None and observed_sha256 != expected_sha256:
        raise RuntimeError(f"sha256 mismatch: observed {observed_sha256}, expected {expected_sha256}")
    if gzip_test:
        gzip_stream_test(path)
    return {
        "observed_size_bytes": int(observed_size),
        "observed_sha256": observed_sha256,
        "gzip_tested": bool(gzip_test),
    }


def invalid_suffix(reason: str) -> str:
    clean = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in reason)[:80]
    return f".invalid.{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.{clean}"


def move_as_invalid(path: Path, reason: str) -> Path:
    destination = path.with_name(path.name + invalid_suffix(reason))
    os.replace(path, destination)
    return destination


def make_read_only(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        path.chmod(mode & ~0o222)
    except Exception:
        pass


def make_writable(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        path.chmod(mode | 0o600)
    except Exception:
        pass


def stream_download_with_resume(target: DownloadTarget, tmp_path: Path) -> dict[str, Any]:
    target.local_path.parent.mkdir(parents=True, exist_ok=True)
    attempts: list[dict[str, Any]] = []
    for attempt in range(1, 41):
        offset = tmp_path.stat().st_size if tmp_path.exists() else 0
        if tmp_path.exists():
            make_writable(tmp_path)
        if offset > target.expected_size_bytes:
            moved = move_as_invalid(tmp_path, "temp_larger_than_expected")
            attempts.append({"attempt": attempt, "status": "restart", "reason": "temp_larger_than_expected", "moved_to": str(moved)})
            offset = 0
        if offset == target.expected_size_bytes:
            break
        headers: dict[str, str] = {}
        mode = "ab"
        if offset:
            headers["Range"] = f"bytes={offset}-"
        else:
            mode = "wb"
        started = time.perf_counter()
        try:
            with url_open(target.official_url, headers=headers, timeout=180) as response:
                status = int(getattr(response, "status", 200))
                if offset and status != 206:
                    offset = 0
                    mode = "wb"
                if status not in {200, 206}:
                    raise RuntimeError(f"GET returned HTTP {status}")
                with tmp_path.open(mode + "") as handle:
                    while True:
                        chunk = response.read(CHUNK_BYTES)
                        if not chunk:
                            break
                        handle.write(chunk)
                    handle.flush()
                    os.fsync(handle.fileno())
            attempts.append(
                {
                    "attempt": attempt,
                    "status": "stream_completed",
                    "http_status": status,
                    "started_offset_bytes": int(offset),
                    "current_size_bytes": int(tmp_path.stat().st_size),
                    "runtime_seconds": round(time.perf_counter() - started, 3),
                }
            )
        except Exception as exc:
            attempts.append(
                {
                    "attempt": attempt,
                    "status": "stream_failed",
                    "started_offset_bytes": int(offset),
                    "current_size_bytes": int(tmp_path.stat().st_size) if tmp_path.exists() else 0,
                    "error": repr(exc),
                    "runtime_seconds": round(time.perf_counter() - started, 3),
                }
            )
            time.sleep(min(10, attempt * 2))
        if tmp_path.exists() and tmp_path.stat().st_size == target.expected_size_bytes:
            break
    if not tmp_path.exists() or tmp_path.stat().st_size != target.expected_size_bytes:
        observed = tmp_path.stat().st_size if tmp_path.exists() else 0
        raise RuntimeError(
            f"download did not reach expected size for {target.filename}: observed {observed}, "
            f"expected {target.expected_size_bytes}"
        )
    return {"attempts": attempts}


def establish_remote_sha256(target: DownloadTarget) -> dict[str, Any]:
    verify_path = target.local_path.with_name(f".{target.local_path.name}.verify.tmp.{os.getpid()}")
    if verify_path.exists():
        move_as_invalid(verify_path, "stale_verify_temp")
    download_info = stream_download_with_resume(target, verify_path)
    validation = validate_download_file(verify_path, target.expected_size_bytes, None, gzip_test=True)
    verify_path.unlink()
    return {
        "expected_sha256": validation["observed_sha256"],
        "remote_verify_download": download_info,
    }


def prepare_resume_temp(target: DownloadTarget, known_sha256: str | None) -> tuple[Path, list[dict[str, Any]], str | None]:
    actions: list[dict[str, Any]] = []
    tmp_path = target.local_path.with_name(f".{target.local_path.name}.download.tmp")
    legacy_part = target.local_path.with_name(target.local_path.name + ".part")

    if tmp_path.exists() and tmp_path.stat().st_size > target.expected_size_bytes:
        moved = move_as_invalid(tmp_path, "stale_download_temp_larger_than_expected")
        actions.append({"action": "moved_stale_temp_invalid", "path": str(moved)})

    if legacy_part.exists():
        if not tmp_path.exists() or legacy_part.stat().st_size > tmp_path.stat().st_size:
            if tmp_path.exists():
                moved = move_as_invalid(tmp_path, "superseded_by_legacy_part")
                actions.append({"action": "moved_smaller_temp_invalid", "path": str(moved)})
            os.replace(legacy_part, tmp_path)
            make_writable(tmp_path)
            actions.append({"action": "adopted_legacy_part_for_resume", "path": str(tmp_path)})
        else:
            moved = move_as_invalid(legacy_part, "superseded_by_download_temp")
            actions.append({"action": "moved_legacy_part_invalid", "path": str(moved)})

    if target.local_path.exists():
        observed_size = target.local_path.stat().st_size
        if observed_size < target.expected_size_bytes:
            if tmp_path.exists() and tmp_path.stat().st_size >= observed_size:
                moved = move_as_invalid(target.local_path, "incomplete_final_superseded_by_temp")
                actions.append({"action": "moved_incomplete_final_invalid", "path": str(moved)})
            else:
                if tmp_path.exists():
                    moved = move_as_invalid(tmp_path, "superseded_by_incomplete_final")
                    actions.append({"action": "moved_temp_invalid", "path": str(moved)})
                os.replace(target.local_path, tmp_path)
                make_writable(tmp_path)
                actions.append({"action": "adopted_incomplete_final_for_resume", "path": str(tmp_path)})
        elif observed_size > target.expected_size_bytes:
            moved = move_as_invalid(target.local_path, "final_larger_than_expected")
            actions.append({"action": "moved_oversized_final_invalid", "path": str(moved)})
        else:
            try:
                validation = validate_download_file(target.local_path, target.expected_size_bytes, known_sha256, gzip_test=True)
                expected_sha256 = known_sha256 or validation["observed_sha256"]
                if known_sha256 is None:
                    remote = establish_remote_sha256(target)
                    expected_sha256 = remote["expected_sha256"]
                    if validation["observed_sha256"] != expected_sha256:
                        moved = move_as_invalid(target.local_path, "checksum_mismatch_against_remote")
                        actions.append(
                            {
                                "action": "moved_checksum_invalid_final",
                                "path": str(moved),
                                "local_sha256": validation["observed_sha256"],
                                "remote_sha256": expected_sha256,
                            }
                        )
                    else:
                        actions.append({"action": "established_sha256_from_official_url", **remote})
                        return tmp_path, actions, expected_sha256
                else:
                    return tmp_path, actions, expected_sha256
            except Exception as exc:
                reason = "same_size_final_validation_failed"
                moved = move_as_invalid(target.local_path, reason)
                actions.append({"action": "moved_invalid_final", "path": str(moved), "error": repr(exc)})
    return tmp_path, actions, known_sha256


def validate_or_download_target(target: DownloadTarget, known_sha256: str | None) -> dict[str, Any]:
    started = time.perf_counter()
    tmp_path, actions, expected_sha256 = prepare_resume_temp(target, known_sha256)
    status = "existing_validated"
    if not target.local_path.exists():
        status = "downloaded"
        if tmp_path.exists() and tmp_path.stat().st_size:
            status = "resumed"
        download_info = stream_download_with_resume(target, tmp_path)
        validation = validate_download_file(tmp_path, target.expected_size_bytes, expected_sha256, gzip_test=True)
        expected_sha256 = expected_sha256 or validation["observed_sha256"]
        os.replace(tmp_path, target.local_path)
        make_read_only(target.local_path)
        final_validation = validate_download_file(target.local_path, target.expected_size_bytes, expected_sha256, gzip_test=True)
        actions.append({"action": "atomic_rename_after_validation", "from": str(tmp_path), "to": str(target.local_path)})
    else:
        final_validation = validate_download_file(target.local_path, target.expected_size_bytes, expected_sha256, gzip_test=True)
        expected_sha256 = expected_sha256 or final_validation["observed_sha256"]
        download_info = {"attempts": []}
        make_read_only(target.local_path)
    return {
        "sample_accession": target.sample_accession,
        "condition": target.condition,
        "file_role": target.role,
        "filename": target.filename,
        "official_url": target.official_url,
        "local_path": str(target.local_path),
        "expected_size_bytes": int(target.expected_size_bytes),
        "expected_sha256": expected_sha256,
        "observed_size_bytes": final_validation["observed_size_bytes"],
        "observed_sha256": final_validation["observed_sha256"],
        "file_type": target.file_type,
        "status": status,
        "actions": actions,
        "download_info": download_info,
        "runtime_seconds": round(time.perf_counter() - started, 3),
    }


def validate_downloads(config: dict[str, Any], data_dir: Path, output_dir: Path) -> dict[str, Any]:
    audit_dir = output_dir / "audit"
    filelist_text = fetch_text(config["geo"]["series_filelist_url"], timeout=60)
    if filelist_text.lstrip().lower().startswith("<html"):
        raise RuntimeError("GEO filelist request returned HTML")
    filelist_rows = parse_geo_filelist(filelist_text)
    targets = expected_targets(config, filelist_rows, data_dir)
    known_sha = load_known_sha256(output_dir)

    rows: list[dict[str, Any]] = []
    for target in targets:
        rows.append(validate_or_download_target(target, known_sha.get(target.filename)))
        atomic_write_json(
            audit_dir / "gse249479_guarded_download_manifest.json",
            {
                "timestamp_utc": now_utc(),
                "note": "NCBI GEO filelist supplies expected byte sizes. SHA-256 values are established from validated official downloads and reused on reruns.",
                "downloads": rows,
            },
        )
        atomic_write_csv(audit_dir / "gse249479_guarded_download_manifest.csv", rows)
    return {
        "timestamp_utc": now_utc(),
        "filelist_url": config["geo"]["series_filelist_url"],
        "filelist_sha256": hashlib.sha256(filelist_text.encode("utf-8")).hexdigest(),
        "downloads": rows,
    }


def read_barcodes(path: Path):
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if frame.shape[1] < 1:
        raise RuntimeError(f"Barcode file has no columns: {path}")
    return frame.iloc[:, 0].astype(str)


def read_features(path: Path):
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if frame.shape[1] < 2:
        raise RuntimeError(f"Feature file must contain at least gene ID and gene symbol columns: {path}")
    columns = ["gene_id", "gene_symbol", "feature_type", "genome"]
    frame = frame.iloc[:, : min(frame.shape[1], len(columns))].copy()
    frame.columns = columns[: frame.shape[1]]
    frame["gene_id"] = frame["gene_id"].astype(str)
    frame["gene_symbol"] = frame["gene_symbol"].astype(str)
    return frame


def unique_feature_index(features) -> tuple[list[str], dict[str, Any]]:
    values = features["gene_id"].astype(str).tolist()
    seen: dict[str, int] = {}
    output: list[str] = []
    duplicates = 0
    for value in values:
        base = value.strip() or "missing_gene_id"
        if base not in seen:
            seen[base] = 0
            output.append(base)
            continue
        duplicates += 1
        seen[base] += 1
        candidate = f"{base}-{seen[base]}"
        while candidate in seen:
            seen[base] += 1
            candidate = f"{base}-{seen[base]}"
        seen[candidate] = 0
        output.append(candidate)
    return output, {
        "preferred_var_name_source": "gene_id",
        "duplicate_gene_id_count": int(duplicates),
        "var_names_modified": bool(duplicates),
    }


def sample_paths_from_config(config: dict[str, Any], data_dir: Path, accession: str) -> dict[str, Path]:
    for sample in config["geo"]["samples"]:
        if sample["accession"] == accession:
            return {
                role: data_dir / "source" / accession / sample["expected_files"][role]
                for role in RNA_FILE_ROLES
            }
    raise KeyError(accession)


def h5_sparse_info(path: Path) -> dict[str, Any]:
    import h5py
    import numpy as np

    with h5py.File(path, "r") as handle:
        value = handle["X"]
        encoding = value.attrs.get("encoding-type", "group")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8", errors="replace")
        shape_attr = value.attrs.get("shape")
        shape = [int(x) for x in shape_attr] if shape_attr is not None else None
        data = value.get("data") if hasattr(value, "get") else None
        indices = value.get("indices") if hasattr(value, "get") else None
        indptr = value.get("indptr") if hasattr(value, "get") else None
        data_bytes = int(data.size * data.dtype.itemsize) if data is not None else None
        indices_bytes = int(indices.size * indices.dtype.itemsize) if indices is not None else None
        indptr_bytes = int(indptr.size * indptr.dtype.itemsize) if indptr is not None else None
        dense_float32_bytes = int(np.prod(shape, dtype=np.int64) * np.dtype("float32").itemsize) if shape else None
        root_provenance = handle.attrs.get("gse249479_guarded_provenance_json")
        if isinstance(root_provenance, bytes):
            root_provenance = root_provenance.decode("utf-8", errors="replace")
        return {
            "encoding_type": encoding,
            "sparse_format": encoding.replace("_matrix", "") if encoding in {"csr_matrix", "csc_matrix"} else None,
            "shape": shape,
            "dtype": str(data.dtype) if data is not None else None,
            "index_dtype": str(indices.dtype) if indices is not None else None,
            "indptr_dtype": str(indptr.dtype) if indptr is not None else None,
            "nnz": int(data.size) if data is not None else None,
            "data_bytes": data_bytes,
            "indices_bytes": indices_bytes,
            "indptr_bytes": indptr_bytes,
            "estimated_sparse_bytes": sum(x for x in [data_bytes, indices_bytes, indptr_bytes] if x is not None),
            "estimated_dense_float32_bytes": dense_float32_bytes,
            "root_write_completed": bool(handle.attrs.get("gse249479_guarded_write_completed", False)),
            "root_provenance_json": root_provenance,
            "top_level_keys": sorted(handle.keys()),
        }


def stamp_h5ad_root_provenance(path: Path, payload: dict[str, Any]) -> None:
    import h5py

    with h5py.File(path, "r+") as handle:
        handle.attrs["gse249479_guarded_write_completed"] = True
        handle.attrs["gse249479_guarded_provenance_json"] = json.dumps(payload, sort_keys=True, default=json_default)


def validate_h5ad(
    path: Path,
    *,
    expected_shape: tuple[int, int] | None = None,
    expected_nnz: int | None = None,
    expected_condition_counts: dict[str, int] | None = None,
    require_write_completed: bool = True,
) -> dict[str, Any]:
    import anndata as ad

    storage = h5_sparse_info(path)
    if storage["sparse_format"] not in {"csr", "csc"}:
        raise RuntimeError(f"X is not stored as CSR/CSC sparse matrix in {path}: {storage['encoding_type']}")
    if expected_shape is not None and tuple(storage["shape"] or []) != tuple(expected_shape):
        raise RuntimeError(f"H5AD shape mismatch for {path}: observed {storage['shape']}, expected {expected_shape}")
    if expected_nnz is not None and int(storage["nnz"] or -1) != int(expected_nnz):
        raise RuntimeError(f"H5AD nnz mismatch for {path}: observed {storage['nnz']}, expected {expected_nnz}")

    adata = ad.read_h5ad(path, backed="r")
    try:
        obs = adata.obs.copy()
        var = adata.var.copy()
        condition_counts = obs["condition"].astype(str).value_counts().sort_index().to_dict() if "condition" in obs else {}
        if expected_condition_counts is not None and condition_counts != expected_condition_counts:
            raise RuntimeError(f"Condition counts mismatch: observed {condition_counts}, expected {expected_condition_counts}")
        uns_provenance = adata.uns.get("gse249479_guarded_reconstruction", {})
        uns_completed = bool(uns_provenance.get("write_completed", False)) if isinstance(uns_provenance, dict) else False
        root_completed = bool(storage.get("root_write_completed", False))
        if require_write_completed and not (uns_completed or root_completed):
            raise RuntimeError(f"H5AD provenance does not mark write completed successfully: {path}")
        identifiers = {
            "obs_names_unique": bool(adata.obs_names.is_unique),
            "var_names_unique": bool(adata.var_names.is_unique),
            "obs_columns": list(obs.columns),
            "var_columns": list(var.columns),
            "obs_name_first": str(adata.obs_names[0]) if adata.n_obs else None,
            "obs_name_last": str(adata.obs_names[-1]) if adata.n_obs else None,
            "var_name_first": str(adata.var_names[0]) if adata.n_vars else None,
            "var_name_last": str(adata.var_names[-1]) if adata.n_vars else None,
        }
    finally:
        adata.file.close()
    if not identifiers["obs_names_unique"] or not identifiers["var_names_unique"]:
        raise RuntimeError(f"H5AD identifiers are not unique: {path}")
    return {
        "path": str(path),
        "file_size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
        "shape": [int(expected_shape[0]), int(expected_shape[1])] if expected_shape else storage["shape"],
        "storage": storage,
        "condition_counts": condition_counts,
        "identifiers": identifiers,
        "write_completed": bool(uns_completed or root_completed),
    }


def h5ad_feature_hashes(path: Path) -> dict[str, str]:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        var = adata.var.copy()
        var_names = [str(x) for x in adata.var_names]
    finally:
        adata.file.close()
    signature_cols = [col for col in ["gene_id", "gene_symbol", "feature_type"] if col in var.columns]
    if signature_cols:
        feature_lines = [
            "\t".join(str(row[col]) for col in signature_cols)
            for _, row in var.iterrows()
        ]
    else:
        feature_lines = var_names
    return {
        "feature_signature_sha256": sha256_text_lines(feature_lines),
        "var_names_sha256": sha256_text_lines(var_names),
    }


def sample_h5ad_path(data_dir: Path, accession: str, condition: str) -> Path:
    return data_dir / "intermediate" / f"{accession}_{condition}_raw.h5ad"


def run_sample_worker(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    report_path = Path(args.report_path).expanduser().resolve()
    try:
        import anndata as ad
        import pandas as pd
        import scipy.io
        import scipy.sparse as sp

        root = Path(args.root).resolve()
        config = load_config(root)
        data_dir = Path(args.data_dir).expanduser().resolve()
        accession = args.accession
        condition = args.condition
        output_h5ad = Path(args.output_h5ad).expanduser().resolve()
        files = sample_paths_from_config(config, data_dir, accession)
        for role, path in files.items():
            if not path.exists():
                raise FileNotFoundError(f"Missing {role} file for {accession}: {path}")

        if output_h5ad.exists():
            try:
                existing = validate_h5ad(output_h5ad, require_write_completed=True)
                feature_hashes = h5ad_feature_hashes(output_h5ad)
                report = {
                    "timestamp_utc": now_utc(),
                    "sample_accession": accession,
                    "condition": condition,
                    "status": "existing_validated",
                    "output_h5ad": str(output_h5ad),
                    "shape": existing["storage"]["shape"],
                    "nnz": existing["storage"]["nnz"],
                    "matrix_dtype": existing["storage"]["dtype"],
                    "sparse_format": existing["storage"]["sparse_format"],
                    "raw_integer_like_counts_retained": bool(str(existing["storage"]["dtype"]).startswith(("int", "uint"))),
                    **feature_hashes,
                    "validation": existing,
                    "current_rss_gb": bytes_to_gb(rss_bytes()),
                    "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
                    "runtime_seconds": round(time.perf_counter() - started, 3),
                }
                atomic_write_json(report_path, report)
                print(json.dumps(report, sort_keys=True), flush=True)
                return 0
            except Exception as exc:
                moved = move_as_invalid(output_h5ad, "existing_sample_h5ad_invalid")
                print(f"[sample-worker] moved invalid existing sample H5AD to {moved}: {exc}", file=sys.stderr, flush=True)

        barcodes = read_barcodes(files["barcodes"])
        features = read_features(files["features"])
        with gzip.open(files["matrix"], "rb") as handle:
            matrix = scipy.io.mmread(handle)
        if not sp.issparse(matrix):
            raise RuntimeError(f"Matrix Market reader returned dense matrix for {files['matrix']}")
        if matrix.shape != (len(features), len(barcodes)):
            raise RuntimeError(
                f"MEX shape mismatch for {accession}: matrix {matrix.shape}, "
                f"features/barcodes {(len(features), len(barcodes))}"
            )
        cell_by_gene = matrix.T.tocsr(copy=False)
        cell_by_gene.sort_indices()
        if not sp.isspmatrix_csr(cell_by_gene):
            raise RuntimeError("Sample matrix is not CSR after transpose")

        var_names, var_handling = unique_feature_index(features)
        var = features.copy()
        var.index = pd.Index(var_names, name="gene_id")
        cell_ids = pd.Index([f"{accession}_{barcode}" for barcode in barcodes.astype(str)], name="cell_id")
        obs = pd.DataFrame(
            {
                "condition": condition,
                "sample_accession": accession,
                "library_id": accession,
                "source_dataset": config["dataset"]["accession"],
                "barcode_original": barcodes.astype(str).to_numpy(copy=False),
            },
            index=cell_ids,
        )
        provenance = {
            "write_completed": True,
            "timestamp_utc": now_utc(),
            "sample_accession": accession,
            "condition": condition,
            "source_files": {
                role: {
                    "path": str(path),
                    "size_bytes": int(path.stat().st_size),
                    "sha256": sha256_file(path),
                }
                for role, path in files.items()
            },
            "raw_counts_preserved": True,
            "full_matrix_densified": False,
            "sparse_format": cell_by_gene.getformat(),
            "matrix_dtype": str(cell_by_gene.dtype),
            "var_name_handling": var_handling,
        }
        adata = ad.AnnData(X=cell_by_gene, obs=obs, var=var)
        adata.uns["gse249479_guarded_reconstruction"] = provenance

        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        tmp_h5ad = output_h5ad.with_name(f".{output_h5ad.name}.tmp.{os.getpid()}")
        if tmp_h5ad.exists():
            move_as_invalid(tmp_h5ad, "stale_sample_temp")
        adata.write_h5ad(tmp_h5ad)
        stamp_h5ad_root_provenance(tmp_h5ad, provenance)
        expected_shape = (int(cell_by_gene.shape[0]), int(cell_by_gene.shape[1]))
        expected_nnz = int(cell_by_gene.nnz)
        tmp_validation = validate_h5ad(
            tmp_h5ad,
            expected_shape=expected_shape,
            expected_nnz=expected_nnz,
            expected_condition_counts={condition: expected_shape[0]},
            require_write_completed=True,
        )
        os.replace(tmp_h5ad, output_h5ad)
        final_validation = validate_h5ad(
            output_h5ad,
            expected_shape=expected_shape,
            expected_nnz=expected_nnz,
            expected_condition_counts={condition: expected_shape[0]},
            require_write_completed=True,
        )
        feature_signature_sha256 = sha256_text_lines(
            [
                "\t".join(str(row[col]) for col in [c for c in ["gene_id", "gene_symbol", "feature_type"] if c in features.columns])
                for _, row in features.iterrows()
            ]
        )
        report = {
            "timestamp_utc": now_utc(),
            "sample_accession": accession,
            "condition": condition,
            "status": "completed",
            "output_h5ad": str(output_h5ad),
            "shape": list(expected_shape),
            "nnz": expected_nnz,
            "matrix_dtype": str(cell_by_gene.dtype),
            "sparse_format": cell_by_gene.getformat(),
            "raw_integer_like_counts_retained": bool(str(cell_by_gene.dtype).startswith(("int", "uint"))),
            "feature_signature_sha256": feature_signature_sha256,
            "var_names_sha256": sha256_text_lines([str(x) for x in var.index.tolist()]),
            "tmp_validation": tmp_validation,
            "validation": final_validation,
            "current_rss_gb": bytes_to_gb(rss_bytes()),
            "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
            "runtime_seconds": round(time.perf_counter() - started, 3),
        }
        atomic_write_json(report_path, report)
        print(json.dumps(report, sort_keys=True), flush=True)
        return 0
    except Exception as exc:
        report = {
            "timestamp_utc": now_utc(),
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
            "current_rss_gb": bytes_to_gb(rss_bytes()),
            "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
            "runtime_seconds": round(time.perf_counter() - started, 3),
        }
        atomic_write_json(report_path, report)
        print(json.dumps(report, sort_keys=True), file=sys.stderr, flush=True)
        return 1


def classify_returncode(returncode: int) -> dict[str, Any]:
    if returncode >= 0:
        return {"classification": "process_exit", "returncode": returncode}
    signum = -returncode
    try:
        name = signal.Signals(signum).name
    except ValueError:
        name = f"SIG{signum}"
    if signum == signal.SIGBUS:
        cause = "native-library SIGBUS"
    elif signum == signal.SIGKILL:
        cause = "possible OOM or external kill"
    else:
        cause = f"terminated by {name}"
        return {"classification": cause, "returncode": returncode, "signal": signum, "signal_name": name}


def run_worker_subprocess(command: list[str], report_path: Path) -> dict[str, Any]:
    started = time.perf_counter()
    completed = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    payload = {
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout[-4000:],
        "stderr": completed.stderr[-4000:],
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "exit_classification": classify_returncode(completed.returncode),
        "report_path": str(report_path),
    }
    if report_path.exists():
        try:
            payload["worker_report"] = json.loads(report_path.read_text(encoding="utf-8"))
        except Exception as exc:
            payload["worker_report_error"] = repr(exc)
    if completed.returncode != 0:
        raise RuntimeError(json.dumps(payload, indent=2, sort_keys=True, default=json_default))
    return payload


def reconstruct_samples(config: dict[str, Any], root: Path, data_dir: Path, output_dir: Path) -> dict[str, Any]:
    audit_dir = output_dir / "audit"
    sample_runs: list[dict[str, Any]] = []
    for sample in config["geo"]["samples"]:
        accession = sample["accession"]
        condition = sample["condition"]
        output_h5ad = sample_h5ad_path(data_dir, accession, condition)
        report_path = audit_dir / f"gse249479_{accession}_{condition}_sample_reconstruction.json"
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            "sample",
            "--root",
            str(root),
            "--data-dir",
            str(data_dir),
            "--accession",
            accession,
            "--condition",
            condition,
            "--output-h5ad",
            str(output_h5ad),
            "--report-path",
            str(report_path),
        ]
        sample_runs.append(run_worker_subprocess(command, report_path))
        atomic_write_json(
            audit_dir / "gse249479_guarded_sample_reconstruction_summary.json",
            {"timestamp_utc": now_utc(), "sample_runs": sample_runs},
        )
    return {"timestamp_utc": now_utc(), "sample_runs": sample_runs}


def sample_reports_from_runs(sample_runs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    reports = []
    for run in sample_runs:
        report = run.get("worker_report")
        if not isinstance(report, dict):
            raise RuntimeError(f"Worker report missing for {run.get('command')}")
        if report.get("status") not in {"completed", "existing_validated"}:
            raise RuntimeError(f"Sample report is not complete: {report}")
        reports.append(report)
    return reports


def compare_feature_axes(sample_reports: list[dict[str, Any]]) -> dict[str, Any]:
    signatures = {report["sample_accession"]: report.get("feature_signature_sha256") for report in sample_reports}
    var_signatures = {report["sample_accession"]: report.get("var_names_sha256") for report in sample_reports}
    if len(set(signatures.values())) != 1 or len(set(var_signatures.values())) != 1:
        raise RuntimeError(f"Feature axes differ across samples: {signatures}, {var_signatures}")
    n_vars = {report["sample_accession"]: int(report["shape"][1]) for report in sample_reports}
    if len(set(n_vars.values())) != 1:
        raise RuntimeError(f"Number of vars differs across samples: {n_vars}")
    return {
        "feature_signature_sha256_by_sample": signatures,
        "var_names_sha256_by_sample": var_signatures,
        "n_vars_by_sample": n_vars,
        "axis_status": "identical_gene_ids_order_and_feature_signature",
    }


def estimate_concat_memory(sample_reports: list[dict[str, Any]]) -> dict[str, Any]:
    total_sparse_bytes = 0
    total_nnz = 0
    total_obs = 0
    n_vars = None
    for report in sample_reports:
        storage = report["validation"]["storage"]
        total_sparse_bytes += int(storage.get("estimated_sparse_bytes") or 0)
        total_nnz += int(storage.get("nnz") or 0)
        total_obs += int(report["shape"][0])
        n_vars = int(report["shape"][1])
    # Conservative multiplier for fallback sparse concatenation and H5AD write buffers.
    estimated_peak_bytes = int(total_sparse_bytes * 3.0 + 1.0 * 1024**3)
    return {
        "total_sparse_bytes": total_sparse_bytes,
        "total_sparse_gb": bytes_to_gb(total_sparse_bytes),
        "total_nnz": total_nnz,
        "total_obs": total_obs,
        "n_vars": n_vars,
        "estimated_peak_rss_bytes": estimated_peak_bytes,
        "estimated_peak_rss_gb": bytes_to_gb(estimated_peak_bytes),
        "stop_threshold_gb": 10.0,
        "status": "within_guard" if estimated_peak_bytes <= 10 * 1024**3 else "stop_exceeds_guard",
    }


def final_output_path(data_dir: Path) -> Path:
    return data_dir / "zeng2026_xenograft_cd34_rna_raw.h5ad"


def run_concat_worker(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    report_path = Path(args.report_path).expanduser().resolve()
    try:
        import anndata as ad

        root = Path(args.root).resolve()
        data_dir = Path(args.data_dir).expanduser().resolve()
        output_h5ad = Path(args.output_h5ad).expanduser().resolve()
        sample_reports = json.loads(Path(args.sample_reports_json).read_text(encoding="utf-8"))
        axis_check = compare_feature_axes(sample_reports)
        estimate = estimate_concat_memory(sample_reports)
        if estimate["status"] != "within_guard":
            raise RuntimeError(f"Concatenation estimate exceeds 10 GB guard: {estimate}")
        expected_shape = (int(estimate["total_obs"]), int(estimate["n_vars"]))
        expected_nnz = int(estimate["total_nnz"])
        expected_condition_counts = {
            str(report["condition"]): int(report["shape"][0])
            for report in sample_reports
        }
        if output_h5ad.exists():
            try:
                existing = validate_h5ad(
                    output_h5ad,
                    expected_shape=expected_shape,
                    expected_nnz=expected_nnz,
                    expected_condition_counts=dict(sorted(expected_condition_counts.items())),
                    require_write_completed=True,
                )
                report = {
                    "timestamp_utc": now_utc(),
                    "status": "existing_validated",
                    "output_h5ad": str(output_h5ad),
                    "axis_check": axis_check,
                    "memory_estimate": estimate,
                    "validation": existing,
                    "current_rss_gb": bytes_to_gb(rss_bytes()),
                    "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
                    "runtime_seconds": round(time.perf_counter() - started, 3),
                }
                atomic_write_json(report_path, report)
                print(json.dumps(report, sort_keys=True), flush=True)
                return 0
            except Exception as exc:
                moved = move_as_invalid(output_h5ad, "existing_combined_h5ad_invalid")
                print(f"[concat-worker] moved invalid existing combined H5AD to {moved}: {exc}", file=sys.stderr, flush=True)

        sample_paths = [report["output_h5ad"] for report in sample_reports]
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        tmp_h5ad = output_h5ad.with_name(f".{output_h5ad.name}.tmp.{os.getpid()}")
        if tmp_h5ad.exists():
            move_as_invalid(tmp_h5ad, "stale_concat_temp")

        method = "anndata.experimental.concat_on_disk"
        concat_on_disk = getattr(getattr(ad, "experimental", None), "concat_on_disk", None)
        if concat_on_disk is not None:
            concat_on_disk(sample_paths, tmp_h5ad, axis=0, join="inner", merge="same", uns_merge=None, index_unique=None)
        else:
            method = "manual_sparse_vstack_fallback"
            import h5py
            import pandas as pd
            import scipy.sparse as sp

            matrices = []
            obs_frames = []
            var_frame = None
            for sample_path in sample_paths:
                backed = ad.read_h5ad(sample_path, backed="r")
                try:
                    obs_frames.append(backed.obs.copy())
                    if var_frame is None:
                        var_frame = backed.var.copy()
                finally:
                    backed.file.close()
                with h5py.File(sample_path, "r") as handle:
                    x = handle["X"]
                    shape = tuple(int(v) for v in x.attrs["shape"])
                    matrices.append(
                        sp.csr_matrix(
                            (x["data"][:], x["indices"][:], x["indptr"][:]),
                            shape=shape,
                        )
                    )
            X = sp.vstack(matrices, format="csr")
            obs = pd.concat(obs_frames, axis=0)
            out = ad.AnnData(X=X, obs=obs, var=var_frame)
            out.write_h5ad(tmp_h5ad)
            del out, X, matrices

        provenance = {
            "write_completed": True,
            "timestamp_utc": now_utc(),
            "method": method,
            "sample_h5ads": sample_paths,
            "axis_check": axis_check,
            "memory_estimate": estimate,
            "full_matrix_densified": False,
            "forbidden_downstream_steps_run": [],
            "root": str(root),
        }
        stamp_h5ad_root_provenance(tmp_h5ad, provenance)
        tmp_validation = validate_h5ad(
            tmp_h5ad,
            expected_shape=expected_shape,
            expected_nnz=expected_nnz,
            expected_condition_counts=dict(sorted(expected_condition_counts.items())),
            require_write_completed=True,
        )
        os.replace(tmp_h5ad, output_h5ad)
        final_validation = validate_h5ad(
            output_h5ad,
            expected_shape=expected_shape,
            expected_nnz=expected_nnz,
            expected_condition_counts=dict(sorted(expected_condition_counts.items())),
            require_write_completed=True,
        )
        report = {
            "timestamp_utc": now_utc(),
            "status": "completed",
            "method": method,
            "output_h5ad": str(output_h5ad),
            "axis_check": axis_check,
            "memory_estimate": estimate,
            "tmp_validation": tmp_validation,
            "validation": final_validation,
            "current_rss_gb": bytes_to_gb(rss_bytes()),
            "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
            "runtime_seconds": round(time.perf_counter() - started, 3),
        }
        atomic_write_json(report_path, report)
        complete_payload = {
            "timestamp_utc": now_utc(),
            "status": "reconstruction_complete",
            "output_h5ad": str(output_h5ad),
            "sha256": final_validation["sha256"],
            "shape": final_validation["shape"],
            "condition_counts": final_validation["condition_counts"],
            "sparse_format": final_validation["storage"]["sparse_format"],
            "nnz": final_validation["storage"]["nnz"],
        }
        atomic_write_json(output_h5ad.with_suffix(output_h5ad.suffix + ".complete.json"), complete_payload)
        print(json.dumps(report, sort_keys=True), flush=True)
        return 0
    except Exception as exc:
        report = {
            "timestamp_utc": now_utc(),
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
            "current_rss_gb": bytes_to_gb(rss_bytes()),
            "peak_rss_gb": bytes_to_gb(peak_rss_bytes()),
            "runtime_seconds": round(time.perf_counter() - started, 3),
        }
        atomic_write_json(report_path, report)
        print(json.dumps(report, sort_keys=True), file=sys.stderr, flush=True)
        return 1


def guarded_concat(root: Path, data_dir: Path, output_dir: Path, sample_runs: list[dict[str, Any]]) -> dict[str, Any]:
    audit_dir = output_dir / "audit"
    reports = sample_reports_from_runs(sample_runs)
    sample_reports_json = audit_dir / "gse249479_guarded_sample_reports_for_concat.json"
    # Keep the worker input as a bare list to avoid accidental schema ambiguity.
    atomic_write_text(sample_reports_json, json.dumps(reports, indent=2, sort_keys=True, default=json_default) + "\n")
    report_path = audit_dir / "gse249479_guarded_concat_reconstruction.json"
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "concat",
        "--root",
        str(root),
        "--data-dir",
        str(data_dir),
        "--output-h5ad",
        str(final_output_path(data_dir)),
        "--sample-reports-json",
        str(sample_reports_json),
        "--report-path",
        str(report_path),
    ]
    return run_worker_subprocess(command, report_path)


def failure_classification(download_summary: dict[str, Any] | None, sample_summary: dict[str, Any] | None, concat_summary: dict[str, Any] | None, error: Exception | None) -> dict[str, Any]:
    if error is None:
        return {
            "previous_failure_operation": "not_inferred_from_available_prior_logs",
            "confirmed_oom": False,
            "file_truncation_or_corruption": "incomplete_downloads_confirmed_before_resume",
            "disk_or_swap_problem": False,
            "native_library_sigbus": False,
            "unknown_cause": True,
            "notes": "The prior Bus error was not accompanied by a core file or 00a execution report in the inspected paths. Existing source files showed incomplete downloads, but that does not prove the Bus error cause.",
        }
    text = repr(error)
    if "SIGBUS" in text or "native-library SIGBUS" in text:
        return {
            "previous_failure_operation": "current_resume_subprocess_failed_with_SIGBUS",
            "confirmed_oom": False,
            "file_truncation_or_corruption": "unknown",
            "disk_or_swap_problem": False,
            "native_library_sigbus": True,
            "unknown_cause": False,
            "notes": text[-1000:],
        }
    if "No space left" in text:
        return {
            "previous_failure_operation": "current_resume_write_or_download",
            "confirmed_oom": False,
            "file_truncation_or_corruption": "unknown",
            "disk_or_swap_problem": True,
            "native_library_sigbus": False,
            "unknown_cause": False,
            "notes": text[-1000:],
        }
    if "did not reach expected size" in text or "size mismatch" in text:
        return {
            "previous_failure_operation": "current_resume_download_validation",
            "confirmed_oom": False,
            "file_truncation_or_corruption": "confirmed_incomplete_or_size_invalid_download",
            "disk_or_swap_problem": False,
            "native_library_sigbus": False,
            "unknown_cause": False,
            "notes": text[-1000:],
        }
    if "exceeds 10 GB guard" in text or "SIGKILL" in text:
        return {
            "previous_failure_operation": "current_resume_memory_guard_or_external_kill",
            "confirmed_oom": "not_confirmed" if "SIGKILL" in text else False,
            "file_truncation_or_corruption": "unknown",
            "disk_or_swap_problem": False,
            "native_library_sigbus": False,
            "unknown_cause": "SIGKILL" in text,
            "notes": text[-1000:],
        }
    return {
        "previous_failure_operation": "current_resume_failed_before_completion",
        "confirmed_oom": False,
        "file_truncation_or_corruption": "unknown",
        "disk_or_swap_problem": False,
        "native_library_sigbus": False,
        "unknown_cause": True,
        "notes": text[-1000:],
    }


def run_orchestrator(args: argparse.Namespace) -> int:
    root = Path(args.root).resolve()
    config = load_config(root)
    data_dir = resolve_path(root, args.data_dir)
    output_dir = resolve_path(root, args.output_dir)
    audit_dir = output_dir / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    report_path = audit_dir / REPORT_NAME

    report: dict[str, Any] = {
        "started_at_utc": now_utc(),
        "root": str(root),
        "data_dir": str(data_dir),
        "output_dir": str(output_dir),
        "phases_requested": args.phase,
        "status": "running",
    }
    error: Exception | None = None
    download_summary = None
    sample_summary = None
    concat_summary = None
    try:
        inventory = collect_environment_inventory(data_dir, output_dir)
        report["phase1_inventory"] = inventory
        atomic_write_json(audit_dir / "gse249479_guarded_phase1_inventory.json", inventory)
        atomic_write_json(report_path, report)

        if args.phase in {"all", "downloads"}:
            download_summary = validate_downloads(config, data_dir, output_dir)
            report["phase1_download_validation"] = download_summary
            atomic_write_json(report_path, report)

        if args.phase in {"all", "reconstruct"}:
            if download_summary is None:
                download_summary = validate_downloads(config, data_dir, output_dir)
                report["phase1_download_validation"] = download_summary
                atomic_write_json(report_path, report)
            sample_summary = reconstruct_samples(config, root, data_dir, output_dir)
            report["phase2_sample_reconstruction"] = sample_summary
            atomic_write_json(report_path, report)
            concat_summary = guarded_concat(root, data_dir, output_dir, sample_summary["sample_runs"])
            report["phase3_guarded_concatenation"] = concat_summary
            atomic_write_json(
                audit_dir / "gse249479_guarded_reconstruction_complete.json",
                concat_summary.get("worker_report", {}),
            )

        report["status"] = "completed"
        report["finished_at_utc"] = now_utc()
        report["failure_classification"] = failure_classification(download_summary, sample_summary, concat_summary, None)
        atomic_write_json(report_path, report)
        return 0
    except Exception as exc:
        error = exc
        report["status"] = "failed"
        report["finished_at_utc"] = now_utc()
        report["error"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc()
        report["failure_classification"] = failure_classification(download_summary, sample_summary, concat_summary, error)
        atomic_write_json(report_path, report)
        print(json.dumps(report["failure_classification"], indent=2, sort_keys=True), file=sys.stderr, flush=True)
        return 1


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    config = load_config(root)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(root))
    parser.add_argument("--data-dir", default=os.environ.get(config["data_dir_env"], config["default_data_dir"]))
    parser.add_argument("--output-dir", default=os.environ.get(config["output_dir_env"], config["default_output_dir"]))
    parser.add_argument("--phase", choices=["all", "downloads", "reconstruct"], default="all")
    parser.add_argument("--worker", choices=["sample", "concat"], default=None)
    parser.add_argument("--accession")
    parser.add_argument("--condition")
    parser.add_argument("--output-h5ad")
    parser.add_argument("--report-path")
    parser.add_argument("--sample-reports-json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.worker == "sample":
        missing = [name for name in ["accession", "condition", "output_h5ad", "report_path"] if getattr(args, name) is None]
        if missing:
            raise SystemExit(f"Missing sample worker arguments: {', '.join(missing)}")
        return run_sample_worker(args)
    if args.worker == "concat":
        missing = [name for name in ["output_h5ad", "report_path", "sample_reports_json"] if getattr(args, name) is None]
        if missing:
            raise SystemExit(f"Missing concat worker arguments: {', '.join(missing)}")
        return run_concat_worker(args)
    return run_orchestrator(args)


if __name__ == "__main__":
    raise SystemExit(main())
