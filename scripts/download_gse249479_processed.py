#!/usr/bin/env python3
"""Download the public GSE249479 processed RNA MEX files needed for Phase 1A."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import stat
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from gse249479_memory_safe import (
    MemoryAudit,
    bytes_to_gb,
    ensure_output_tree,
    load_config,
    memory_threshold_bytes,
    configure_temp_environment,
    relative_or_absolute,
    require_repo_local_dataset_paths,
    repo_root,
    require_active_branch,
    resolve_path,
    sha256_file,
    write_json,
)


USER_AGENT = "scgeo-gse249479-phase1a/1.0"
CHUNK_BYTES = 1024 * 1024
RNA_FILE_ROLES = ("barcodes", "features", "matrix")


@dataclass(frozen=True)
class DownloadTarget:
    sample_accession: str
    condition: str
    role: str
    filename: str
    official_url: str
    expected_size_bytes: int | None
    file_type: str | None
    local_path: Path


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sample_bucket(accession: str) -> str:
    match = re.fullmatch(r"([A-Za-z]+)(\d+)", accession)
    if not match:
        raise ValueError(f"Unsupported GEO accession format: {accession}")
    prefix, digits = match.groups()
    return f"{prefix}{digits[:-3]}nnn"


def url_open(url: str, *, method: str = "GET", headers: dict[str, str] | None = None, timeout: int = 60):
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)
    request = urllib.request.Request(url, headers=request_headers, method=method)
    return urllib.request.urlopen(request, timeout=timeout)


def fetch_text(url: str, *, timeout: int = 60) -> str:
    with url_open(url, timeout=timeout) as response:
        payload = response.read()
    return payload.decode("utf-8", errors="replace")


def parse_geo_filelist(text: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        name = (row.get("Name") or "").strip()
        if not name:
            continue
        rows.append(
            {
                "record_kind": (row.get("#Archive/File") or "").strip(),
                "filename": name,
                "modified_time_geo": (row.get("Time") or "").strip(),
                "size_bytes": int(row["Size"]) if (row.get("Size") or "").strip().isdigit() else None,
                "file_type": (row.get("Type") or "").strip() or None,
            }
        )
    return pd.DataFrame(rows)


def expected_file_records(config: dict[str, Any], filelist: pd.DataFrame, data_dir: Path) -> list[DownloadTarget]:
    targets: list[DownloadTarget] = []
    by_name = {
        str(row.filename): row
        for row in filelist.itertuples(index=False)
        if str(row.record_kind).lower() == "file"
    }
    template = config["geo"]["sample_url_template"]
    for sample in config["geo"]["samples"]:
        accession = sample["accession"]
        for role in RNA_FILE_ROLES:
            filename = sample["expected_files"][role]
            if filename not in by_name:
                raise RuntimeError(f"Expected GEO file is absent from filelist.txt: {filename}")
            record = by_name[filename]
            url = template.format(
                sample_bucket=sample_bucket(accession),
                accession=accession,
                filename=urllib.parse.quote(filename),
            )
            targets.append(
                DownloadTarget(
                    sample_accession=accession,
                    condition=sample["condition"],
                    role=role,
                    filename=filename,
                    official_url=url,
                    expected_size_bytes=int(record.size_bytes) if pd.notna(record.size_bytes) else None,
                    file_type=record.file_type,
                    local_path=data_dir / "source" / accession / filename,
                )
            )
    return targets


def remote_head(url: str) -> dict[str, Any]:
    try:
        with url_open(url, method="HEAD", timeout=60) as response:
            headers = dict(response.headers.items())
            return {
                "status": int(getattr(response, "status", 200)),
                "content_length": int(headers["Content-Length"]) if headers.get("Content-Length") else None,
                "content_type": headers.get("Content-Type"),
                "last_modified": headers.get("Last-Modified"),
                "accept_ranges": headers.get("Accept-Ranges"),
            }
    except urllib.error.HTTPError as exc:
        return {
            "status": int(exc.code),
            "content_length": None,
            "content_type": exc.headers.get("Content-Type") if exc.headers else None,
            "last_modified": exc.headers.get("Last-Modified") if exc.headers else None,
            "accept_ranges": exc.headers.get("Accept-Ranges") if exc.headers else None,
        }


def looks_like_html(path: Path) -> bool:
    with path.open("rb") as handle:
        head = handle.read(512).lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<title>404" in head


def validate_download(path: Path, expected_size: int | None) -> None:
    if not path.exists():
        raise FileNotFoundError(path)
    observed_size = path.stat().st_size
    if observed_size <= 0:
        raise RuntimeError(f"Downloaded file is empty: {path}")
    if expected_size is not None and observed_size != expected_size:
        raise RuntimeError(
            f"Downloaded file size mismatch for {path.name}: observed {observed_size}, expected {expected_size}"
        )
    if looks_like_html(path):
        raise RuntimeError(f"Downloaded file looks like an HTML error page: {path}")
    if path.suffix == ".gz":
        with path.open("rb") as handle:
            magic = handle.read(2)
        if magic != b"\x1f\x8b":
            raise RuntimeError(f"Downloaded gzip file has invalid magic bytes: {path}")


def make_read_only(path: Path) -> None:
    current_mode = path.stat().st_mode
    path.chmod(current_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def download_one(target: DownloadTarget) -> dict[str, Any]:
    target.local_path.parent.mkdir(parents=True, exist_ok=True)
    head = remote_head(target.official_url)
    if head["status"] != 200:
        raise RuntimeError(f"Official GEO URL returned HTTP {head['status']}: {target.official_url}")
    if target.expected_size_bytes is not None and head["content_length"] not in {None, target.expected_size_bytes}:
        raise RuntimeError(
            f"Remote size mismatch for {target.filename}: HEAD {head['content_length']}, "
            f"filelist {target.expected_size_bytes}"
        )

    status = "existing_validated"
    if target.local_path.exists():
        validate_download(target.local_path, target.expected_size_bytes)
    else:
        status = "downloaded"
        part_path = target.local_path.with_name(f"{target.local_path.name}.part")
        offset = part_path.stat().st_size if part_path.exists() else 0
        headers: dict[str, str] = {}
        mode = "ab" if offset else "wb"
        if offset:
            headers["Range"] = f"bytes={offset}-"
        try:
            with url_open(target.official_url, headers=headers, timeout=120) as response:
                response_status = int(getattr(response, "status", 200))
                if offset and response_status != 206:
                    offset = 0
                    mode = "wb"
                with part_path.open(mode + "") as handle:
                    while True:
                        chunk = response.read(CHUNK_BYTES)
                        if not chunk:
                            break
                        handle.write(chunk)
        except Exception:
            raise
        validate_download(part_path, target.expected_size_bytes)
        part_path.replace(target.local_path)
        validate_download(target.local_path, target.expected_size_bytes)
    make_read_only(target.local_path)
    stat_result = target.local_path.stat()
    return {
        "sample_accession": target.sample_accession,
        "condition": target.condition,
        "file_role": target.role,
        "filename": target.filename,
        "official_url": target.official_url,
        "local_path": str(target.local_path),
        "expected_size_bytes": target.expected_size_bytes,
        "observed_size_bytes": int(stat_result.st_size),
        "content_type": head["content_type"],
        "remote_last_modified": head["last_modified"],
        "sha256": sha256_file(target.local_path),
        "download_status": status,
        "downloaded_at_utc": now_utc(),
        "immutable_source": True,
    }


def sample_record_text(accession: str) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    try:
        text = fetch_text(url, timeout=60)
    except Exception as exc:  # noqa: BLE001 - metadata inventory should remain available.
        return f"unavailable: {exc!r}"
    text = re.sub(r"<script.*?</script>", " ", text, flags=re.IGNORECASE | re.DOTALL)
    text = re.sub(r"<style.*?</style>", " ", text, flags=re.IGNORECASE | re.DOTALL)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text)
    return re.sub(r"\s+", " ", text).strip()


def classify_metadata_signal(name: str, text: str = "") -> dict[str, str]:
    combined = f"{name} {text}".lower()
    return {
        "donor_identity_status": "not_observed_in_public_inventory",
        "xenograft_recipient_status": (
            "xenograft_context_mentioned_recipient_identity_not_confirmed"
            if "xenograft" in combined
            else "not_observed_in_public_inventory"
        ),
        "cord_blood_pool_status": (
            "cord_blood_context_mentioned_pool_identity_not_confirmed"
            if "cord blood" in combined or "cord-blood" in combined
            else "not_observed_in_public_inventory"
        ),
        "souporcell_genetic_clade_status": (
            "souporcell_or_genetic_clade_metadata_mentioned"
            if "soupor" in combined or "genetic clade" in combined or "genotyping" in combined
            else "not_observed_in_public_inventory"
        ),
        "original_cell_annotations_status": (
            "possible_in_seurat_or_annotation_file_not_loaded"
            if "seurat" in combined or "annotation" in combined or "metadata" in combined
            else "not_observed_in_public_inventory"
        ),
    }


def metadata_inventory(config: dict[str, Any], filelist: pd.DataFrame) -> pd.DataFrame:
    keywords = [str(x).lower() for x in config["geo"]["metadata_inventory_keywords"]]
    rows: list[dict[str, Any]] = []
    series_url = config["geo"]["series_supplementary_url"]
    for record in filelist.itertuples(index=False):
        filename = str(record.filename)
        low = filename.lower()
        if not any(keyword in low for keyword in keywords):
            continue
        signals = classify_metadata_signal(filename)
        rows.append(
            {
                "inventory_source": "GEO_series_filelist",
                "accession": None,
                "filename": filename,
                "official_url": urllib.parse.urljoin(series_url, urllib.parse.quote(filename)),
                "file_type": record.file_type,
                "size_bytes": int(record.size_bytes) if pd.notna(record.size_bytes) else None,
                "category": "public_supplementary_metadata_or_container",
                "load_status": "inventoried_not_loaded",
                "notes": "RDS/genotyping/annotation-like public file inventoried only; expression matrices were not converted from RDS.",
                **signals,
            }
        )
    for sample in config["geo"]["samples"]:
        accession = sample["accession"]
        text = sample_record_text(accession)
        signals = classify_metadata_signal(accession, text)
        rows.append(
            {
                "inventory_source": "GEO_sample_record_text",
                "accession": accession,
                "filename": None,
                "official_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "file_type": "HTML",
                "size_bytes": None,
                "category": "sample_record_metadata_text",
                "load_status": "record_text_inspected_no_supplementary_container_loaded",
                "notes": text[:500],
                **signals,
            }
        )
    return pd.DataFrame(rows)


def write_checksums(manifest: pd.DataFrame, output_dir: Path, data_dir: Path) -> Path:
    checksum_path = output_dir / "audit" / "checksums.sha256"
    lines = []
    for row in manifest.sort_values(["sample_accession", "file_role"]).itertuples(index=False):
        local_path = Path(row.local_path)
        lines.append(f"{row.sha256}  {relative_or_absolute(local_path, data_dir)}")
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return checksum_path


def run(
    *,
    root: Path,
    data_dir: Path,
    output_dir: Path,
    memory_threshold_gb: float | None = None,
) -> dict[str, Any]:
    config = load_config(root)
    configure_temp_environment(root, config, create=True)
    require_repo_local_dataset_paths(root, data_dir)
    require_active_branch(root, config["required_git_branch"])
    ensure_output_tree(output_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    threshold_bytes = int((memory_threshold_gb or float(config["default_memory_threshold_gb"])) * 1024**3)
    audit = MemoryAudit(threshold_bytes)

    audit_dir = output_dir / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    with audit.section("discover_geo_filelist"):
        filelist_text = fetch_text(config["geo"]["series_filelist_url"], timeout=60)
        if filelist_text.lstrip().lower().startswith("<!doctype") or filelist_text.lstrip().lower().startswith("<html"):
            raise RuntimeError("GEO filelist request returned HTML rather than filelist.txt")
        filelist = parse_geo_filelist(filelist_text)
        targets = expected_file_records(config, filelist, data_dir)

    rows: list[dict[str, Any]] = []
    with audit.section("download_processed_rna_mex_files"):
        for target in targets:
            rows.append(download_one(target))
            time.sleep(0.2)
    manifest = pd.DataFrame(rows)
    manifest_path = audit_dir / "download_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    checksum_path = write_checksums(manifest, output_dir, data_dir)

    with audit.section("inventory_public_metadata"):
        inventory = metadata_inventory(config, filelist)
        inventory_path = audit_dir / "public_metadata_inventory.csv"
        inventory.to_csv(inventory_path, index=False)

    memory_path = audit_dir / "reconstruction_memory.csv"
    memory = audit.dataframe()
    memory.insert(0, "script", Path(__file__).name)
    memory.to_csv(memory_path, index=False)

    summary = {
        "timestamp_utc": now_utc(),
        "data_dir": str(data_dir),
        "manifest_path": str(manifest_path),
        "checksum_path": str(checksum_path),
        "metadata_inventory_path": str(inventory_path),
        "n_downloaded_files": int(len(manifest)),
        "downloaded_bytes": int(manifest["observed_size_bytes"].sum()),
        "peak_rss_gb": bytes_to_gb(audit.peak_rss_bytes),
        "official_geo_sources": {
            "series_filelist_url": config["geo"]["series_filelist_url"],
            "sample_url_template": config["geo"]["sample_url_template"],
        },
    }
    write_json(audit_dir / "download_summary.json", summary)
    return summary


def parse_args() -> argparse.Namespace:
    root = repo_root()
    config = load_config(root)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default=os.environ.get(config["data_dir_env"], config["default_data_dir"]),
        help="Directory for public GEO downloads.",
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get(config["output_dir_env"], config["default_output_dir"]),
        help="Directory for audit artifacts.",
    )
    parser.add_argument(
        "--memory-threshold-gb",
        type=float,
        default=float(os.environ.get(config["memory_threshold_gb_env"], config["default_memory_threshold_gb"])),
        help="RSS threshold in GB.",
    )
    return parser.parse_args()


def main() -> int:
    root = repo_root()
    args = parse_args()
    summary = run(
        root=root,
        data_dir=resolve_path(root, args.data_dir),
        output_dir=resolve_path(root, args.output_dir),
        memory_threshold_gb=float(args.memory_threshold_gb),
    )
    json.dump(summary, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
