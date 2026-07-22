#!/usr/bin/env python3
"""Download only official processed GSE211713 MEX files.

`--all` deliberately launches one fresh subprocess per GSM. Existing files are
reused only after exact-size validation, and every file receives SHA-256.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "configs/gse211713_dataset_c_v1.json").read_text())
DATA_DIR = Path(os.environ.get("SCGEO_GSE211713_DATA_DIR", "/home/liuyuchen/data/gse211713")).resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_record(gsm: str) -> dict[str, Any]:
    columns = CONFIG["sample_columns"]
    for values in CONFIG["samples"]:
        record = dict(zip(columns, values))
        if record["geo_accession"] == gsm:
            return record
    raise KeyError(f"Unknown GSM: {gsm}")


def expected_files(record: dict[str, Any]) -> list[tuple[str, str, int]]:
    stem = f"{record['geo_accession']}_{record['sample_title']}"
    return [
        ("barcodes", f"{stem}_barcodes.tsv.gz", int(record["barcode_file_bytes"])),
        ("features", f"{stem}_genes.tsv.gz", int(CONFIG["gene_file_bytes_per_sample"])),
        ("matrix", f"{stem}_count_matrix.mtx.gz", int(record["matrix_file_bytes"])),
    ]


def download_one(url: str, destination: Path, expected_bytes: int) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and destination.stat().st_size == expected_bytes:
        return "reused_validated"
    if destination.exists():
        quarantine = destination.with_name(destination.name + f".invalid.{int(time.time())}")
        os.replace(destination, quarantine)
    partial = destination.with_name("." + destination.name + ".part")
    command = [
        "curl", "-L", "--fail", "--silent", "--show-error", "--retry", "5",
        "--retry-all-errors", "--continue-at", "-", "--output", str(partial), url,
    ]
    subprocess.run(command, check=True)
    observed = partial.stat().st_size
    if observed != expected_bytes:
        raise RuntimeError(f"Size mismatch for {destination.name}: expected {expected_bytes}, observed {observed}")
    os.replace(partial, destination)
    return "downloaded"


def download_sample(gsm: str) -> list[dict[str, Any]]:
    record = sample_record(gsm)
    directory = DATA_DIR / "processed_mex" / gsm
    rows: list[dict[str, Any]] = []
    for role, filename, expected_bytes in expected_files(record):
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6499nnn/{gsm}/suppl/{filename}"
        destination = directory / filename
        started = time.perf_counter()
        status = download_one(url, destination, expected_bytes)
        rows.append({
            "gsm": gsm, "sample_title": record["sample_title"], "role": role,
            "filename": filename, "url": url, "expected_bytes": expected_bytes,
            "observed_bytes": destination.stat().st_size, "sha256": sha256_file(destination),
            "status": status, "runtime_seconds": time.perf_counter() - started,
            "validated_utc": datetime.now(timezone.utc).isoformat(),
        })
    manifest_dir = DATA_DIR / "manifests/samples"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    path = manifest_dir / f"{gsm}_download_manifest.json"
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)
    print(json.dumps({"gsm": gsm, "files": len(rows), "status": "passed"}))
    return rows


def consolidate_manifests() -> Path:
    rows: list[dict[str, Any]] = []
    for values in CONFIG["samples"]:
        gsm = values[0]
        path = DATA_DIR / "manifests/samples" / f"{gsm}_download_manifest.json"
        rows.extend(json.loads(path.read_text()))
    output = DATA_DIR / "manifests/gse211713_processed_download_manifest.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    os.replace(temporary, output)
    return output


def run_all() -> None:
    for values in CONFIG["samples"]:
        gsm = values[0]
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--gsm", gsm], check=True)
    manifest = consolidate_manifests()
    print(json.dumps({"status": "passed", "samples": 20, "files": 60, "manifest": str(manifest)}))


def main() -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gsm")
    group.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if args.all:
        run_all()
    else:
        download_sample(args.gsm)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
