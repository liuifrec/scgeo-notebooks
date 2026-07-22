"""Phase C0 design, replication, file, and memory audit for GSE211713.

This module only transforms checksumable public-study metadata encoded in the
versioned configuration. It never reads expression values and cannot run QC,
normalization, feature selection, representations, integration, ScGeo, or a
comparator.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import psutil


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/gse211713_dataset_c_v1.json"
NOTEBOOK_PATH = ROOT / "notebooks/public_validation/gse211713/00_study_design_replication_and_file_audit.ipynb"
SCRIPT_PATH = Path(__file__).resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_output(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def validate_checkout(config: dict[str, Any]) -> tuple[Path, str]:
    branch = git_output(ROOT, "branch", "--show-current")
    if branch != config["required_git_branch"]:
        raise RuntimeError(f"Required branch {config['required_git_branch']}; observed {branch}")
    scgeo = (ROOT / config["frozen_scgeo_repository"]).resolve()
    observed = git_output(scgeo, "rev-parse", "HEAD")
    if observed != config["frozen_scgeo_commit"]:
        raise RuntimeError(f"Frozen ScGeo mismatch: {observed}")
    if git_output(scgeo, "status", "--short"):
        raise RuntimeError("Frozen ScGeo checkout is not clean")
    return scgeo, branch


def build_sample_manifest(config: dict[str, Any]) -> pd.DataFrame:
    frame = pd.DataFrame(config["samples"], columns=config["sample_columns"])
    frame["condition"] = frame["dose_gy"].map({0: "non_irradiated_control", 10: "IR_10Gy", 17: "IR_17Gy"})
    frame["time_point"] = frame["time_month"].map(lambda x: "NI_time_not_mapped" if pd.isna(x) else f"{int(x)}M")
    frame["biological_replicate_id"] = "mouse_" + frame["geo_accession"]
    frame["animal_identifier"] = frame["geo_accession"] + "_unique_mouse_proxy"
    frame["animal_identifier_basis"] = "paper states each scRNA-seq measurement/library is an independent mouse; no separate animal ID published"
    frame["technical_batch"] = "not_reported"
    frame["chemistry"] = "10x 3-prime V3/V3 NextGem; sample-level chemistry mapping unavailable"
    frame["tissue_compartment"] = "whole_lung_unsorted"
    frame["pooling_status"] = "one mouse lung per library; no biological pooling reported"
    frame["biological_replicate_eligible"] = True
    frame["eligibility_reason"] = "eligible mouse-level biological sample; contrast eligibility still depends on group sample count and design"
    frame["cell_count_basis"] = "GEO filtered barcode count before publication-level QC"
    frame["publication_retained_cells_per_sample"] = "not reported"
    frame["gene_file_bytes"] = int(config["gene_file_bytes_per_sample"])
    frame["processed_files_total_bytes"] = frame["barcode_file_bytes"] + frame["matrix_file_bytes"] + frame["gene_file_bytes"]
    frame["geo_sample_url"] = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + frame["geo_accession"]
    return frame


def build_study_design(samples: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (dose, time_point), group in samples.groupby(["dose_gy", "time_point"], sort=False):
        published = config["dataset"]["publication_retained_cells_by_dose"][str(int(dose))]
        rows.append({
            "dose_gy": int(dose),
            "time_point": time_point,
            "independent_mouse_libraries": int(len(group)),
            "deposited_filtered_cells": int(group["deposited_filtered_cells"].sum()),
            "publication_retained_cells_dose_total": int(published),
            "control_time_mapping": "unavailable" if dose == 0 else "not_applicable",
            "longitudinal_pairing": False,
            "pooling_status": "none at library preparation",
            "replication_interpretation": (
                "five independent controls, but their individual month/age mapping is not published"
                if dose == 0 else
                ("one independent mouse at this dose-time; no within-time replication" if dose == 10 else "two independent mice at this dose-time")
            ),
        })
    return pd.DataFrame(rows)


def build_candidate_replicates(samples: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "geo_accession", "sample_title", "condition", "dose_gy", "time_point",
        "biological_replicate_id", "animal_identifier", "animal_identifier_basis",
        "technical_batch", "tissue_compartment", "pooling_status",
        "deposited_filtered_cells", "biological_replicate_eligible", "eligibility_reason",
    ]
    return samples[columns].copy()


def sample_file_url(accession: str, filename: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6499nnn/{accession}/suppl/{filename}"


def build_file_inventory(samples: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in samples.itertuples(index=False):
        filenames = {
            "barcodes": (f"{row.geo_accession}_{row.sample_title}_barcodes.tsv.gz", row.barcode_file_bytes),
            "count_matrix": (f"{row.geo_accession}_{row.sample_title}_count_matrix.mtx.gz", row.matrix_file_bytes),
            "genes": (f"{row.geo_accession}_{row.sample_title}_genes.tsv.gz", row.gene_file_bytes),
        }
        for role, (filename, size) in filenames.items():
            rows.append({
                "source": "GEO", "geo_accession": row.geo_accession, "file_role": role,
                "filename": filename, "compressed_bytes": int(size),
                "url": sample_file_url(row.geo_accession, filename),
                "needed_for_phase_c1": True, "download_recommendation": "download processed file only",
            })
    rows.extend([
        {
            "source": "GEO", "geo_accession": "GSE211713", "file_role": "processed_archive",
            "filename": "GSE211713_RAW.tar", "compressed_bytes": int(config["series_archive_bytes"]),
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211713/suppl/GSE211713_RAW.tar",
            "needed_for_phase_c1": False, "download_recommendation": "optional convenience archive; individual processed files permit resumable download",
        },
        {
            "source": "SRA", "geo_accession": "PRJNA871884", "file_role": "raw_reads",
            "filename": "FASTQ/SRA runs", "compressed_bytes": pd.NA,
            "url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA871884",
            "needed_for_phase_c1": False, "download_recommendation": "do not download; processed filtered MEX matrices are sufficient",
        },
    ])
    for source in config["official_sources"]:
        if source["source_id"] in {"geo_filelist", "geo_family_soft", "supplementary_figures", "reporting_summary", "source_data", "interactive_atlas"}:
            rows.append({
                "source": source["type"], "geo_accession": "GSE211713", "file_role": source["source_id"],
                "filename": source["url"].rsplit("/", 1)[-1], "compressed_bytes": pd.NA,
                "url": source["url"], "needed_for_phase_c1": source["source_id"] == "interactive_atlas",
                "download_recommendation": "metadata/annotation provenance only; verify an official cell-level mapping before label reuse",
            })
    return pd.DataFrame(rows)


def contrast_row(identifier: str, comparison: str, scope: str, n_a: int, n_b: int, aware: bool, priority: str, reason: str) -> dict[str, Any]:
    return {
        "contrast_id": identifier, "comparison": comparison, "analysis_scope": scope,
        "independent_samples_group_a": n_a, "independent_samples_group_b": n_b,
        "replicate_aware_inference_eligible": aware, "priority": priority,
        "design_interpretation": reason, "longitudinal": False,
    }


def build_proposed_contrasts() -> pd.DataFrame:
    rows = [
        contrast_row("late17_vs_control", "17 Gy months 4-5 vs non-irradiated controls", "whole lung and adequately covered cell types", 4, 5, True, "primary", "replicate-aware cross-sectional late fibrogenic-phase comparison; controls are independent but lack individual month mapping"),
        contrast_row("early17_vs_control", "17 Gy months 1-2 vs non-irradiated controls", "whole lung and adequately covered cell types", 4, 5, True, "primary", "replicate-aware cross-sectional early response; controls are independent but lack individual month mapping"),
        contrast_row("early_vs_late17", "17 Gy months 1-2 vs months 4-5", "whole lung and adequately covered cell types", 4, 4, True, "primary", "replicate-aware cross-sectional temporal remodeling; month 3 held out as transition"),
        contrast_row("all17_vs_control", "all 17 Gy samples vs non-irradiated controls", "whole lung and adequately covered cell types", 10, 5, True, "secondary", "replicate-aware overall radiation contrast, but it averages heterogeneous post-irradiation times"),
        contrast_row("all10_vs_control", "all 10 Gy samples vs non-irradiated controls", "whole lung and adequately covered cell types", 5, 5, True, "secondary", "replicate-aware only after collapsing five distinct months; cannot estimate a within-time 10 Gy effect"),
        contrast_row("early_vs_late10", "10 Gy months 1-2 vs months 4-5", "whole lung and adequately covered cell types", 2, 2, True, "exploratory", "minimal replicate-aware cross-sectional phase comparison; low n and heterogeneous months"),
    ]
    for month in range(1, 6):
        rows.extend([
            contrast_row(f"control_vs_17gy_{month}m", f"non-irradiated controls vs 17 Gy at {month} month", "whole lung/cell type", 5, 2, True, "secondary", "mouse-level inference possible, but the same control pool is reused and control month mapping is unavailable"),
            contrast_row(f"control_vs_10gy_{month}m", f"non-irradiated controls vs 10 Gy at {month} month", "whole lung/cell type", 5, 1, False, "descriptive", "single 10 Gy mouse at this month; no radiation-group variance estimate"),
            contrast_row(f"dose_10_vs_17_{month}m", f"10 Gy vs 17 Gy at {month} month", "whole lung/cell type", 1, 2, False, "descriptive", "10 Gy arm has one mouse; dose-response inference is unsupported within month"),
        ])
    rows.extend([
        contrast_row("fibroblast_late17_vs_control", "late 17 Gy vs control in fibroblasts", "fibroblast and myofibroblast states", 4, 5, True, "primary_if_coverage", "paper reports 3,488 fibroblasts overall; require per-mouse coverage and verified labels; dissociation bias is documented"),
        contrast_row("abundance_late17_vs_control", "late 17 Gy vs control cell-type proportions", "mouse-level abundance", 4, 5, True, "secondary_cautious", "sample-level proportions permit replicate-aware modeling, but fibrotic-lung dissociation bias can dominate abundance"),
        contrast_row("three_dose_within_time", "0, 10 and 17 Gy dose response within each month", "whole lung/cell type", 5, 1, False, "not_inferential", "control month mapping is unavailable and the 10 Gy arm has one mouse per month"),
    ])
    return pd.DataFrame(rows)


def build_memory_feasibility(samples: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    cells = int(config["dataset"]["deposited_filtered_barcode_count"])
    retained = int(config["dataset"]["publication_retained_cell_count"])
    genes = int(config["dataset"]["feature_count"])
    dense_f32 = cells * genes * 4
    dense_f64 = cells * genes * 8
    retained_f32 = retained * genes * 4
    sparse_low = cells * 1000 * 8 + (cells + 1) * 4
    sparse_high = cells * 2500 * 8 + (cells + 1) * 4
    processed_sum = int(samples["processed_files_total_bytes"].sum())
    gib = 1024**3
    rows = [
        ["GEO processed archive", config["series_archive_bytes"], config["series_archive_bytes"] / gib, "download", "969 MiB archive; raw SRA reads unnecessary"],
        ["individual compressed processed files", processed_sum, processed_sum / gib, "download", "60 MEX/TSV gzip files; prefer resumable per-GSM transfer"],
        ["estimated CSR float32 low", sparse_low, sparse_low / gib, "memory", "1,000 nonzeros per cell; data+indices+indptr only"],
        ["estimated CSR float32 high", sparse_high, sparse_high / gib, "memory", "2,500 nonzeros per cell; data+indices+indptr only"],
        ["complete deposited dense float32", dense_f32, dense_f32 / gib, "danger", "prohibited full-matrix densification"],
        ["complete deposited dense float64", dense_f64, dense_f64 / gib, "danger", "prohibited and likely exceeds safe workstation memory"],
        ["publication-retained dense float32", retained_f32, retained_f32 / gib, "danger", "still prohibited; quoted only for planning"],
        ["recommended free disk", 12 * gib, 12.0, "disk_plan", "allows archive, extracted compressed files, sparse intermediates, H5AD and atomic-write headroom"],
        ["recommended peak RSS ceiling", 18 * gib, 18.0, "memory_plan", "stream one sample at a time; sparse concatenate; stop before 20 GiB warning threshold"],
    ]
    return pd.DataFrame(rows, columns=["item", "estimated_bytes", "estimated_gib", "category", "basis_and_action"])


def package_versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for package in ["pandas", "psutil", "nbformat", "nbclient"]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    return versions


def run_audit() -> dict[str, Any]:
    started = time.perf_counter()
    config = load_config()
    scgeo, branch = validate_checkout(config)
    output_dir = (ROOT / config["output_dir"]).resolve()
    audit_dir = output_dir / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    samples = build_sample_manifest(config)
    if len(samples) != 20 or samples["biological_replicate_id"].nunique() != 20:
        raise RuntimeError("Expected 20 unique independent mouse libraries")
    if int(samples["deposited_filtered_cells"].sum()) != config["dataset"]["deposited_filtered_barcode_count"]:
        raise RuntimeError("Deposited barcode total mismatch")

    study_design = build_study_design(samples, config)
    candidates = build_candidate_replicates(samples)
    inventory = build_file_inventory(samples, config)
    contrasts = build_proposed_contrasts()
    memory = build_memory_feasibility(samples, config)

    atomic_csv(study_design, audit_dir / "study_design.csv")
    atomic_csv(samples, audit_dir / "gsm_sample_manifest.csv")
    atomic_csv(candidates, audit_dir / "candidate_replicates.csv")
    atomic_csv(inventory, audit_dir / "public_file_inventory.csv")
    atomic_csv(contrasts, audit_dir / "proposed_contrasts.csv")
    atomic_csv(memory, audit_dir / "memory_feasibility.csv")

    summary = {
        "status": "passed",
        "phase": "C0_only",
        "dataset": config["dataset"],
        "replication": {
            "valid_biological_replicate_unit": "individual mouse lung represented by one GSM library",
            "valid_biological_replicates": 20,
            "control_mice": 5,
            "ten_gy_mice": 5,
            "seventeen_gy_mice": 10,
            "within_time_replicates": {"control": "time mapping unavailable", "10Gy": 1, "17Gy": 2},
            "replicate_aware_inference": "partially_supported",
            "eligible_inference": "mouse-level contrasts with at least two independent mice per group; strongest after prespecified early/late phase aggregation for 17 Gy",
            "ineligible_inference": "10 Gy within-month effects and within-month three-dose trends",
            "cells_are_replicates": False,
            "longitudinal": False,
            "pairing": "none; different mice were killed at each irradiated time point",
            "pooling": "no biological pooling reported for the 20 scRNA-seq libraries; later CellChat analyses pooled samples computationally and do not redefine replicate identity",
        },
        "confounding_and_missing_metadata": [
            "Non-irradiated controls are described as age-matched, but GEO records time as NI and provide no one-to-one control-month mapping.",
            "The 10 Gy arm has one mouse per time point, so dose/time-specific biological variance is unavailable for that arm.",
            "Sample-specific 10x V3 versus V3 NextGem chemistry and other technical-batch assignments are not provided.",
            "Dose and time must never be used as integration batch covariates.",
            "Fibrotic-lung dissociation bias is documented and limits abundance interpretation.",
        ],
        "available_annotations": {
            "geo_cell_level_mapping": False,
            "paper_main_cell_types": 18,
            "main_types": ["AT2", "AT1", "club", "ciliated", "fibroblast", "smooth muscle", "mesotheliocyte", "endothelial", "monocyte", "alveolar macrophage", "interstitial macrophage", "dendritic", "plasmacytoid dendritic", "neutrophil", "basophil", "T", "NK", "B"],
            "paper_subtypes": {
                "fibroblast": ["Col13a1-positive matrix fibroblast", "Col14a1-positive matrix fibroblast", "myofibroblast"],
                "endothelial": ["lymphatic", "artery", "vein", "gCap", "aCap"],
                "macrophage": ["AM_C1", "AM_C2", "IM_C1", "IM_C2", "IM_C3"],
                "proliferating": ["dendritic", "alveolar macrophage", "T cell"],
            },
            "fibroblast_cells_publication_total": 3488,
            "reuse_rule": "do not assign official labels to GEO barcodes unless an official cell-level mapping is recovered and validated",
        },
        "analysis_support": {
            "replicate_aware_control_vs_radiation": "yes for aggregated or 17 Gy contrasts; not for 10 Gy within a single month",
            "dose_response": "descriptive within month; no conventional replicate-aware trend because 10 Gy n=1 and control time mapping is missing",
            "time_analysis": "cross-sectional only; 17 Gy early-versus-late is replicate-aware; no longitudinal claims",
            "within_cell_type": "feasible after verified annotation and per-mouse coverage checks",
            "fibroblast_focused": "scientifically relevant and feasible in principle, but only 3,488 publication-retained fibroblasts overall and dissociation bias is substantial",
            "abundance": "mouse-level proportions are computable, but must be secondary and explicitly conditioned on dissociation/capture bias",
        },
        "recommended_representation_ensemble_for_future_phase": [
            "one sparse-normalized 50-component PCA basis on prespecified HVGs, exposed as PCA20/PCA30/PCA50 sensitivity views",
            "optional scVI without dose, condition, or time as a batch covariate",
            "diffusion map as exploratory temporal-geometry sensitivity",
            "UMAP for display only",
            "no Harmony or Scanorama unless a genuine non-condition technical batch field is recovered",
        ],
        "download_plan": "download only the 60 processed per-GSM MEX/TSV gzip files (about 0.94 GiB) or the 969 MiB archive; verify sizes, stream one library at a time, and skip SRA raw reads",
        "resource_plan": "reserve 12 GiB disk and cap future sparse ingestion near 18 GiB RSS; estimated 30-90 minutes for download, sparse ingest, merge and atomic validation; Phase C0 itself is under one minute",
        "suitability": "suitable_with_design_constraints",
        "major_risks": [
            "unbalanced replication (10 Gy n=1 versus 17 Gy n=2 per month)",
            "unmapped control time/age labels",
            "cross-sectional rather than longitudinal sampling",
            "no GEO per-cell annotation file",
            "documented dissociation bias in fibrotic lungs",
            "unknown sample-specific chemistry/technical batch",
        ],
        "forbidden_steps_executed": [],
    }
    atomic_json(summary, audit_dir / "audit_summary.json")

    runtime = time.perf_counter() - started
    provenance = {
        "status": "passed", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "repository_branch": branch, "repository_commit_before_audit": git_output(ROOT, "rev-parse", "HEAD"),
        "frozen_scgeo_repository": str(scgeo), "frozen_scgeo_commit": git_output(scgeo, "rev-parse", "HEAD"),
        "frozen_scgeo_clean": not bool(git_output(scgeo, "status", "--short")),
        "source_hashes": {
            str(CONFIG_PATH.relative_to(ROOT)): sha256_file(CONFIG_PATH),
            str(SCRIPT_PATH.relative_to(ROOT)): sha256_file(SCRIPT_PATH),
            str(NOTEBOOK_PATH.relative_to(ROOT)): sha256_file(NOTEBOOK_PATH),
        },
        "official_sources": config["official_sources"],
        "software_versions": package_versions(),
        "runtime_seconds": runtime,
        "rss_gib": psutil.Process().memory_info().rss / 1024**3,
        "numerical_expression_analysis_performed": False,
        "network_access_during_notebook": False,
    }
    atomic_json(provenance, audit_dir / "provenance.json")
    return {"summary": summary, "provenance": provenance, "output_dir": str(output_dir)}


if __name__ == "__main__":
    print(json.dumps(run_audit(), indent=2, sort_keys=True))
