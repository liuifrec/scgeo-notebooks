#!/usr/bin/env python3
"""Strict acceptance gate for the repaired GSE211713 C6-C7 v2 run."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import nbformat
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "results/public_validation/gse211713_dataset_c"
C6 = RESULT_ROOT / "c6_v2"
C7 = RESULT_ROOT / "scgeo_v2"
REP = Path("/home/liuyuchen/data/gse211713/gse211713_revision_representations_v2.h5ad")
COMPACT = Path("/home/liuyuchen/data/gse211713/gse211713_revision_qc_hvg.h5ad")
ANNOTATED = Path("/home/liuyuchen/data/gse211713/gse211713_revision_annotated.h5ad")
SCGEO = Path("/home/liuyuchen/Github/scgeo")
EXPECTED = {
    COMPACT: "340504d5c88e60a3fc70f0c044ec9f65010226e99fa976e447fa8e6698418560",
    ANNOTATED: "6f7dd2c10112960a4dca4a41b56fcc670a99e16232e0b9d86365ad5fdd66a6ca",
}
SCGEO_COMMIT = "9a0ed16cbaa57f935f9c9bc87d1643a25b51012c"
MAJOR = {"Endothelial", "Epithelial", "Fibroblast/stromal", "Lymphoid", "Myeloid", "Proliferating"}
FIBRO = {"Col13a1-like matrix fibroblast", "Col14a1-like matrix fibroblast", "Myofibroblast"}
PRIMARY_REPS = {"X_pca30", "X_pca50", "X_scvi"}
ALL_REPS = {"X_pca20", "X_pca30", "X_pca50", "X_diffmap", "X_scvi"}
PRIMARY_COUNTS = {
    "control_vs_17gy_early": (5, 4, 126),
    "control_vs_17gy_late": (5, 4, 126),
    "17gy_early_vs_late": (4, 4, 70),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args: str, repo: Path = ROOT) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def require(condition: bool, message: str, checks: list[dict]) -> None:
    checks.append({"check": message, "passed": bool(condition)})
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    checks: list[dict] = []
    for path, expected in EXPECTED.items():
        require(path.is_file() and sha256(path) == expected, f"frozen input checksum: {path}", checks)
    require(git("rev-parse", "HEAD", repo=SCGEO) == SCGEO_COMMIT, "frozen ScGeo commit", checks)
    require(not git("status", "--short", repo=SCGEO), "frozen ScGeo clean", checks)

    c6_report = json.loads((C6 / "execution/gse211713_phase_c6_v2_execution_report.json").read_text())
    c7_report = json.loads((C7 / "execution/gse211713_phase_c7_v2_execution_report.json").read_text())
    require(c6_report.get("status") == "passed", "final C6 report passed", checks)
    require(c7_report.get("status") == "passed", "final C7 report passed", checks)
    require(not c6_report.get("deviations_from_approved_protocol"), "C6 has no protocol deviations", checks)
    require(not c7_report.get("deviations_from_approved_protocol"), "C7 has no protocol deviations", checks)

    rep_hash = sha256(REP)
    assembly = json.loads((C6 / "metadata/assembly_quality_stage.json").read_text())
    require(rep_hash == assembly["representation_sha256"], "representation checksum matches C6 report", checks)
    backed = ad.read_h5ad(REP, backed="r")
    try:
        require(backed.shape == (131157, 0), "representation object is 131157 by zero genes", checks)
        require(ALL_REPS.issubset(backed.obsm.keys()), "all five quantitative representations persisted", checks)
        require("X_umap_display_only" in backed.obsm, "display-only UMAP persisted separately", checks)
        dims = {key: int(backed.obsm[key].shape[1]) for key in ALL_REPS}
        nonfinite = {key: int((~np.isfinite(np.asarray(backed.obsm[key]))).sum()) for key in ALL_REPS}
    finally:
        backed.file.close()
    require(dims == {"X_pca20": 20, "X_pca30": 30, "X_pca50": 50, "X_diffmap": 20, "X_scvi": 20}, "representation dimensions", checks)
    require(all(value == 0 for value in nonfinite.values()), "all representation coordinates finite", checks)
    config6 = json.loads((ROOT / "configs/gse211713_representations_v1.json").read_text())
    require(set(config6["representation_roles"]["primary"]) == PRIMARY_REPS, "primary consensus representations exact", checks)
    require(config6["representation_roles"]["dimensional_sensitivity"] == ["X_pca20"], "PCA20 dimensional sensitivity role", checks)
    require(config6["representation_roles"]["exploratory_sensitivity"] == ["X_diffmap"], "diffusion exploratory sensitivity role", checks)
    require(config6["scvi"]["max_epochs"] == 100 and config6["scvi"]["batch_covariate"] is None, "scVI ceiling and no batch covariate", checks)
    scvi = json.loads((C6 / "metadata/scvi_stage.json").read_text())
    history = pd.read_csv(C6 / "artifacts/scvi_training_history.csv")
    require(scvi["max_epochs"] == 100 and scvi["actual_epochs"] == len(history), "scVI actual epochs recorded", checks)
    require(scvi["accelerator"] == "gpu" and scvi["cuda_available"], "scVI GPU execution", checks)
    require(all(pd.api.types.is_numeric_dtype(history[column]) for column in history), "scVI history is tidy numeric", checks)

    notebooks = sorted((ROOT / "notebooks/public_validation/gse211713").glob("0[45]*.ipynb"))
    require(len(notebooks) == 9, "nine C6-C7 source notebooks present", checks)
    for path in notebooks:
        notebook = nbformat.read(path, as_version=4)
        clean = all(not cell.get("outputs") and cell.get("execution_count") is None for cell in notebook.cells if cell.cell_type == "code")
        require(clean, f"source notebook output-free: {path.name}", checks)
    executed = sorted((RESULT_ROOT / "executed_notebooks_v2").rglob("*.ipynb"))
    require(len(executed) == 9, "nine executed review notebooks retained", checks)
    for path in executed:
        ignored = subprocess.run(["git", "-C", str(ROOT), "check-ignore", "-q", str(path.relative_to(ROOT))]).returncode == 0
        require(ignored, f"executed notebook ignored: {path.name}", checks)

    staged = git("diff", "--cached", "--name-only").splitlines()
    forbidden_suffixes = (".h5ad", ".h5", ".pt", ".ckpt", ".tar.gz", ".mtx", ".mtx.gz")
    require(not any(path.endswith(forbidden_suffixes) or path.startswith("results/") or "/models/" in path for path in staged), "no data/model/result artifact staged", checks)

    schema = json.loads((C7 / "output_schema.json").read_text())
    full = pd.read_csv(C7 / "full_state_evidence.csv")
    major = pd.read_csv(C7 / "major_primary_state_evidence.csv")
    fibro = pd.read_csv(C7 / "fibroblast_primary_state_evidence.csv")
    consensus = pd.read_csv(C7 / "primary_consensus_state_evidence.csv")
    eligibility = pd.read_csv(C7 / "eligibility_and_coverage_used.csv")
    bootstrap = pd.read_csv(C7 / "biological_mouse_bootstrap_intervals.csv")
    centers = pd.read_csv(C7 / "individual_mouse_centers.csv")
    permutations = pd.read_csv(C7 / "exact_mouse_permutations.csv")
    permutation_values = pd.read_csv(C7 / "exact_mouse_permutation_values.csv")

    key = ["contrast", "hierarchy", "state", "representation"]
    consensus_key = ["contrast", "hierarchy", "state"]
    require(full.shape[0] == 195 and not full.duplicated(key).any(), "full evidence 195 unique canonical keys", checks)
    require(consensus.shape[0] == 27 and not consensus.duplicated(consensus_key).any(), "primary consensus 27 unique canonical keys", checks)
    require(schema["duplicate_full_evidence_keys"] == 0 and schema["duplicate_primary_consensus_keys"] == 0, "schema duplicate counts zero", checks)
    text_values = " ".join(full.astype(str).to_numpy().ravel()).lower()
    require("smoke_test" not in text_values and "smoke test" not in text_values, "no smoke test in canonical evidence", checks)
    require(not full["state"].astype(str).str.lower().eq("nan").any(), "no nan state", checks)
    require(set(major["state"]) == MAJOR and major.shape[0] == 90, "major output exact six states", checks)
    require(set(fibro["state"]) == FIBRO and fibro.shape[0] == 45, "fibroblast output exact three states", checks)
    require(set(full["representation"]) == ALL_REPS, "all five representations in evidence", checks)
    require(set(full.loc[full["representation_role"].eq("primary"), "representation"]) == PRIMARY_REPS, "only PCA30/PCA50/scVI primary", checks)

    config7 = json.loads((ROOT / "configs/gse211713_scgeo_v1.json").read_text())
    require(config7["coverage"] == {"major_min_cells_per_mouse": 30, "fibroblast_min_cells_per_mouse": 20, "minimum_eligible_mice_per_group": 3}, "frozen coverage criteria unchanged", checks)
    require(config7["estimator"]["n_boot"] == 500 and config7["estimator"]["resampling_label"] == "biological_mouse_bootstrap", "500 biological-mouse bootstraps", checks)
    require(config7["sample_key"] == "mouse_id", "sample key is mouse_id", checks)
    eligible_rows = eligibility[eligibility["eligible"].astype(bool)]
    require((eligible_rows["state_cells"] >= eligible_rows["minimum_cells_per_mouse"]).all(), "per-mouse state thresholds enforced", checks)
    require((eligibility.groupby(["contrast", "hierarchy", "state", "group"])["eligible"].sum() >= 3).all(), "at least three eligible mice per state/group", checks)

    for contrast, (n0, n1, assignments) in PRIMARY_COUNTS.items():
        rows = full[full["contrast"].eq(contrast)]
        require(set(rows["n_mice_group0"].astype(int)) == {n0} and set(rows["n_mice_group1"].astype(int)) == {n1}, f"primary mouse counts: {contrast}", checks)
        perm = permutations[permutations["contrast"].eq(contrast)]
        require(set(perm["assignment_count"].astype(int)) == {assignments}, f"exact permutation assignments: {contrast}", checks)
        require(set(perm["representation"]) == PRIMARY_REPS, f"permutations by each primary representation: {contrast}", checks)
        for group_key, frame in permutation_values[permutation_values["contrast"].eq(contrast)].groupby(["hierarchy", "state", "representation"]):
            require(frame.shape[0] == assignments and int(frame["is_observed_assignment"].sum()) == 1, f"full permutation values: {contrast}/{group_key}", checks)

    interval_columns = {"raw_ci95_low", "raw_ci95_high", "normalized_ci95_low", "normalized_ci95_high", "direction_stability", "sign_stability"}
    require(interval_columns.issubset(bootstrap.columns) and bootstrap.shape[0] == 195, "raw/normalized bootstrap intervals and stability fields", checks)
    require(set(bootstrap["resampling_status"]) == {"biological_mouse_bootstrap"}, "bootstrap label is biological_mouse_bootstrap", checks)
    center_key = ["contrast", "hierarchy", "state", "representation", "group", "mouse_id"]
    require(centers.shape[0] == 1920 and not centers.duplicated(center_key).any(), "1920 unique individual geometric-median mouse centers", checks)
    coordinate_columns = [column for column in centers if column.startswith("coordinate_")]
    require(len(coordinate_columns) == 50, "mouse center coordinate schema supports PCA50", checks)
    for row in centers.itertuples(index=False):
        values = np.asarray([getattr(row, f"coordinate_{index}") for index in range(1, int(row.dimensions) + 1)], dtype=float)
        if not np.isfinite(values).all():
            raise RuntimeError("Individual mouse center contains nonfinite coordinates")
    checks.append({"check": "all individual mouse center coordinates finite", "passed": True})

    require(set(full["inference_status"]) == {"replicate_aware_primary", "replicate_aware_secondary_time_heterogeneous"}, "canonical inference-status classes", checks)
    require(set(consensus["inference_status"]) == {"replicate_aware_primary"}, "primary consensus inference status", checks)
    require(not any("summary" in path.lower() for path in schema["canonical_sources"]), "canonical manifest excludes summaries", checks)

    diff_check = subprocess.run(["git", "-C", str(ROOT), "diff", "--check"], capture_output=True, text=True)
    require(diff_check.returncode == 0, "git diff --check", checks)
    payload = {
        "status": "passed",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "checks_passed": len(checks),
        "checks": checks,
        "input_checksums": {str(path): sha256(path) for path in EXPECTED},
        "representation_h5ad": str(REP),
        "representation_sha256": rep_hash,
        "representation_dimensions": dims,
        "representation_nonfinite_coordinates": nonfinite,
        "canonical_row_counts": schema["files"],
        "frozen_scgeo_commit": SCGEO_COMMIT,
        "frozen_scgeo_clean": True,
        "git_status_short": git("status", "--short").splitlines(),
        "git_diff_check_output": diff_check.stdout,
    }
    target = C7 / "final_acceptance_gate.json"
    temporary = target.with_name(f".{target.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(target)
    print(json.dumps({"status": "passed", "checks_passed": len(checks), "report": str(target), "representation_sha256": rep_hash}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
