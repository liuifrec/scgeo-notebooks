"""Small, fully synthetic deployment validation; no downloads or biological inputs."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import runpy

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scgeo as sg


def synthetic_design():
    """Four independent synthetic samples; unequal cell coverage is intentional."""
    coverage = {"shared": [12, 8, 10, 6], "sparse": [8, 7, 6, 0],
                "control_only": [6, 5, 0, 0], "tiny": [1, 1, 1, 1]}
    rows, coordinates = [], []
    for j, (state, counts) in enumerate(coverage.items()):
        for sample, n in enumerate(counts):
            condition = "control" if sample < 2 else "treated"
            offset = np.array([int(sample >= 2) + (-.15 if sample % 2 == 0 else .15), j * 2.0])
            t = np.linspace(-.1, .1, n) if n > 1 else np.zeros(n)
            coordinates.extend(offset + np.column_stack([t, t ** 2 - np.mean(t ** 2)]) if n else [])
            rows.extend({"state": state, "condition": condition, "sample": f"sample_{sample}"} for _ in range(n))
    obs = pd.DataFrame(rows, index=[f"synthetic_{i:03d}" for i in range(len(rows))])
    data = ad.AnnData(X=None, obs=obs)
    data.obsm["X_primary"] = np.array(coordinates)
    data.obsm["X_rotated"] = data.obsm["X_primary"] @ np.array([[0., -1.], [1., 0.]])
    data.obsm["X_display"] = data.obsm["X_primary"].copy()
    data.obsm["velocity_display"] = np.tile([.1, 0.], (data.n_obs, 1))
    data.uns["synthetic_projection"] = {"basis": "X_display", "velocity_key": "velocity_display",
                                        "source": "hand-supplied synthetic vectors; no fitted velocity"}
    data.uns["paga"] = {"groups": "state", "node_order": list(coverage),
                        "connectivities": np.ones((4, 4)) - np.eye(4),
                        "source": "hand-supplied synthetic undirected graph; no fitted PAGA"}
    config = {"purpose": "sample_inference", "design": "independent", "state_key": "state",
              "condition_key": "condition", "group0": "control", "group1": "treated",
              "sample_key": "sample", "reps": ["X_primary", "X_rotated"],
              "methods": ["robust_shift", "representation_stability", "local_geometry_stability", "plotting"],
              "min_cells": 4, "min_samples_per_condition": 2, "bootstrap_unit": "sample", "k_values": [3],
              "embedding_key": "X_display",
              "representation_provenance": {"X_primary": "deterministic synthetic coordinates",
                                             "X_rotated": "orthogonal rotation of X_primary; dependent sensitivity check"},
              "projected_velocity": {"key": "velocity_display", "basis": "X_display",
                                     "provenance_uns_key": "synthetic_projection"},
              "paga": {"key": "paga", "node_order": list(coverage)}}
    return data, config


def analyze(data):
    """Explicit synthetic analysis, separate from all read-only validation."""
    common = dict(condition_key="condition", group0="control", group1="treated",
                  sample_key="sample", center="mean", n_boot=12, seed=7)
    sg.tl.robust_shift(data, rep="X_primary", by="state", bootstrap_unit="sample", **common)
    sg.tl.representation_stability(data, reps=["X_primary", "X_rotated"], node_key="state",
                                   min_cells=4, **common)
    sg.tl.local_geometry_stability(data, reps=["X_primary", "X_rotated"], node_key="state",
                                  sample_key="sample", k_values=[3], n_boot=12, seed=7)
    return sg.get.state_report(data, node_key="state", strict=True)


def exact_assignments(data):
    """Enumerate the six 2-vs-2 allocations; this is not a new ScGeo test method."""
    selected = data.obs.state.eq("shared").to_numpy()
    centers = pd.DataFrame(data.obsm["X_primary"][selected], columns=["x", "y"])
    centers["sample"] = data.obs.loc[selected, "sample"].to_numpy()
    centers = centers.groupby("sample", sort=True)[["x", "y"]].mean()
    rows = []
    for treated in itertools.combinations(centers.index, 2):
        control = [s for s in centers.index if s not in treated]
        delta = centers.loc[list(treated)].mean() - centers.loc[control].mean()
        rows.append({"treated_samples": list(treated), "delta_norm": float(np.linalg.norm(delta))})
    return {"n_assignments": len(rows), "assignments": rows,
            "interpretation": "Synthetic enumeration only. Six assignments do not support a fine-grained tail probability; "
                              "complementary allocations have the same unsigned statistic. No significance claim."}


def run(output_dir, *, expected_commit, source_root, formats=("png", "svg", "pdf")):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    verify = runpy.run_path(str(Path(source_root) / "scripts/deployment_preflight.py"))["environment_report"]
    environment = verify(expected_commit=expected_commit, source_root=source_root)
    if not environment["ok"]:
        (output_dir / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
        raise RuntimeError("Environment preflight failed; inspect environment.json before running analysis.")
    data, config = synthetic_design()
    contract = sg.get.input_contract(data, config=config)
    assert contract["ok"] and contract["eligible_states"] == ["shared"]
    digest = hashlib.sha256(data.obs.to_csv().encode())
    for key in sorted(data.obsm):
        digest.update(key.encode() + data.obsm[key].tobytes())
    state_report = analyze(data)
    expected_params = {
        "robust_shift": {"rep": "X_primary", "by": "state", "sample_key": "sample", "resolved_bootstrap_unit": "sample",
                         "condition_key": "condition", "group0": "control", "group1": "treated"},
        "representation_stability": {"reps": ["X_primary", "X_rotated"], "sample_key": "sample", "node_key": "state",
                                     "condition_key": "condition", "group0": "control", "group1": "treated", "min_cells": 4},
        "local_geometry_stability": {"reps": ["X_primary", "X_rotated"], "sample_key": "sample", "node_key": "state", "k_values": [3]}}
    config["evidence"] = {key: {"method": key, "params": params} for key, params in expected_params.items()}
    post = sg.get.input_contract(data, config=config)
    assert post["ok"]
    shared_delta = data.uns["scgeo"]["robust_shift"]["by"]["shared"]["delta"]
    np.testing.assert_allclose(shared_delta, [1., 0.], atol=1e-6)
    state_report.to_csv(output_dir / "state_report.csv", index=False)
    exports = []
    with sg.pl.style_context("paper"):
        for name, result in [
            ("state_evidence", sg.pl.state_evidence_panel(state_report, show_sample_counts=True,
                 show_inference=True, title="Synthetic deployment check · coverage is not significance",
                 return_result=True, show=False)),
            ("sample_effects", sg.pl.replicate_effect(data, coordinate=0, scope="state",
                 title="Synthetic sample centers · sparse states remain descriptive",
                 return_result=True, show=False)),
        ]:
            result.metadata.update(dataset="fully_synthetic", source_commit=expected_commit,
                                   input_sha256=digest.hexdigest(), input_contract=contract)
            files = result.save(output_dir / name, formats=formats, dpi=150)
            manifest = json.loads(Path(files["metadata"]).read_text())
            for kind in ("figures", "tables"):
                for item in manifest[kind].values():
                    assert hashlib.sha256((output_dir / item["file"]).read_bytes()).hexdigest() == item["sha256"]
                    if kind == "tables":
                        assert len(pd.read_csv(output_dir / item["file"])) == item["rows"]
            exports.append(manifest)
            plt.close(result.fig)
    report = {"schema_version": 1, "dataset": "fully_synthetic", "environment": environment,
              "input_sha256": digest.hexdigest(), "input_contract": contract, "evidence_contract": post,
              "shared_delta": np.asarray(shared_delta).tolist(), "permutations": exact_assignments(data),
              "state_report": json.loads(state_report.to_json(orient="records")), "exports": exports,
              "acceptance": {"environment_verified": True, "contract_passed": True,
                             "coverage_exclusions_preserved": contract["excluded_states"],
                             "source_hashes_verified": True, "biological_validation": False}}
    (output_dir / "validation.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    # A portable, concise record suitable for review. Full paths, the complete
    # package inventory and complete table metadata remain in validation.json.
    core_packages = {"scgeo", "anndata", "numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "plotly"}
    summary = {"schema_version": 1, "dataset": report["dataset"], "acceptance": report["acceptance"],
               "source_commit": expected_commit, "package_sha256": environment["source"]["package_sha256"],
               "python": environment["python"]["version"],
               "packages": {p["name"]: p["version"] for p in environment["packages"] if p["name"] in core_packages},
               "input_sha256": report["input_sha256"], "shared_delta": report["shared_delta"],
               "state_coverage": [{k: row[k] for k in ("state", "n_cells", "n_samples", "status")}
                                  for row in contract["state_coverage"]],
               "permutations": report["permutations"],
               "exports": [{"plot": name, "figures": manifest["figures"],
                            "tables": {key: {k: item[k] for k in ("file", "sha256", "rows")}
                                       for key, item in manifest["tables"].items()}}
                           for name, manifest in zip(("state_evidence", "sample_effects"), exports)],
               "reproduce": "python examples/deployment_smoke.py --source-root CHECKOUT "
                            f"--expected-commit {expected_commit} --output-dir NEW_OUTPUT_DIRECTORY"}
    (output_dir / "validation-summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--expected-commit", required=True)
    args = parser.parse_args()
    result = run(args.output_dir, source_root=args.source_root, expected_commit=args.expected_commit)
    print(json.dumps(result["acceptance"], indent=2))
