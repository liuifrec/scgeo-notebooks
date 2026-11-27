"""Read-only public application helpers for ScGeo 0.2.0.

No original estimator is rerun here. The only new analysis is an explicitly
labeled, undirected state adjacency derived from a cell-neighbor graph.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
import scgeo as sg

ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.2.0"
FROZEN_COMMIT = "9a0ed16cbaa57f935f9c9bc87d1643a25b51012c"
REPS = ["X_pca20", "X_pca30", "X_pca50", "X_scvi", "X_diffmap"]
FORMATS = ("png", "svg", "pdf")
INPUTS = {
    "B": {
        "representation": ("gse249479_revision_representations.h5ad",
            "5d91ecd63f62b14d77fab04943fb5cd0f530e795ccbc0403616036568e3d9979", (34432, 0)),
        "compact": ("gse249479_revision_qc_hvg.h5ad",
            "86eed5f139e2d40ed2560d476198c11dbae423321a828d91e467013c5757b759", (34432, 27413)),
    },
    "C": {
        "representation": ("gse211713_revision_representations_v2.h5ad",
            "e8728846aff10091166b3f7d82b34afaf207142957b4daaafb8d6c496ca24cf2", (131157, 0)),
        "compact": ("gse211713_revision_qc_hvg.h5ad",
            "340504d5c88e60a3fc70f0c044ec9f65010226e99fa976e447fa8e6698418560", (131157, 23933)),
        "annotated": ("gse211713_revision_annotated.h5ad",
            "6f7dd2c10112960a4dca4a41b56fcc670a99e16232e0b9d86365ad5fdd66a6ca", (131157, 23933)),
    },
}
PRIMARY = [
    ("control_vs_17gy_early", "control", "17Gy_early"),
    ("control_vs_17gy_late", "control", "17Gy_late"),
    ("17gy_early_vs_late", "17Gy_early", "17Gy_late"),
]
OTHER = "Other retained cells"


def sha256(path):
    with Path(path).open("rb") as f:
        h = hashlib.sha256()
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(test, message):
    if not test:
        raise ValueError(message)


def env_path(name):
    value = os.environ.get(name)
    require(value, f"Set {name} to an existing local directory.")
    path = Path(value).expanduser().resolve()
    require(path.is_dir(), f"Directory absent: {path}")
    return path


class Application:
    def __init__(self, dataset):
        require(dataset in INPUTS, "Expected dataset B or C.")
        require(sg.__version__ == VERSION and importlib.metadata.version("scgeo") == VERSION,
                "Use the published ScGeo 0.2.0 environment.")
        self.dataset = dataset
        self.input_root = env_path(f"SCGEO_APPLICATION_{dataset}_ROOT")
        self.frozen = env_path("SCGEO_APPLICATION_FROZEN_RESULTS")
        base = env_path("SCGEO_APPLICATION_OUTPUT_ROOT")
        for source in (self.input_root, self.frozen):
            require(base != source and source not in base.parents,
                    "Output must be outside all source directories.")
        self.out = Path(tempfile.mkdtemp(prefix=f"dataset-{dataset}-", dir=base))
        self.records = {}
        self.plots = []
        self.graphs = []
        self.config = json.loads((ROOT / "configs" / (
            "gse249479_scgeo_v1.json" if dataset == "B" else "gse211713_scgeo_v1.json"
        )).read_text())
        self.scope = ("descriptive_only; no biological sample identity"
                      if dataset == "B" else "replicate-aware cross-sectional mouse contrasts")
        self.paths = {}
        for name, (filename, expected, shape) in INPUTS[dataset].items():
            path = self.input_root / filename
            self.record(path, expected)
            self.paths[name] = path
        print("Verified", len(self.paths), "frozen inputs; outputs:", self.out)

    def record(self, path, expected=None):
        path = Path(path).resolve()
        digest = sha256(path)
        if expected is not None:
            require(digest == expected, f"Frozen checksum mismatch: {path}")
        key = str(path)
        if key in self.records:
            require(self.records[key]["sha256"] == digest, f"Input changed: {path}")
        self.records[key] = {"sha256": digest, "expected_sha256": expected,
                             "verification": "matched_expected" if expected else "observed_only"}
        return digest

    def load(self):
        a = ad.read_h5ad(self.paths["representation"])
        require(a.shape == INPUTS[self.dataset]["representation"][2] and a.obs_names.is_unique,
                "Unexpected representation shape or duplicate cell IDs.")
        for rep in REPS:
            require(rep in a.obsm and np.isfinite(a.obsm[rep]).all(), f"Missing/nonfinite {rep}")
        annotations = None
        for name in self.paths:
            if name == "representation":
                continue
            b = ad.read_h5ad(self.paths[name], backed="r")
            try:
                require(b.shape == INPUTS[self.dataset][name][2], f"Unexpected {name} shape")
                require(b.obs_names.is_unique and a.obs_names.equals(b.obs_names),
                        f"Cell identities/order differ in {name}; stop rather than align implicitly.")
                if name == "annotated":
                    annotations = b.obs.copy()
            finally:
                b.file.close()
        return a, annotations

    def csv(self, path, expected=None):
        self.record(path, expected)
        return pd.read_csv(path)

    def b_evidence(self, a, treated):
        path = self.frozen / "gse249479_dataset_b/scgeo" / (
            f"04_{treated.lower()}_detailed_representation_state_evidence.csv")
        frame = self.csv(path)
        return adapt_b(frame, a, treated, self.config["detailed_state_order"])

    def c_evidence(self, a, contrast, group0, group1, hierarchy):
        directory = self.frozen / "gse211713_dataset_c/scgeo_v2/per_contrast" / (
            f"{contrast}__{hierarchy}")
        report_path = directory / "report.json"
        self.record(report_path)
        report = json.loads(report_path.read_text())
        require(report["status"] == "passed" and report["contrast"] == contrast
                and report["hierarchy"] == hierarchy, "Unexpected frozen contrast report")
        hashes = report["output_checksums"]
        raw = self.csv(directory / "state_evidence.csv", hashes["state_evidence.csv"])
        centers = self.csv(directory / "individual_mouse_centers.csv",
                           hashes["individual_mouse_centers.csv"])
        states = self.config["states"][hierarchy]
        threshold = self.config["coverage"]["major_min_cells_per_mouse" if hierarchy == "major"
                                           else "fibroblast_min_cells_per_mouse"]
        adapted = adapt_c(raw, a, contrast, group0, group1, hierarchy, states, threshold,
                          self.config["coverage"]["minimum_eligible_mice_per_group"])
        return adapted, adapt_centers(centers, adapted, contrast, hierarchy, group0, group1)

    def save(self, result, name, note):
        from IPython.display import Image, display
        result.metadata.update({"application_dataset": self.dataset,
                                "application_scope": self.scope,
                                "frozen_estimator_source": FROZEN_COMMIT,
                                "note": note})
        try:
            result.save(self.out / name, formats=FORMATS, dpi=165)
            self.plots.append(name)
            display(Image(filename=str(self.out / (name + ".png")), width=950))
        finally:
            plt.close(result.fig)

    def finish(self):
        for path, record in self.records.items():
            require(sha256(path) == record["sha256"], f"Source changed during execution: {path}")
        manifests = audit_exports(self.out)
        require(manifests["plot_manifests"] == len(self.plots), "Plot manifest count mismatch")
        try:
            commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                             text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            commit = None
        result = {"status": "success", "dataset": self.dataset, "scope": self.scope,
                  "package": VERSION, "package_path": sg.__file__,
                  "frozen_estimator_source": FROZEN_COMMIT, "companion_commit": commit,
                  "python": sys.version, "python_executable": sys.executable,
                  "sources": self.records, "graphs": self.graphs, "plots": self.plots,
                  "packages": {name: importlib.metadata.version(name) for name in
                               ("scgeo", "anndata", "scanpy", "numpy", "pandas", "scipy",
                                "matplotlib", "scikit-learn")}, **manifests}
        (self.out / "application_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
        print("PASS:", manifests, "outputs:", self.out)
        return result


def audit_exports(directory):
    counts = {"plot_manifests": 0, "figures": 0, "source_tables": 0}
    for path in sorted(Path(directory).rglob("*.json")):
        obj = json.loads(path.read_text())
        if not isinstance(obj, dict) or obj.get("schema_version") != "scgeo.plot_export.v1":
            continue
        require(obj["software"]["scgeo"] == VERSION, "Wrong plot software version")
        require(set(obj["figures"]) == set(FORMATS), "Missing figure format")
        counts["plot_manifests"] += 1
        for kind, key in (("figures", "figures"), ("tables", "source_tables")):
            for record in obj[kind].values():
                filename = Path(record["file"])
                require(filename.name == str(filename), "Unexpected export path")
                require(sha256(path.parent / filename) == record["sha256"], "Export hash mismatch")
                counts[key] += 1
    return counts


def check_rows(frame, states, contrast, hierarchy):
    require(not frame.duplicated(["state", "representation"]).any(), "Duplicate evidence keys")
    require(set(frame["state"]) == set(states), "Unexpected evidence state set")
    require(set(frame["contrast"]) == {contrast} and set(frame["hierarchy"]) == {hierarchy},
            "Unexpected contrast or hierarchy")
    require(set(frame["representation"]) == set(REPS) and len(frame) == len(states) * len(REPS),
            "Incomplete representation-by-state evidence")


def adapt_b(raw, a, treated, states):
    frame = raw.copy()
    check_rows(frame, states, f"PBS_vs_{treated}", "detailed")
    require(set(frame["inference_status"]) == {"descriptive_only"}, "Wrong inference scope")
    require(set(frame["status"]).issubset({"usable", "insufficient_coverage"}),
            "Unknown descriptive coverage")
    frame["n_cells0"] = frame["n_cells_PBS"].astype(int)
    frame["n_cells1"] = frame[f"n_cells_{treated}"].astype(int)
    for group, key in (("PBS", "n_cells0"), (treated, "n_cells1")):
        counts = a.obs.loc[a.obs["condition"].eq(group), "application_state"].value_counts()
        require(np.array_equal(frame[key].to_numpy(), frame["state"].map(counts).fillna(0).to_numpy()),
                "Frozen descriptive cell counts disagree with representation")
    eligible = frame["status"].eq("usable") & frame[["n_cells0", "n_cells1"]].min(axis=1).ge(20)
    frame["delta_norm"] = frame["robust_displacement_magnitude"].where(eligible)
    frame["normalized_delta_norm"] = frame["normalized_effect"].where(eligible)
    frame["effect_status"] = np.where(eligible & frame["normalized_delta_norm"].notna(),
                                     "assessed", "insufficient_coverage")
    frame["coverage_status"] = frame["status"]
    frame["descriptive_only"] = True
    frame["inference_level"] = "cell_descriptive"
    frame["group0"], frame["group1"] = "PBS", treated
    # Do not rename cell-resampling intervals to biological uncertainty.
    return frame


def adapt_c(raw, a, contrast, group0, group1, hierarchy, states, threshold, minimum_mice):
    frame = raw.copy()
    check_rows(frame, states, contrast, hierarchy)
    require(set(frame["inference_status"]) == {"replicate_aware_primary"}, "Not a primary contrast")
    require(set(frame["coverage_status"]) == {"qualified_per_mouse_before_scgeo"},
            "Unexpected coverage policy")
    require(a.obs["mouse_id"].notna().all(), "Missing mouse identity")
    require(a.obs.groupby("mouse_id", observed=True)["application_group"].nunique().le(1).all(),
            "A mouse appears in multiple cross-sectional groups")
    for group, suffix in ((group0, "0"), (group1, "1")):
        sub = a.obs.loc[a.obs["application_group"].eq(group)]
        roster = set(sub["mouse_id"].astype(str))
        numbers = []
        for row in frame.itertuples(index=False):
            counts = sub.loc[sub["application_state"].eq(row.state), "mouse_id"].astype(str).value_counts()
            eligible = {mouse: int(counts.get(mouse, 0)) for mouse in roster
                        if int(counts.get(mouse, 0)) >= threshold}
            ids = set(json.loads(getattr(row, f"eligible_mouse_ids_group{suffix}")))
            excluded = set(json.loads(getattr(row, f"excluded_mouse_ids_group{suffix}")))
            recorded = json.loads(getattr(row, f"cells_per_mouse_group{suffix}"))
            require(ids == set(eligible) and excluded == roster - ids and recorded == eligible
                    and getattr(row, f"n_mice_group{suffix}") == len(ids)
                    and len(ids) >= minimum_mice,
                    f"Frozen mouse qualification mismatch: {contrast}, {row.state}, {group}")
            numbers.append(sum(eligible.values()))
        frame[f"n_cells{suffix}"] = numbers
        frame[f"n_samples{suffix}"] = frame[f"n_mice_group{suffix}"].astype(int)
    require(np.isfinite(frame[["delta_norm", "normalized_delta_norm"]]).all().all()
            and (frame[["delta_norm", "normalized_delta_norm"]] >= 0).all().all(),
            "Invalid frozen displacement")
    frame["normalized_magnitude_ci95_low"] = frame["normalized_ci95_low"]
    frame["normalized_magnitude_ci95_high"] = frame["normalized_ci95_high"]
    frame["magnitude_ci95_low"] = frame["raw_ci95_low"]
    frame["magnitude_ci95_high"] = frame["raw_ci95_high"]
    frame["effect_status"] = "assessed"
    frame["descriptive_only"] = False
    frame["inference_level"] = "biological_sample"
    frame["bootstrap_unit"] = "sample"
    frame["group0"], frame["group1"] = group0, group1
    return frame


def adapt_centers(raw, evidence, contrast, hierarchy, group0, group1):
    frame = raw.loc[raw["representation"].eq("X_pca50")].copy()
    require(set(frame["contrast"]) == {contrast} and set(frame["hierarchy"]) == {hierarchy},
            "Wrong mouse center source")
    require(not frame.duplicated(["state", "group", "mouse_id"]).any(), "Duplicate mouse center")
    require(set(frame["state"]) == set(evidence["state"]), "Unexpected mouse center states")
    require(set(frame["group"]) == {group0, group1} and frame["dimensions"].eq(50).all()
            and np.isfinite(frame["coordinate_1"]).all(), "Invalid PCA50 mouse center rows")
    for row in evidence.loc[evidence["representation"].eq("X_pca50")].itertuples(index=False):
        for group, suffix in ((group0, "0"), (group1, "1")):
            actual = frame.loc[frame["state"].eq(row.state) & frame["group"].eq(group)]
            require(dict(zip(actual["mouse_id"].astype(str), actual["state_cells"].astype(int)))
                    == json.loads(getattr(row, f"cells_per_mouse_group{suffix}")),
                    "Mouse center membership/count mismatch")
    frame = frame.rename(columns={"group": "condition", "mouse_id": "sample",
                                  "coordinate_1": "value"})
    frame["coordinate"] = 0
    return frame[["state", "condition", "sample", "coordinate", "value",
                  "representation", "state_cells"]]


def adjacency(a, graph, state_key="application_state"):
    """H.T @ W @ H, divided by sqrt(total incident weights per state).

    W is the symmetrized nonnegative cell graph (no self edges). This is a
    descriptive normalized edge mass, NOT the Scanpy PAGA estimator.
    """
    labels = a.obs[state_key].astype(str)
    names = sorted(labels.unique())
    codes = pd.Categorical(labels, categories=names).codes
    require(len(names) > 0 and (codes >= 0).all(), "Missing graph state")
    w = sparse.csr_matrix(graph, dtype=np.float64).copy()
    require(w.shape == (a.n_obs, a.n_obs) and np.isfinite(w.data).all()
            and (w.data >= 0).all(), "Invalid cell graph")
    w = (w + w.T) * 0.5
    w.setdiag(0)
    w.eliminate_zeros()
    h = sparse.csr_matrix((np.ones(a.n_obs), (np.arange(a.n_obs), codes)),
                          shape=(a.n_obs, len(names)))
    mass = (h.T @ w @ h).toarray()
    degree = mass.sum(axis=1)
    scale = np.sqrt(np.outer(degree, degree))
    strength = np.divide(mass, scale, out=np.zeros_like(mass), where=scale > 0)
    np.fill_diagonal(strength, 0)
    strength = (strength + strength.T) * 0.5
    a.uns["application_adjacency"] = {"groups": state_key,
        "node_order": names, "connectivities": pd.DataFrame(strength, index=names, columns=names)}
    return names, pd.DataFrame(strength, index=names, columns=names)


def install_graph(app, a, graph, name, origin):
    names, table = adjacency(a, graph)
    table.to_csv(app.out / f"{name}_normalized_adjacency.csv")
    app.graphs.append({"name": name, "source": origin, "nodes": names,
                       "scope": "pooled supplied graph; descriptive normalized edge mass, not PAGA",
                       "formula": "H.T @ sym(W) @ H / sqrt(row_sum outer product), diagonal zero"})
    return names


def report_for(frame):
    report = frame.loc[frame["representation"].eq("X_pca50")].copy()
    require(not report["state"].duplicated().any(), "Duplicate report state")
    return report


def plot_evidence(app, a, frame, group0, group1, condition, basis, name, names, sample=None):
    report = report_for(frame)
    report.to_csv(app.out / f"{name}_adapted_pca50_report.csv", index=False)
    frame.to_csv(app.out / f"{name}_adapted_all_representations.csv", index=False)
    scope = "Frozen descriptive displacement; no biological replication" if sample is None else (
        "Frozen qualified-mouse displacement; cross-sectional")
    title = f"{group0} vs {group1}"
    app.save(sg.pl.state_evidence_panel(
        report, effect_metric="normalized_delta_norm", show_ci=sample is not None,
        show_sample_counts=sample is not None, show_inference=True,
        max_states=None, title=f"{title} | PCA50", return_result=True, show=False),
        f"{name}_evidence", scope)
    app.save(sg.pl.robustness_matrix(
        frame, row_key="state", col_key="representation", value_key="normalized_delta_norm",
        col_order=REPS, summary=None, cbar_label="Normalized displacement",
        title=f"{title} | Representation sensitivity", return_result=True, show=False),
        f"{name}_representation_matrix", "Stored scalar magnitudes; representations not independent")
    common = dict(node_key="application_state", condition_key=condition,
                  group0=group0, group1=group1, basis=basis, paga_key="application_adjacency",
                  graph_node_order=names, node_position="pooled", connectivity_threshold=.01,
                  max_points=6000, random_state=1729, return_result=True, show=False)
    pie = dict(pie_key=condition, pie_mode="sample_abundance" if sample else "cell_fraction")
    if sample:
        pie["pie_sample_key"] = sample
    app.save(sg.pl.paga_shift_map(
        a, **common, **pie, show_shift=False, show_endpoints=False, show_velocity=False,
        node_scale_by_n=False, title=f"{title} | Neighbor adjacency and composition"),
        f"{name}_composition_pies",
        "PAGA-style renderer of descriptive pooled neighbor adjacency, not PAGA inference; "
        + ("equal-mouse state frequencies within condition" if sample else "captured-cell composition"))
    app.save(sg.pl.recovery_compass(
        a, **common, mode="stored", shift_df=report, analysis_rep="X_pca50",
        effect_metric="normalized_delta_norm", alignment_df=None, sample_key=sample,
        min_cells=20, layout="summary", fill_color_mode="constant", ring_mode="none",
        node_size_mode="constant", title=f"{title} | State response compass"),
        f"{name}_compass", scope + "; alignment/OOD absent; connectors show display means, "
        "not quantitative displacements; graph is descriptive adjacency")
    return report
