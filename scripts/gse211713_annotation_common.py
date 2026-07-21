#!/usr/bin/env python3
"""Condition-blind coarse annotation and mouse/state coverage for GSE211713."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import igraph as ig
import leidenalg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/public_validation/gse211713_dataset_c"
DATA_DIR = Path(os.environ.get("SCGEO_GSE211713_DATA_DIR", "/home/liuyuchen/data/gse211713")).resolve()
COMPACT = DATA_DIR / "gse211713_revision_qc_hvg.h5ad"
ANNOTATED = DATA_DIR / "gse211713_revision_annotated.h5ad"
MARKER_PATH = ROOT / "configs/gse211713_annotation_markers_v1.json"
MARKERS = json.loads(MARKER_PATH.read_text())
SEED = int(MARKERS["random_seed"])
WARNING_GIB = 20.0
STOP_GIB = 24.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(records: list[dict[str, Any]], stage: str, started: float) -> None:
    rss = psutil.Process().memory_info().rss / 1024**3
    records.append({"stage": stage, "rss_gib": rss, "elapsed_seconds": time.perf_counter() - started, "timestamp_utc": datetime.now(timezone.utc).isoformat()})
    if rss >= STOP_GIB:
        raise MemoryError(f"Annotation hard stop at {rss:.2f} GiB during {stage}")
    if rss >= WARNING_GIB:
        print(f"WARNING: annotation RSS {rss:.2f} GiB during {stage}")


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    frame.to_csv(temporary, index=False, compression=compression)
    os.replace(temporary, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def robust_z(values: np.ndarray) -> np.ndarray:
    median = np.median(values); mad = np.median(np.abs(values - median)); scale = 1.4826 * mad
    if scale <= 0:
        scale = values.std()
    return ((values - median) / scale if scale > 0 else np.zeros_like(values)).astype(np.float32)


def module_scores(adata: ad.AnnData, modules: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    symbol_to_indices: dict[str, list[int]] = {}
    for index, symbol in enumerate(symbols):
        symbol_to_indices.setdefault(symbol, []).append(index)
    total = adata.obs["total_counts"].to_numpy(np.float64)
    scale = np.divide(1e4, total, out=np.zeros_like(total), where=total > 0)
    scores: dict[str, np.ndarray] = {}; fractions: dict[str, np.ndarray] = {}; definitions = []
    for name, spec in modules.items():
        available = [gene for gene in spec["genes"] if gene in symbol_to_indices]
        indices = [symbol_to_indices[gene][0] for gene in available]
        if not indices:
            values = np.zeros(adata.n_obs, dtype=np.float32); fraction = values.copy()
        else:
            sub = adata.X[:, indices].tocsr(copy=True)
            sub = sparse.diags(scale.astype(np.float32)) @ sub
            np.log1p(sub.data, out=sub.data)
            values = np.asarray(sub.sum(axis=1)).ravel().astype(np.float32) / len(spec["genes"])
            fraction = adata.X[:, indices].getnnz(axis=1).astype(np.float32) / len(spec["genes"])
        scores[name] = robust_z(values); fractions[name] = fraction
        definitions.append({
            "module": name, "parent": spec.get("parent", "major"), "source": spec["source"],
            "requested_genes": ";".join(spec["genes"]), "available_genes": ";".join(available),
            "n_requested": len(spec["genes"]), "n_available": len(available),
            "score": "mean log1p counts per 10,000 across requested markers, then robust-z standardized; absent markers contribute zero",
        })
    return pd.DataFrame(scores, index=adata.obs_names), pd.DataFrame(fractions, index=adata.obs_names), pd.DataFrame(definitions)


def annotation_pca(adata: ad.AnnData, memory: list[dict[str, Any]], started: float) -> tuple[np.ndarray, np.ndarray, sparse.csr_matrix]:
    hvg = adata.var["highly_variable"].to_numpy(bool)
    matrix = adata.X[:, hvg].tocsr(copy=True)
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0).astype(np.float32)
    matrix = sparse.diags(scale) @ matrix
    np.log1p(matrix.data, out=matrix.data)
    record(memory, "hvg_sparse_log_matrix", started)
    model = PCA(n_components=int(MARKERS["annotation_pca_components"]), svd_solver="arpack", random_state=SEED)
    coordinates = model.fit_transform(matrix).astype(np.float32)
    record(memory, "annotation_pca", started)
    return coordinates, model.explained_variance_ratio_.astype(np.float32), matrix


def cluster_graph(coordinates: np.ndarray, memory: list[dict[str, Any]], started: float) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    os.environ.setdefault("NUMBA_CACHE_DIR", str(RESULTS / "_numba_cache"))
    from pynndescent import NNDescent
    index = NNDescent(coordinates, n_neighbors=int(MARKERS["annotation_neighbors"]) + 1, metric="euclidean", random_state=SEED, n_jobs=1, low_memory=True)
    neighbors, distances = index.neighbor_graph
    neighbors = neighbors[:, 1:]; distances = distances[:, 1:]
    sources = np.repeat(np.arange(coordinates.shape[0], dtype=np.int32), neighbors.shape[1])
    targets = neighbors.reshape(-1).astype(np.int32)
    graph = ig.Graph(n=coordinates.shape[0], edges=list(zip(sources.tolist(), targets.tolist())), directed=False)
    graph.simplify(multiple=True, loops=True)
    record(memory, "annotation_neighbor_graph", started)
    clusters: dict[str, np.ndarray] = {}
    for resolution in MARKERS["leiden_resolutions"]:
        partition = leidenalg.find_partition(
            graph, leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=float(resolution), seed=SEED, n_iterations=-1,
        )
        clusters[f"leiden_{resolution}"] = np.asarray(partition.membership, dtype=np.int32)
        record(memory, f"leiden_{resolution}", started)
    return clusters, neighbors, distances


def assign_annotations(adata: ad.AnnData, clusters: dict[str, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    major_scores, major_fraction, major_defs = module_scores(adata, MARKERS["major_modules"])
    refined_scores, refined_fraction, refined_defs = module_scores(adata, MARKERS["refined_modules"])
    primary_key = f"leiden_{MARKERS['primary_leiden_resolution']}"
    primary = clusters[primary_key].astype(str)
    score_table = major_scores.copy(); score_table["cluster"] = primary
    cluster_means = score_table.groupby("cluster", observed=False).mean()
    lineage_names = [name for name in MARKERS["major_modules"] if name != "proliferating"]
    rules = MARKERS["assignment_rules"]
    cluster_labels: dict[str, str] = {}; cluster_margins: dict[str, float] = {}
    for cluster, row in cluster_means.iterrows():
        ranked = row[lineage_names].sort_values(ascending=False)
        margin = float(ranked.iloc[0] - ranked.iloc[1])
        if ranked.iloc[0] < rules["cluster_major_minimum_z"] or margin < rules["cluster_major_minimum_margin"]:
            label = rules["ambiguous_label"]
        else:
            label = str(ranked.index[0])
        cluster_labels[str(cluster)] = label; cluster_margins[str(cluster)] = margin
    obs = pd.DataFrame(index=adata.obs_names)
    obs["annotation_cluster"] = primary
    obs["major_annotation"] = pd.Series(primary, index=obs.index).map(cluster_labels).to_numpy()
    obs["major_cluster_margin"] = pd.Series(primary, index=obs.index).map(cluster_margins).to_numpy(float)
    proliferating = (major_scores["proliferating"] >= rules["proliferating_override_z"]) & (major_fraction["proliferating"] >= rules["minimum_detected_marker_fraction_for_moderate"])
    obs.loc[proliferating, "major_annotation"] = "proliferating"
    obs["annotation_state"] = obs["major_annotation"]
    obs["refined_marker_margin"] = np.nan
    for parent in ["epithelial", "endothelial", "myeloid", "lymphoid", "fibroblast/stromal"]:
        modules = [name for name, spec in MARKERS["refined_modules"].items() if spec["parent"] == parent]
        mask = obs["major_annotation"].eq(parent).to_numpy()
        if not mask.any(): continue
        values = refined_scores.loc[mask, modules].to_numpy()
        order = np.argsort(values, axis=1)
        top_index = order[:, -1]; second_index = order[:, -2] if len(modules) > 1 else order[:, -1]
        top = values[np.arange(values.shape[0]), top_index]; second = values[np.arange(values.shape[0]), second_index]
        margin = top - second if len(modules) > 1 else top
        labels = np.asarray(modules, dtype=object)[top_index]
        insufficient = (top < rules["refined_minimum_z"]) | (margin < rules["refined_minimum_margin"])
        labels[insufficient] = rules["ambiguous_fibroblast_label"] if parent == "fibroblast/stromal" else parent + "_ambiguous"
        obs.loc[mask, "annotation_state"] = labels
        obs.loc[mask, "refined_marker_margin"] = margin
    top_fraction = np.zeros(adata.n_obs, dtype=np.float32)
    for name in MARKERS["major_modules"]:
        mask = obs["major_annotation"].eq(name).to_numpy()
        if mask.any(): top_fraction[mask] = major_fraction.loc[mask, name]
    obs["detected_major_marker_fraction"] = top_fraction
    confidence = np.full(adata.n_obs, "low", dtype=object)
    moderate = (obs["major_cluster_margin"].to_numpy() >= rules["moderate_confidence_margin"]) & (top_fraction >= rules["minimum_detected_marker_fraction_for_moderate"])
    high = (obs["major_cluster_margin"].to_numpy() >= rules["high_confidence_margin"]) & (top_fraction >= rules["minimum_detected_marker_fraction_for_moderate"])
    confidence[moderate] = "moderate"; confidence[high] = "high"
    confidence[obs["major_annotation"].eq(rules["ambiguous_label"]).to_numpy()] = "ambiguous"
    obs["annotation_confidence"] = confidence
    for name in major_scores:
        safe_name = name.replace("/", "_")
        obs[f"major_score_{safe_name}"] = major_scores[name].to_numpy(np.float32)
    for name in refined_scores: obs[f"refined_score_{name}"] = refined_scores[name].to_numpy(np.float32)
    cluster_output = cluster_means.reset_index()
    cluster_output["major_annotation"] = cluster_output["cluster"].map(cluster_labels)
    cluster_output["major_margin"] = cluster_output["cluster"].map(cluster_margins)
    definitions = pd.concat([major_defs.assign(level="major"), refined_defs.assign(level="refined")], ignore_index=True)
    score_summary = pd.concat([major_scores.add_prefix("major_score_"), refined_scores.add_prefix("refined_score_")], axis=1)
    score_summary["annotation_state"] = obs["annotation_state"]
    state_scores = score_summary.groupby("annotation_state", observed=False).mean().reset_index()
    return obs, cluster_output, definitions, state_scores


def build_coverage(obs: pd.DataFrame) -> dict[str, pd.DataFrame]:
    count_rows = []
    for hierarchy, column in [("major", "major_annotation"), ("refined", "annotation_state")]:
        grouped = obs.groupby(["mouse_id", "irradiation_group", "month_post_irradiation", column], dropna=False, observed=False).size().reset_index(name="retained_cells")
        grouped = grouped.rename(columns={column: "state"}); grouped.insert(0, "hierarchy", hierarchy)
        count_rows.append(grouped)
    counts = pd.concat(count_rows, ignore_index=True)
    contrasts = {
        "control_vs_17Gy_early": ("control", "17Gy_early", "primary"),
        "control_vs_17Gy_late": ("control", "17Gy_late", "primary"),
        "17Gy_early_vs_late": ("17Gy_early", "17Gy_late", "primary"),
        "control_vs_all_17Gy": ("control", "all_17Gy", "secondary_time_heterogeneous"),
        "control_vs_all_10Gy": ("control", "all_10Gy", "secondary_time_heterogeneous"),
        "control_vs_17Gy_month3": ("control", "17Gy_month3", "descriptive_only"),
    }
    for month in range(1, 6):
        contrasts[f"control_vs_10Gy_{month}M"] = ("control", f"10Gy_{month}M", "descriptive_only")
        contrasts[f"10Gy_vs_17Gy_{month}M"] = (f"10Gy_{month}M", f"17Gy_{month}M", "descriptive_only")
    memberships: dict[str, pd.Series] = {
        "control": obs["irradiation_group"].eq("control"),
        "17Gy_early": obs["irradiation_group"].eq("17Gy") & obs["month_post_irradiation"].isin([1,2]),
        "17Gy_late": obs["irradiation_group"].eq("17Gy") & obs["month_post_irradiation"].isin([4,5]),
        "all_17Gy": obs["irradiation_group"].eq("17Gy"), "all_10Gy": obs["irradiation_group"].eq("10Gy"),
        "17Gy_month3": obs["irradiation_group"].eq("17Gy") & obs["month_post_irradiation"].eq(3),
    }
    for month in range(1,6):
        memberships[f"10Gy_{month}M"] = obs["irradiation_group"].eq("10Gy") & obs["month_post_irradiation"].eq(month)
        memberships[f"17Gy_{month}M"] = obs["irradiation_group"].eq("17Gy") & obs["month_post_irradiation"].eq(month)
    eligibility_rows = []
    for hierarchy, column in [("major", "major_annotation"), ("refined", "annotation_state")]:
        states = sorted(obs[column].unique())
        for contrast, (group_a, group_b, scope) in contrasts.items():
            for state in states:
                row = {"contrast": contrast, "scope": scope, "hierarchy": hierarchy, "state": state, "minimum_cells_per_mouse": 30, "minimum_eligible_mice_per_group": 3}
                eligible_counts = []
                for label, group in [("a", group_a), ("b", group_b)]:
                    mask = memberships[group] & obs[column].eq(state)
                    per_mouse = obs.loc[mask].groupby("mouse_id", observed=False).size()
                    eligible = int((per_mouse >= 30).sum())
                    row[f"group_{label}"] = group; row[f"total_mice_group_{label}"] = int(obs.loc[memberships[group], "mouse_id"].nunique())
                    row[f"eligible_mice_group_{label}"] = eligible; eligible_counts.append(eligible)
                row["annotation_qualified"] = "ambiguous" not in state.lower() and "low-quality" not in state.lower()
                row["eligible_standard_coverage"] = min(eligible_counts) >= 3
                row["eligible_standard"] = row["eligible_standard_coverage"] and row["annotation_qualified"]
                is_fibro = hierarchy == "refined" and ("fibroblast" in state.lower())
                if is_fibro:
                    exploratory = []
                    for group in [group_a, group_b]:
                        per_mouse = obs.loc[memberships[group] & obs[column].eq(state)].groupby("mouse_id", observed=False).size()
                        exploratory.append(int((per_mouse >= 20).sum()))
                    row["fibroblast_eligible_mice_group_a_at20"] = exploratory[0]; row["fibroblast_eligible_mice_group_b_at20"] = exploratory[1]
                    row["eligible_exploratory_fibroblast20"] = (
                        min(exploratory) >= 3
                        and row["annotation_qualified"]
                        and not row["eligible_standard"]
                    )
                else:
                    row["fibroblast_eligible_mice_group_a_at20"] = np.nan; row["fibroblast_eligible_mice_group_b_at20"] = np.nan; row["eligible_exploratory_fibroblast20"] = False
                if not row["annotation_qualified"]: reason = "excluded despite numeric coverage: ambiguous marker-inferred state"
                elif row["eligible_standard"]: reason = "eligible: at least 30 cells in at least 3 mice per group"
                elif row["eligible_exploratory_fibroblast20"]: reason = "exploratory fibroblast only: fails 30-cell standard but passes prespecified 20-cell threshold"
                else: reason = "insufficient mouse/state coverage under prespecified thresholds"
                row["status_reason"] = reason; eligibility_rows.append(row)
    eligibility = pd.DataFrame(eligibility_rows)
    excluded = eligibility[~eligibility["eligible_standard"]].copy()
    fibro = counts[(counts.hierarchy == "refined") & counts.state.str.contains("fibroblast", case=False, regex=False)].copy()
    fibro["at_least_20_cells"] = fibro.retained_cells >= 20; fibro["at_least_30_cells"] = fibro.retained_cells >= 30
    return {"counts": counts, "eligibility": eligibility, "excluded": excluded, "fibro": fibro}


def save_figures(obs: pd.DataFrame, cluster_table: pd.DataFrame, coverage: dict[str, pd.DataFrame], umap: np.ndarray) -> None:
    directory = RESULTS / "annotation/figures"; directory.mkdir(parents=True, exist_ok=True)
    major_order = sorted(obs.major_annotation.unique()); palette = dict(zip(major_order, plt.cm.tab10(np.linspace(0,1,len(major_order)))))
    fig, axes = plt.subplots(1,2,figsize=(14,6),constrained_layout=True)
    for label in major_order:
        mask = obs.major_annotation.eq(label).to_numpy(); axes[0].scatter(umap[mask,0],umap[mask,1],s=.3,alpha=.35,color=palette[label],rasterized=True,label=label)
    axes[0].set(title="Conservative major annotation",xlabel="display UMAP 1",ylabel="display UMAP 2"); axes[0].legend(frameon=False,markerscale=8,fontsize=7)
    fibro = obs.major_annotation.eq("fibroblast/stromal").to_numpy(); fibro_labels = sorted(obs.loc[fibro,"annotation_state"].unique()); fibro_palette=dict(zip(fibro_labels,plt.cm.Set2(np.linspace(0,1,len(fibro_labels)))))
    axes[1].scatter(umap[~fibro,0],umap[~fibro,1],s=.2,color="#DDDDDD",alpha=.12,rasterized=True)
    for label in fibro_labels:
        mask=obs.annotation_state.eq(label).to_numpy(); axes[1].scatter(umap[mask,0],umap[mask,1],s=1,alpha=.55,color=fibro_palette[label],rasterized=True,label=label)
    axes[1].set(title="Fibroblast-state marker inference",xlabel="display UMAP 1",ylabel="display UMAP 2"); axes[1].legend(frameon=False,markerscale=6,fontsize=7)
    fig.suptitle("GSE211713 condition-blind annotation; UMAP is display only")
    for ext in ("png","svg"): fig.savefig(directory/f"gse211713_annotation_umap.{ext}",dpi=250,bbox_inches="tight",metadata={"Creator":"GSE211713 condition-blind annotation","Date":None})
    plt.close(fig)
    primary = coverage["eligibility"][coverage["eligibility"].scope.eq("primary") & coverage["eligibility"].hierarchy.eq("refined")]
    matrix = primary.pivot(index="state",columns="contrast",values="eligible_standard").astype(float)
    fig,ax=plt.subplots(figsize=(8,max(5,.28*len(matrix))),constrained_layout=True); image=ax.imshow(matrix,aspect="auto",cmap="RdYlGn",vmin=0,vmax=1)
    ax.set_xticks(range(len(matrix.columns)),matrix.columns,rotation=30,ha="right"); ax.set_yticks(range(len(matrix.index)),matrix.index,fontsize=7); ax.set_title("Primary-contrast mouse/state coverage (30 cells in at least 3 mice per group)")
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]): ax.text(x,y,"eligible" if matrix.iloc[y,x] else "no",ha="center",va="center",fontsize=6)
    for ext in ("png","svg"): fig.savefig(directory/f"gse211713_primary_coverage.{ext}",dpi=250,bbox_inches="tight",metadata={"Creator":"GSE211713 mouse-level coverage","Date":None})
    plt.close(fig)
    (RESULTS/"annotation/gse211713_annotation_umap_alt.txt").write_text("Two display-only UMAP panels from a condition-blind annotation PCA: conservative major lung compartments and marker-inferred fibroblast states. Grey cells in the fibroblast panel are other compartments. Publication labels were not transferred because no barcode mapping was available.\n")
    (RESULTS/"coverage/gse211713_primary_coverage_alt.txt").write_text("Heatmap marking whether each refined state meets at least 30 retained cells in at least three independent mouse libraries per group for each of the three primary contrasts.\n")


def run() -> dict[str, Any]:
    started=time.perf_counter(); memory=[]
    adata=ad.read_h5ad(COMPACT); record(memory,"compact_loaded",started)
    coordinates,variance,hvg_log=annotation_pca(adata,memory,started)
    clusters,neighbors,distances=cluster_graph(coordinates,memory,started)
    obs,cluster_table,definitions,state_scores=assign_annotations(adata,clusters); record(memory,"marker_annotation",started)
    for key, values in clusters.items(): obs[key]=values.astype(str)
    adata.obs=adata.obs.join(obs)
    adata.obsm["X_pca_annotation_only"]=coordinates
    adata.uns["annotation_pca_explained_variance_ratio"]=variance
    os.environ.setdefault("NUMBA_CACHE_DIR",str(RESULTS/"_numba_cache"))
    import umap
    display=umap.UMAP(n_neighbors=15,min_dist=.3,n_components=2,metric="euclidean",random_state=SEED,n_jobs=1,init="random",precomputed_knn=(neighbors,distances,None)).fit_transform(coordinates).astype(np.float32)
    adata.obsm["X_umap_annotation_display_only"]=display; record(memory,"display_umap",started)
    coverage=build_coverage(adata.obs)
    atomic_csv(coverage["counts"],RESULTS/"coverage/mouse_state_cell_counts.csv")
    atomic_csv(coverage["eligibility"],RESULTS/"coverage/contrast_state_eligibility.csv")
    atomic_csv(coverage["excluded"],RESULTS/"coverage/excluded_states_with_reasons.csv")
    atomic_csv(coverage["fibro"],RESULTS/"coverage/fibroblast_coverage.csv")
    atomic_csv(definitions,RESULTS/"annotation/marker_definitions.csv")
    atomic_csv(cluster_table,RESULTS/"annotation/cluster_marker_scores.csv")
    atomic_csv(state_scores,RESULTS/"annotation/state_marker_score_summary.csv")
    confidence=adata.obs.groupby(
        ["major_annotation", "annotation_state", "annotation_confidence"], observed=True
    ).size().reset_index(name="cells")
    atomic_csv(confidence,RESULTS/"annotation/annotation_confidence_table.csv")
    cell_table=adata.obs[["mouse_id","gsm","dose_gy","month_post_irradiation","irradiation_group","major_annotation","annotation_state","annotation_confidence","annotation_cluster"]].reset_index()
    atomic_csv(cell_table,RESULTS/"annotation/cell_annotations.csv.gz",compression="gzip")
    sensitivity=[]; keys=list(clusters)
    for i,key_a in enumerate(keys):
        sensitivity.append({"resolution":key_a,"clusters":int(np.unique(clusters[key_a]).size),"primary":key_a==f"leiden_{MARKERS['primary_leiden_resolution']}"})
        for key_b in keys[i+1:]: sensitivity.append({"resolution":f"{key_a}_vs_{key_b}","clusters":np.nan,"primary":False,"adjusted_rand_index":adjusted_rand_score(clusters[key_a],clusters[key_b])})
    atomic_csv(pd.DataFrame(sensitivity),RESULTS/"annotation/resolution_sensitivity.csv")
    save_figures(adata.obs,cluster_table,coverage,display)
    adata.uns["annotation"]={
        "status":"marker_inferred_not_publication_barcode_mapped","condition_blind":True,
        "marker_config_sha256":sha256_file(MARKER_PATH),"primary_cluster_resolution":MARKERS["primary_leiden_resolution"],
        "annotation_pca_role":"clustering and annotation only; not a final ScGeo representation",
        "umap_role":"display_only","normalized_expression_matrix_stored":False,
    }
    temporary=ANNOTATED.with_name(f".{ANNOTATED.name}.tmp.{os.getpid()}.h5ad")
    adata.write_h5ad(temporary,compression="gzip",compression_opts=4)
    backed=ad.read_h5ad(temporary,backed="r")
    try:
        if backed.shape!=adata.shape or "CSRDataset" not in type(backed.X).__name__: raise RuntimeError("Annotated backed validation failed")
    finally: backed.file.close()
    os.replace(temporary,ANNOTATED); record(memory,"annotated_written_validated",started)
    atomic_csv(pd.DataFrame(memory),RESULTS/"annotation/runtime_memory_log.csv")
    counts=adata.obs.groupby(
        ["major_annotation", "annotation_state", "annotation_confidence"], observed=True
    ).size().reset_index(name="cells")
    primary=coverage["eligibility"][coverage["eligibility"].scope.eq("primary") & coverage["eligibility"].eligible_standard]
    report={
        "status":"passed","shape":list(adata.shape),"raw_counts_sparse_format":"csr","raw_counts_dtype":str(adata.X.dtype),
        "annotation_counts":counts.to_dict(orient="records"),"fibroblast_cells":int(adata.obs.major_annotation.eq("fibroblast/stromal").sum()),
        "eligible_primary_states":primary[["contrast","hierarchy","state"]].to_dict(orient="records"),
        "annotated_h5ad":str(ANNOTATED),"sha256":sha256_file(ANNOTATED),"peak_rss_gib":max(x["rss_gib"] for x in memory),
        "runtime_seconds":time.perf_counter()-started,"versions":{name:importlib.metadata.version(name) for name in ["anndata","scikit-learn","igraph","leidenalg","pynndescent","umap-learn"]},
        "final_representation_ensemble_executed":False,
        "future_representation_plan":{"primary":["PCA30","PCA50","scVI without dose/time/condition batch correction"],"sensitivity":["PCA20","diffusion map"],"display_only":["UMAP"]},
    }
    atomic_json(report,RESULTS/"annotation/annotation_validation.json")
    print(json.dumps({k:report[k] for k in ["status","shape","fibroblast_cells","sha256","peak_rss_gib","runtime_seconds"]}))
    return report


def main() -> int:
    parser=argparse.ArgumentParser(); parser.parse_args(); run(); return 0


if __name__=="__main__": raise SystemExit(main())
