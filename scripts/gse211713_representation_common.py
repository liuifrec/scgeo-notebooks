"""Repaired, resource-guarded representation workflow for GSE211713 C6 v2."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
import subprocess
import threading
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/gse211713_representations_v1.json"
SCVI_PYTHON = Path("/home/liuyuchen/micromamba/envs/sc_atac/bin/python")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def configured_paths(config: dict[str, Any] | None = None) -> dict[str, Path]:
    config = load_config() if config is None else config
    output = Path(os.environ.get("SCGEO_GSE211713_C6_V2_OUTPUT_DIR", config["default_output_dir"]))
    if not output.is_absolute():
        output = ROOT / output
    return {
        "compact": Path(os.environ.get("SCGEO_GSE211713_COMPACT_H5AD", config["default_compact_h5ad"])).resolve(),
        "annotated": Path(os.environ.get("SCGEO_GSE211713_ANNOTATED_H5AD", config["default_annotated_h5ad"])).resolve(),
        "representation": Path(os.environ.get("SCGEO_GSE211713_REPRESENTATION_V2_H5AD", config["default_representation_h5ad"])).resolve(),
        "output": output.resolve(),
        "source_repo": Path(os.environ.get("SCGEO_SOURCE_REPO", ROOT.parent / "scgeo")).resolve(),
    }


def ensure_tree(output: Path) -> None:
    for name in ["artifacts", "models", "metadata", "figures", "figure_sources", "alt_text", "execution"]:
        (output / name).mkdir(parents=True, exist_ok=True)


def git_commit(repo: Path) -> str:
    return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()


def git_status(repo: Path) -> list[str]:
    return subprocess.check_output(["git", "-C", str(repo), "status", "--short"], text=True).splitlines()


def active_branch() -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), "branch", "--show-current"], text=True).strip()


def validate_environment(config: dict[str, Any], paths: dict[str, Path]) -> None:
    if active_branch() != config["required_branch"]:
        raise RuntimeError(f"Required branch {config['required_branch']}; observed {active_branch()}")
    if git_commit(paths["source_repo"]) != config["frozen_scgeo_commit"]:
        raise RuntimeError("Frozen ScGeo commit mismatch")
    if git_status(paths["source_repo"]):
        raise RuntimeError("Frozen ScGeo worktree is not clean")
    expected = {
        "compact": config["expected_compact_sha256"],
        "annotated": config["expected_annotated_sha256"],
    }
    for key, checksum in expected.items():
        if not paths[key].is_file():
            raise FileNotFoundError(paths[key])
        observed = sha256(paths[key])
        if observed != checksum:
            raise RuntimeError(f"{key} checksum mismatch: {observed}")
    manifest = ROOT / config["invalidation_manifest"]
    if not manifest.is_file():
        raise RuntimeError("Invalidated-run manifest is required before C6 v2")


def _tmp(path: Path) -> Path:
    return path.with_name(f".{path.name}.tmp.{os.getpid()}")


def atomic_json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _tmp(path)
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _tmp(path)
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _tmp(path)
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(temporary, path)


def atomic_sparse_npz(matrix: sparse.spmatrix, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _tmp(path)
    with temporary.open("wb") as handle:
        sparse.save_npz(handle, matrix.tocsr(), compressed=True)
    os.replace(temporary, path)


def package_versions(names: list[str]) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for name in names:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


class ResourceMonitor:
    def __init__(self, stage: str, config: dict[str, Any]) -> None:
        self.stage = stage
        self.warning = float(config["warning_rss_gib"])
        self.hard_stop = float(config["hard_stop_rss_gib"])
        self.peak_rss_gib = 0.0
        self.peak_gpu_allocated_gib = 0.0
        self.peak_gpu_reserved_gib = 0.0
        self.rows: list[dict[str, Any]] = []
        self._stop = threading.Event()
        self._error: str | None = None

    def _snapshot(self) -> tuple[float, float, float]:
        rss = psutil.Process(os.getpid()).memory_info().rss / 1024**3
        allocated = reserved = 0.0
        try:
            import torch
            if torch.cuda.is_available():
                allocated = torch.cuda.memory_allocated() / 1024**3
                reserved = torch.cuda.memory_reserved() / 1024**3
        except Exception:
            pass
        self.peak_rss_gib = max(self.peak_rss_gib, float(rss))
        self.peak_gpu_allocated_gib = max(self.peak_gpu_allocated_gib, float(allocated))
        self.peak_gpu_reserved_gib = max(self.peak_gpu_reserved_gib, float(reserved))
        if rss > self.hard_stop:
            self._error = f"RSS {rss:.3f} GiB exceeded hard stop {self.hard_stop:.3f} GiB"
            self._stop.set()
        return float(rss), float(allocated), float(reserved)

    def _poll(self) -> None:
        while not self._stop.wait(0.25):
            self._snapshot()

    @contextmanager
    def operation(self, operation: str):
        self._stop.clear()
        self._error = None
        rss0, ga0, gr0 = self._snapshot()
        started = time.perf_counter()
        thread = threading.Thread(target=self._poll, daemon=True)
        thread.start()
        try:
            yield
        finally:
            self._stop.set()
            thread.join(timeout=2)
            rss1, ga1, gr1 = self._snapshot()
            self.rows.append(
                {
                    "stage": self.stage,
                    "operation": operation,
                    "runtime_seconds": time.perf_counter() - started,
                    "rss_before_gib": rss0,
                    "rss_after_gib": rss1,
                    "peak_rss_gib": self.peak_rss_gib,
                    "gpu_allocated_before_gib": ga0,
                    "gpu_allocated_after_gib": ga1,
                    "peak_gpu_allocated_gib": self.peak_gpu_allocated_gib,
                    "gpu_reserved_before_gib": gr0,
                    "gpu_reserved_after_gib": gr1,
                    "peak_gpu_reserved_gib": self.peak_gpu_reserved_gib,
                    "warning_threshold_gib": self.warning,
                    "hard_stop_gib": self.hard_stop,
                    "timestamp_utc": utc_now(),
                }
            )
            if self._error:
                raise MemoryError(self._error)

    def summary(self, runtime_seconds: float) -> dict[str, Any]:
        return {
            "runtime_seconds": float(runtime_seconds),
            "peak_rss_gib": self.peak_rss_gib,
            "peak_gpu_allocated_gib": self.peak_gpu_allocated_gib,
            "peak_gpu_reserved_gib": self.peak_gpu_reserved_gib,
            "warning_threshold_crossed": self.peak_rss_gib >= self.warning,
            "hard_stop_crossed": self.peak_rss_gib >= self.hard_stop,
        }


def load_inputs(paths: dict[str, Path]) -> ad.AnnData:
    compact = ad.read_h5ad(paths["compact"])
    annotated = ad.read_h5ad(paths["annotated"], backed="r")
    try:
        if compact.shape != annotated.shape or not compact.obs_names.equals(annotated.obs_names):
            raise RuntimeError("Compact and annotated inputs are not aligned")
        for column in ["major_annotation", "annotation_state", "annotation_confidence"]:
            if column not in annotated.obs:
                raise RuntimeError(f"Annotated input lacks {column}")
            compact.obs[column] = annotated.obs[column].astype(str).to_numpy()
    finally:
        annotated.file.close()
    if not sparse.isspmatrix_csr(compact.X):
        raise RuntimeError("Frozen counts must remain CSR sparse")
    if compact.X.dtype != np.float32:
        compact.X = compact.X.astype(np.float32)
    return compact


def normalized_hvg_matrix(adata: ad.AnnData, n_hvg: int, target_sum: float) -> tuple[sparse.csr_matrix, list[str]]:
    if "highly_variable" not in adata.var:
        raise RuntimeError("Frozen compact input lacks highly_variable")
    selected = np.flatnonzero(adata.var["highly_variable"].to_numpy(bool))
    if selected.size < n_hvg:
        raise RuntimeError(f"Expected {n_hvg} HVGs; observed {selected.size}")
    selected = selected[:n_hvg]
    matrix = adata[:, selected].X.tocsr(copy=True).astype(np.float32)
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals), where=totals > 0)
    matrix = (sparse.diags(scale.astype(np.float32)) @ matrix).tocsr()
    np.log1p(matrix.data, out=matrix.data)
    return matrix, [str(adata.var_names[i]) for i in selected]


def run_pca_stage() -> dict[str, Any]:
    config, paths = load_config(), configured_paths()
    ensure_tree(paths["output"])
    validate_environment(config, paths)
    started = time.perf_counter()
    monitor = ResourceMonitor("pca", config)
    adata = load_inputs(paths)
    with monitor.operation("sparse_normalization_log1p_pca50"):
        matrix, hvg_names = normalized_hvg_matrix(adata, int(config["pca"]["n_hvg"]), float(config["pca"]["normalization_target_sum"]))
        model = PCA(n_components=50, svd_solver="arpack", random_state=int(config["random_seed"]))
        pca50 = model.fit_transform(matrix).astype(np.float32)
    atomic_npz(paths["output"] / "artifacts/pca_embeddings.npz", X_pca20=pca50[:, :20], X_pca30=pca50[:, :30], X_pca50=pca50)
    atomic_csv(pd.DataFrame({"component": np.arange(1, 51), "explained_variance_ratio": model.explained_variance_ratio_.astype(float)}), paths["output"] / "artifacts/pca_explained_variance.csv")
    report = {
        "status": "passed", "stage": "pca", "dimensions": {"X_pca20": 20, "X_pca30": 30, "X_pca50": 50},
        "hvg_count": len(hvg_names), "hvg_sha256": hashlib.sha256("\n".join(hvg_names).encode()).hexdigest(),
        "algorithm": config["pca"]["algorithm"], **monitor.summary(time.perf_counter() - started),
    }
    atomic_csv(pd.DataFrame(monitor.rows), paths["output"] / "metadata/pca_resource_log.csv")
    atomic_json(report, paths["output"] / "metadata/pca_stage.json")
    del adata, matrix, pca50
    gc.collect()
    return report


def run_diffusion_stage() -> dict[str, Any]:
    import scanpy as sc
    config, paths = load_config(), configured_paths()
    ensure_tree(paths["output"])
    validate_environment(config, paths)
    started = time.perf_counter()
    monitor = ResourceMonitor("diffusion", config)
    adata = load_inputs(paths)
    payload = np.load(paths["output"] / "artifacts/pca_embeddings.npz", allow_pickle=False)
    pca50 = payload["X_pca50"].astype(np.float32)
    work = ad.AnnData(X=sparse.csr_matrix((adata.n_obs, 0), dtype=np.float32), obs=adata.obs.copy())
    work.obsm["X_pca50"] = pca50
    with monitor.operation("pca50_neighbors_diffmap"):
        sc.pp.neighbors(work, use_rep="X_pca50", n_neighbors=int(config["diffusion"]["n_neighbors"]), metric=config["diffusion"]["metric"], random_state=int(config["random_seed"]))
        sc.tl.diffmap(work, n_comps=int(config["diffusion"]["n_components"]))
    diff = np.asarray(work.obsm["X_diffmap"], dtype=np.float32)
    atomic_npz(paths["output"] / "artifacts/diffmap_embedding.npz", X_diffmap=diff)
    atomic_sparse_npz(work.obsp["connectivities"], paths["output"] / "artifacts/pca50_connectivities.npz")
    atomic_sparse_npz(work.obsp["distances"], paths["output"] / "artifacts/pca50_distances.npz")
    report = {
        "status": "passed", "stage": "diffusion", "dimensions": {"X_diffmap": int(diff.shape[1])},
        "n_neighbors": int(config["diffusion"]["n_neighbors"]), "random_seed": int(config["random_seed"]),
        "connectivities_nnz": int(work.obsp["connectivities"].nnz), **monitor.summary(time.perf_counter() - started),
    }
    atomic_csv(pd.DataFrame(monitor.rows), paths["output"] / "metadata/diffusion_resource_log.csv")
    atomic_json(report, paths["output"] / "metadata/diffusion_stage.json")
    return report


def tidy_scvi_history(history: Any) -> pd.DataFrame:
    if not isinstance(history, dict) or not history:
        raise RuntimeError("scVI history is not a non-empty metric dictionary")
    columns: list[pd.Series] = []
    for metric, value in history.items():
        if isinstance(value, pd.DataFrame):
            if value.shape[1] == 1:
                series = value.iloc[:, 0]
                series.name = str(metric)
                columns.append(series)
            else:
                for column in value.columns:
                    series = value[column].copy()
                    series.name = str(metric) if str(column) == str(metric) else f"{metric}_{column}"
                    columns.append(series)
        elif isinstance(value, pd.Series):
            series = value.copy()
            series.name = str(metric)
            columns.append(series)
        elif np.isscalar(value):
            columns.append(pd.Series([value], name=str(metric)))
    if not columns:
        raise RuntimeError("No tabular scVI history metrics were recoverable")
    frame = pd.concat(columns, axis=1)
    frame = frame.loc[:, ~frame.columns.duplicated()].copy()
    frame.insert(0, "epoch", np.arange(1, frame.shape[0] + 1, dtype=int))
    for column in frame.columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    metric_columns = [column for column in frame.columns if column != "epoch"]
    frame = frame.dropna(how="all", subset=metric_columns).reset_index(drop=True)
    frame["epoch"] = np.arange(1, frame.shape[0] + 1, dtype=int)
    if frame.empty or any(not pd.api.types.is_numeric_dtype(frame[column]) for column in frame.columns):
        raise RuntimeError("scVI history validation failed: expected scalar numeric columns")
    return frame


def run_scvi_stage() -> dict[str, Any]:
    import scvi
    import torch
    config, paths = load_config(), configured_paths()
    ensure_tree(paths["output"])
    validate_environment(config, paths)
    if bool(config["scvi"]["require_gpu"]) and not torch.cuda.is_available():
        raise RuntimeError("Approved repair requires an available CUDA GPU")
    started = time.perf_counter()
    monitor = ResourceMonitor("scvi", config)
    adata = load_inputs(paths)
    hvg = np.flatnonzero(adata.var["highly_variable"].to_numpy(bool))[: int(config["scvi"]["n_hvg"])]
    adata = adata[:, hvg].copy()
    adata.X = adata.X.tocsr().astype(np.float32, copy=False)
    scvi.settings.seed = int(config["random_seed"])
    torch.manual_seed(int(config["random_seed"]))
    torch.cuda.manual_seed_all(int(config["random_seed"]))
    torch.set_float32_matmul_precision("high")
    torch.cuda.reset_peak_memory_stats()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(
        adata, n_latent=int(config["scvi"]["n_latent"]), n_hidden=int(config["scvi"]["n_hidden"]),
        n_layers=int(config["scvi"]["n_layers"]), dropout_rate=float(config["scvi"]["dropout_rate"]),
        gene_likelihood=config["scvi"]["gene_likelihood"], dispersion=config["scvi"]["dispersion"],
    )
    with monitor.operation("train_100_epoch_ceiling"):
        model.train(
            max_epochs=int(config["scvi"]["max_epochs"]), batch_size=int(config["scvi"]["batch_size"]),
            early_stopping=bool(config["scvi"]["early_stopping"]), early_stopping_patience=int(config["scvi"]["early_stopping_patience"]),
            train_size=float(config["scvi"]["train_size"]), accelerator="gpu", devices=1,
            datasplitter_kwargs={"num_workers": int(config["scvi"]["num_workers"])}, plan_kwargs={"lr": 1e-3}, check_val_every_n_epoch=1,
        )
    with monitor.operation("save_model_and_latent"):
        latent = model.get_latent_representation().astype(np.float32)
        model_dir = paths["output"] / "models/scvi"
        model.save(model_dir, overwrite=False)
        atomic_npz(
            paths["output"] / "artifacts/scvi_latent.npz",
            X_scvi=latent,
            obs_names=adata.obs_names.astype(str).to_numpy(dtype=str),
        )
    history = tidy_scvi_history(model.history)
    atomic_csv(history, paths["output"] / "artifacts/scvi_training_history.csv")
    parsed = pd.read_csv(paths["output"] / "artifacts/scvi_training_history.csv")
    if parsed.shape != history.shape or any(not pd.api.types.is_numeric_dtype(parsed[column]) for column in parsed.columns):
        raise RuntimeError("Written scVI history failed round-trip numeric validation")
    callback = getattr(getattr(model, "trainer", None), "early_stopping_callback", None)
    stopped_epoch = int(getattr(callback, "stopped_epoch", 0) or 0)
    actual_epochs = int(parsed["epoch"].max())
    peak_allocated = torch.cuda.max_memory_allocated() / 1024**3
    peak_reserved = torch.cuda.max_memory_reserved() / 1024**3
    monitor.peak_gpu_allocated_gib = max(monitor.peak_gpu_allocated_gib, float(peak_allocated))
    monitor.peak_gpu_reserved_gib = max(monitor.peak_gpu_reserved_gib, float(peak_reserved))
    report = {
        "status": "passed", "stage": "scvi", "latent_dimensions": int(latent.shape[1]), "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
        "max_epochs": int(config["scvi"]["max_epochs"]), "actual_epochs": actual_epochs,
        "early_stopping_enabled": bool(config["scvi"]["early_stopping"]), "early_stopped": bool(stopped_epoch > 0 or actual_epochs < int(config["scvi"]["max_epochs"])),
        "early_stopping_stopped_epoch": stopped_epoch, "batch_covariate": None, "accelerator": "gpu",
        "cuda_available": bool(torch.cuda.is_available()), "cuda_device": torch.cuda.get_device_name(0),
        "versions": {"python": os.sys.version.split()[0], "scvi-tools": scvi.__version__, "torch": torch.__version__, "cuda_runtime": torch.version.cuda, **package_versions(["anndata", "scanpy", "lightning", "pytorch-lightning"])},
        "history_rows": int(parsed.shape[0]), "history_columns": parsed.columns.tolist(), **monitor.summary(time.perf_counter() - started),
    }
    atomic_csv(pd.DataFrame(monitor.rows), paths["output"] / "metadata/scvi_resource_log.csv")
    atomic_json(report, paths["output"] / "metadata/scvi_stage.json")
    return report


def _stratified_indices(labels: np.ndarray, maximum: int, seed: int) -> np.ndarray:
    if labels.size <= maximum:
        return np.arange(labels.size)
    rng = np.random.default_rng(seed)
    chosen: list[int] = []
    for label in pd.unique(labels):
        idx = np.flatnonzero(labels == label)
        quota = min(idx.size, max(1, round(maximum * idx.size / labels.size)))
        chosen.extend(rng.choice(idx, size=quota, replace=False).tolist())
    chosen_array = np.array(sorted(set(chosen)), dtype=int)
    if chosen_array.size > maximum:
        chosen_array = np.sort(rng.choice(chosen_array, size=maximum, replace=False))
    return chosen_array


def _knn(coords: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    model = NearestNeighbors(n_neighbors=min(k + 1, coords.shape[0]), metric="euclidean", n_jobs=1).fit(coords)
    distances, indices = model.kneighbors(coords)
    return indices[:, 1:], distances[:, 1:]


def _largest_component(indices: np.ndarray) -> float:
    n, k = indices.shape
    rows = np.repeat(np.arange(n), k)
    cols = indices.reshape(-1)
    graph = sparse.coo_matrix((np.ones(rows.size * 2), (np.r_[rows, cols], np.r_[cols, rows])), shape=(n, n)).tocsr()
    _, labels = connected_components(graph, directed=False, return_labels=True)
    return float(np.bincount(labels).max() / labels.size)


def _overlap(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean([len(set(x).intersection(y)) / a.shape[1] for x, y in zip(a, b, strict=True)]))


def _distortion(reference: np.ndarray, candidate: np.ndarray) -> float:
    ref_med, cand_med = np.median(reference[:, 0]), np.median(candidate[:, 0])
    ratio = np.log((candidate[:, 0] / max(cand_med, 1e-8) + 1e-8) / (reference[:, 0] / max(ref_med, 1e-8) + 1e-8))
    return float(np.median(np.abs(ratio)))


def standardize_obs(adata: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    obs = adata.obs.copy()
    major_map = {
        "endothelial": "Endothelial", "epithelial": "Epithelial", "fibroblast/stromal": "Fibroblast/stromal",
        "lymphoid": "Lymphoid", "myeloid": "Myeloid", "proliferating": "Proliferating",
    }
    obs["state_major"] = obs["major_annotation"].astype(str).map(major_map)
    reverse_fibro = {source: display for display, source in config["source_fibroblast_labels"].items()}
    obs["fibroblast_state"] = obs["annotation_state"].astype(str).map(reverse_fibro)
    obs["sample_key"] = obs["mouse_id"].astype(str)
    return obs


def _quality_metrics(rep: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    keys = ["X_pca20", "X_pca30", "X_pca50", "X_diffmap", "X_scvi"]
    k = int(config["quality"]["n_neighbors"])
    global_labels = rep.obs["annotation_state"].astype(str).to_numpy()
    global_idx = _stratified_indices(global_labels, int(config["quality"]["silhouette_sample_size"]), int(config["random_seed"]))
    major_mask = rep.obs["state_major"].isin(config["approved_major_states"]).to_numpy()
    major_positions = np.flatnonzero(major_mask)
    major_local = _stratified_indices(rep.obs.iloc[major_positions]["state_major"].astype(str).to_numpy(), int(config["quality"]["silhouette_sample_size"]), int(config["random_seed"]))
    major_idx = major_positions[major_local]
    fibro_mask = rep.obs["fibroblast_state"].isin(config["approved_fibroblast_states"]).to_numpy()
    fibro_positions = np.flatnonzero(fibro_mask)
    fibro_local = _stratified_indices(rep.obs.iloc[fibro_positions]["fibroblast_state"].astype(str).to_numpy(), int(config["quality"]["fibroblast_silhouette_sample_size"]), int(config["random_seed"]))
    fibro_idx = fibro_positions[fibro_local]
    reference_idx, reference_dist = _knn(np.asarray(rep.obsm["X_pca50"])[global_idx], k)
    rows = []
    roles = config["representation_roles"]
    for key in keys:
        full = np.asarray(rep.obsm[key], dtype=np.float32)
        subset = full[global_idx]
        indices, distances = _knn(subset, k)
        mouse = rep.obs.iloc[global_idx]["mouse_id"].astype(str).to_numpy()
        same_mouse = float(np.mean([np.mean(mouse[row] == mouse[i]) for i, row in enumerate(indices)]))
        major_sil = float(silhouette_score(full[major_idx], rep.obs.iloc[major_idx]["state_major"].astype(str), metric="euclidean"))
        fibro_sil = float(silhouette_score(full[fibro_idx], rep.obs.iloc[fibro_idx]["fibroblast_state"].astype(str), metric="euclidean"))
        dose_sil = float(silhouette_score(subset, rep.obs.iloc[global_idx]["dose_gy"].astype(str), metric="euclidean"))
        time_sil = float(silhouette_score(subset, rep.obs.iloc[global_idx]["month_post_irradiation"].astype(str), metric="euclidean"))
        role = "primary" if key in roles["primary"] else "dimensional_sensitivity" if key in roles["dimensional_sensitivity"] else "exploratory_sensitivity"
        rows.append({
            "representation": key, "role": role, "dimensions": int(full.shape[1]), "nonfinite_coordinates": int((~np.isfinite(full)).sum()),
            "mouse_neighborhood_fraction": same_mouse, "major_state_silhouette": major_sil, "fibroblast_subtype_silhouette": fibro_sil,
            "dose_silhouette_descriptive": dose_sil, "time_silhouette_descriptive": time_sil,
            "neighborhood_overlap_vs_pca50": 1.0 if key == "X_pca50" else _overlap(indices, reference_idx),
            "local_distortion_vs_pca50": 0.0 if key == "X_pca50" else _distortion(reference_dist, distances),
            "largest_connected_component_fraction": _largest_component(indices), "quality_sample_cells": int(global_idx.size),
            "major_metric_cells": int(major_idx.size), "fibroblast_metric_cells": int(fibro_idx.size),
        })
    return pd.DataFrame(rows)


def run_assembly_quality_stage() -> dict[str, Any]:
    import scanpy as sc
    config, paths = load_config(), configured_paths()
    ensure_tree(paths["output"])
    validate_environment(config, paths)
    started = time.perf_counter()
    monitor = ResourceMonitor("assembly_quality", config)
    adata = load_inputs(paths)
    obs = standardize_obs(adata, config)
    pca = np.load(paths["output"] / "artifacts/pca_embeddings.npz", allow_pickle=False)
    diff = np.load(paths["output"] / "artifacts/diffmap_embedding.npz", allow_pickle=False)["X_diffmap"].astype(np.float32)
    # Repair-pass compatibility: attempt 2 wrote only obs_names with NumPy object
    # dtype. The artifact is locally generated and checksummed; embeddings are
    # never modified. Future writes use fixed-width Unicode above.
    scvi_payload = np.load(paths["output"] / "artifacts/scvi_latent.npz", allow_pickle=True)
    if not np.array_equal(scvi_payload["obs_names"].astype(str), adata.obs_names.to_numpy().astype(str)):
        raise RuntimeError("scVI latent observation order mismatch")
    rep = ad.AnnData(X=sparse.csr_matrix((adata.n_obs, 0), dtype=np.float32), obs=obs, var=pd.DataFrame(index=pd.Index([], dtype=str)))
    for key in ["X_pca20", "X_pca30", "X_pca50"]:
        rep.obsm[key] = pca[key].astype(np.float32)
    rep.obsm["X_diffmap"] = diff
    rep.obsm["X_scvi"] = scvi_payload["X_scvi"].astype(np.float32)
    rep.obsp["connectivities_pca50"] = sparse.load_npz(paths["output"] / "artifacts/pca50_connectivities.npz").tocsr()
    rep.obsp["distances_pca50"] = sparse.load_npz(paths["output"] / "artifacts/pca50_distances.npz").tocsr()
    display = ad.AnnData(X=sparse.csr_matrix((rep.n_obs, 0), dtype=np.float32), obs=rep.obs.copy())
    display.obsm["X_pca50"] = rep.obsm["X_pca50"]
    with monitor.operation("display_only_umap"):
        sc.pp.neighbors(display, use_rep="X_pca50", n_neighbors=int(config["quality"]["umap_n_neighbors"]), metric=config["quality"]["umap_metric"], random_state=int(config["random_seed"]))
        sc.tl.umap(display, min_dist=float(config["quality"]["umap_min_dist"]), random_state=int(config["random_seed"]))
    rep.obsm["X_umap_display_only"] = np.asarray(display.obsm["X_umap"], dtype=np.float32)
    with monitor.operation("representation_quality"):
        quality = _quality_metrics(rep, config)
    stage_reports = {name: json.loads((paths["output"] / f"metadata/{name}_stage.json").read_text()) for name in ["pca", "diffusion", "scvi"]}
    quality["runtime_seconds"] = quality["representation"].map({
        "X_pca20": stage_reports["pca"]["runtime_seconds"], "X_pca30": stage_reports["pca"]["runtime_seconds"], "X_pca50": stage_reports["pca"]["runtime_seconds"],
        "X_diffmap": stage_reports["diffusion"]["runtime_seconds"], "X_scvi": stage_reports["scvi"]["runtime_seconds"],
    })
    quality["peak_rss_gib"] = quality["representation"].map({
        "X_pca20": stage_reports["pca"]["peak_rss_gib"], "X_pca30": stage_reports["pca"]["peak_rss_gib"], "X_pca50": stage_reports["pca"]["peak_rss_gib"],
        "X_diffmap": stage_reports["diffusion"]["peak_rss_gib"], "X_scvi": stage_reports["scvi"]["peak_rss_gib"],
    })
    quality["peak_gpu_allocated_gib"] = quality["representation"].map({"X_pca20": 0.0, "X_pca30": 0.0, "X_pca50": 0.0, "X_diffmap": 0.0, "X_scvi": stage_reports["scvi"]["peak_gpu_allocated_gib"]})
    quality["peak_gpu_reserved_gib"] = quality["representation"].map({"X_pca20": 0.0, "X_pca30": 0.0, "X_pca50": 0.0, "X_diffmap": 0.0, "X_scvi": stage_reports["scvi"]["peak_gpu_reserved_gib"]})
    atomic_csv(quality, paths["output"] / "metadata/representation_metrics.csv")
    rep.uns["representation_v2"] = {
        "status": "passed", "sample_key": "mouse_id", "biological_replicate_unit": "mouse/GSM",
        "roles": config["representation_roles"], "batch_covariate": None, "umap_use": "display_only",
        "compact_sha256": sha256(paths["compact"]), "annotated_sha256": sha256(paths["annotated"]),
        "frozen_scgeo_commit": config["frozen_scgeo_commit"], "config_sha256": sha256(CONFIG_PATH),
    }
    temporary = _tmp(paths["representation"])
    paths["representation"].parent.mkdir(parents=True, exist_ok=True)
    with monitor.operation("atomic_h5ad_write"):
        rep.write_h5ad(temporary, compression="gzip", compression_opts=4)
        backed = ad.read_h5ad(temporary, backed="r")
        try:
            if backed.shape != rep.shape:
                raise RuntimeError("Backed representation shape mismatch")
            required = set(sum(config["representation_roles"].values(), []))
            if required.difference(backed.obsm.keys()):
                raise RuntimeError("Backed representation lacks required embeddings")
        finally:
            backed.file.close()
        os.replace(temporary, paths["representation"])
    checksum = sha256(paths["representation"])
    report = {
        "status": "passed", "stage": "assembly_quality", "representation_h5ad": str(paths["representation"]),
        "representation_sha256": checksum, "shape": list(rep.shape), "all_finite": bool((quality["nonfinite_coordinates"] == 0).all()),
        "all_five_representations_persisted": set(config["representation_roles"]["primary"] + config["representation_roles"]["dimensional_sensitivity"] + config["representation_roles"]["exploratory_sensitivity"]).issubset(rep.obsm.keys()),
        **monitor.summary(time.perf_counter() - started),
    }
    atomic_csv(pd.DataFrame(monitor.rows), paths["output"] / "metadata/assembly_quality_resource_log.csv")
    atomic_json(report, paths["output"] / "metadata/assembly_quality_stage.json")
    atomic_json({"status": "passed", "representation_sha256": checksum, "metrics_rows": int(quality.shape[0]), "timestamp_utc": utc_now()}, paths["output"] / "metadata/c6_summary.json")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["pca", "diffusion", "scvi", "assembly_quality"])
    args = parser.parse_args()
    result = {"pca": run_pca_stage, "diffusion": run_diffusion_stage, "scvi": run_scvi_stage, "assembly_quality": run_assembly_quality_stage}[args.stage]()
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
