"""Resource-guarded, descriptive-only representation helpers for GSE249479.

Raw expression is read from the immutable compact H5AD. Only the 3,000-HVG
sparse count submatrix is copied for PCA/scVI. No complete expression matrix is
densified, and every heavy representation is designed to run in a fresh process.
"""

from __future__ import annotations

import _thread
import gc
import hashlib
import importlib.metadata as im
import json
import os
import platform
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/gse249479_representations_v1.json"


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def configured_paths(config: dict[str, Any] | None = None) -> dict[str, Path]:
    config = load_config() if config is None else config
    output = Path(os.environ.get("SCGEO_GSE249479_OUTPUT_DIR", config["default_output_dir"]))
    if not output.is_absolute():
        output = ROOT / output
    return {
        "compact": Path(os.environ.get("SCGEO_GSE249479_COMPACT_H5AD", config["default_compact_h5ad"])).resolve(),
        "representation": Path(os.environ.get("SCGEO_GSE249479_REPRESENTATION_H5AD", config["default_representation_h5ad"])).resolve(),
        "output": output.resolve(),
        "source_repo": Path(os.environ.get("SCGEO_SOURCE_REPO", ROOT.parent / "scgeo")).resolve(),
    }


def ensure_tree(output: Path) -> None:
    for name in ["representations", "audit", "figure_sources", "figures", "alt_text", "metadata", "version_records", "execution", "executed_notebooks", "models"]:
        (output / name).mkdir(parents=True, exist_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def atomic_json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(frame: pd.DataFrame, path: Path, **kwargs: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".csv.gz" if path.name.endswith(".csv.gz") else ".csv"
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}{suffix}")
    frame.to_csv(tmp, index=False, **kwargs)
    os.replace(tmp, path)


def atomic_npz(path: Path, **arrays: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


def atomic_sparse_npz(matrix: sparse.spmatrix, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}.npz")
    sparse.save_npz(tmp, matrix, compressed=True)
    os.replace(tmp, path)


def git_commit(path: Path) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()
    except Exception:
        return None


def active_branch() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "branch", "--show-current"], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()
    except Exception:
        return None


def gpu_snapshot() -> dict[str, Any]:
    result = {"gpu_available": False, "gpu_memory_allocated_gib": 0.0, "gpu_memory_reserved_gib": 0.0, "gpu_name": None}
    try:
        import torch

        available = bool(torch.cuda.is_available())
        result["gpu_available"] = available
        if available:
            result.update({
                "gpu_memory_allocated_gib": float(torch.cuda.memory_allocated() / 1024**3),
                "gpu_memory_reserved_gib": float(torch.cuda.memory_reserved() / 1024**3),
                "gpu_name": torch.cuda.get_device_name(0),
            })
    except Exception:
        pass
    return result


class ResourceLog:
    def __init__(self, stage: str, config: dict[str, Any]):
        self.stage = stage
        self.warning = float(config["warning_rss_gib"])
        self.hard_stop = float(config["hard_stop_rss_gib"])
        self.rows: list[dict[str, Any]] = []
        self.peak_rss_gib = psutil.Process().memory_info().rss / 1024**3
        self.peak_gpu_allocated_gib = 0.0
        self.peak_gpu_reserved_gib = 0.0
        self._stop = threading.Event()
        self._hard_exceeded = False

    def snapshot(self) -> dict[str, Any]:
        rss = psutil.Process().memory_info().rss / 1024**3
        gpu = gpu_snapshot()
        self.peak_rss_gib = max(self.peak_rss_gib, rss)
        self.peak_gpu_allocated_gib = max(self.peak_gpu_allocated_gib, float(gpu["gpu_memory_allocated_gib"]))
        self.peak_gpu_reserved_gib = max(self.peak_gpu_reserved_gib, float(gpu["gpu_memory_reserved_gib"]))
        return {"rss_gib": rss, **gpu}

    def _monitor(self) -> None:
        while not self._stop.wait(0.5):
            snap = self.snapshot()
            if snap["rss_gib"] >= self.hard_stop:
                self._hard_exceeded = True
                _thread.interrupt_main()
                return

    @contextmanager
    def operation(self, name: str):
        close_and_collect()
        before = self.snapshot()
        self._stop.clear()
        monitor = threading.Thread(target=self._monitor, daemon=True)
        monitor.start()
        started = time.perf_counter()
        status = "passed"
        try:
            yield
        except BaseException:
            status = "hard_stop" if self._hard_exceeded else "failed"
            raise
        finally:
            self._stop.set()
            monitor.join(timeout=2)
            close_and_collect()
            after = self.snapshot()
            self.rows.append({
                "stage": self.stage, "operation": name, "status": status,
                "runtime_seconds": time.perf_counter() - started,
                "rss_before_gib": before["rss_gib"], "rss_after_gib": after["rss_gib"],
                "peak_stage_rss_gib": self.peak_rss_gib,
                "gpu_available": after["gpu_available"], "gpu_name": after["gpu_name"],
                "gpu_allocated_before_gib": before["gpu_memory_allocated_gib"],
                "gpu_allocated_after_gib": after["gpu_memory_allocated_gib"],
                "peak_stage_gpu_allocated_gib": self.peak_gpu_allocated_gib,
                "peak_stage_gpu_reserved_gib": self.peak_gpu_reserved_gib,
                "warning_rss_gib": self.warning, "hard_stop_rss_gib": self.hard_stop,
                "timestamp_utc": utc_now(),
            })
            if self._hard_exceeded:
                raise MemoryError(f"Stage {self.stage} exceeded {self.hard_stop:.1f} GiB RSS")

    def write(self, output: Path, reset_combined: bool = False) -> None:
        frame = pd.DataFrame(self.rows)
        atomic_csv(frame, output / "audit" / f"{self.stage}_resource_log.csv")
        combined = output / "audit" / "03_representation_resource_log.csv"
        if combined.exists() and not reset_combined:
            frame = pd.concat([pd.read_csv(combined), frame], ignore_index=True)
        atomic_csv(frame, combined)


def close_and_collect() -> None:
    try:
        import matplotlib.pyplot as plt

        plt.close("all")
    except Exception:
        pass
    gc.collect()


def validate_environment(config: dict[str, Any], p: dict[str, Path]) -> None:
    if active_branch() != config["required_branch"]:
        raise RuntimeError(f"Required branch {config['required_branch']!r}; observed {active_branch()!r}")
    if git_commit(p["source_repo"]) != config["frozen_scgeo_commit"]:
        raise RuntimeError("Frozen ScGeo commit mismatch")
    if not p["compact"].is_file():
        raise FileNotFoundError(p["compact"])


def package_versions(extra: list[str] | None = None) -> dict[str, str | None]:
    names = ["python", "anndata", "numpy", "pandas", "scipy", "scikit-learn", "scanpy", "umap-learn", "psutil", "nbclient", "nbformat"]
    if extra:
        names.extend(extra)
    result = {}
    for name in names:
        if name == "python":
            result[name] = platform.python_version()
        else:
            try:
                result[name] = im.version(name)
            except im.PackageNotFoundError:
                result[name] = None
    return result


def write_stage_metadata(stage: str, config: dict[str, Any], p: dict[str, Path], resources: ResourceLog, payload: dict[str, Any], extra_packages: list[str] | None = None) -> None:
    metadata = {
        "stage": stage, "timestamp_utc": utc_now(), "inference_status": "descriptive_only",
        "compact_h5ad": str(p["compact"]), "compact_sha256": sha256(p["compact"]),
        "notebook_repository_commit": git_commit(ROOT), "frozen_scgeo_commit": git_commit(p["source_repo"]),
        "scgeo_modified": False, "forbidden_methods_run": [],
        "peak_cpu_rss_gib": resources.peak_rss_gib,
        "peak_gpu_allocated_gib": resources.peak_gpu_allocated_gib,
        "peak_gpu_reserved_gib": resources.peak_gpu_reserved_gib,
        "packages": package_versions(extra_packages),
    }
    metadata.update(payload)
    atomic_json(metadata, p["output"] / "metadata" / f"{stage}_metadata.json")
    atomic_json({"stage": stage, "packages": metadata["packages"]}, p["output"] / "version_records" / f"{stage}_versions.json")


def load_compact_hvg(p: dict[str, Path], expected_hvg: int) -> tuple[ad.AnnData, sparse.csr_matrix, np.ndarray]:
    adata = ad.read_h5ad(p["compact"])
    if not sparse.issparse(adata.X):
        raise RuntimeError("Compact X must remain sparse")
    if adata.uns.get("dataset_b_inference_status") != "descriptive_only":
        raise RuntimeError("Compact object lost descriptive_only status")
    hvg = adata.var["highly_variable"].to_numpy(dtype=bool)
    if int(hvg.sum()) != expected_hvg:
        raise RuntimeError(f"Expected {expected_hvg} HVGs; observed {int(hvg.sum())}")
    X = adata.X[:, hvg].tocsr().astype(np.float32, copy=True)
    return adata, X, hvg


def obs_name_hash(names: pd.Index) -> str:
    return hashlib.sha256("\n".join(names.astype(str)).encode()).hexdigest()


def run_pca() -> dict[str, Any]:
    from sklearn.decomposition import TruncatedSVD
    from sklearn.utils.sparsefuncs import inplace_row_scale

    config, p = load_config(), configured_paths()
    validate_environment(config, p)
    ensure_tree(p["output"])
    stage = "03a_pca_representations"
    resources = ResourceLog(stage, config)
    source_sha_before = sha256(p["compact"])
    started = time.perf_counter()
    with resources.operation("load_and_subset_3000_hvgs"):
        adata, X, hvg = load_compact_hvg(p, int(config["pca"]["n_hvg"]))
        totals = adata.obs["total_counts"].to_numpy(dtype=np.float32)

    with resources.operation("sparse_normalize_log1p_hvgs"):
        scale = np.divide(np.float32(config["pca"]["normalization_target_sum"]), totals, out=np.zeros_like(totals), where=totals > 0)
        inplace_row_scale(X, scale)
        X.data = np.log1p(X.data).astype(np.float32, copy=False)

    seed = int(config["random_seed"])
    n_components = int(config["pca"]["n_components"])
    with resources.operation("shared_sparse_pca50"):
        model = TruncatedSVD(n_components=n_components, algorithm="randomized", n_iter=int(config["pca"]["n_iter"]), random_state=seed)
        embedding = model.fit_transform(X).astype(np.float32)
        loadings = model.components_.astype(np.float32)
        variance = model.explained_variance_.astype(np.float32)
        variance_ratio = model.explained_variance_ratio_.astype(np.float32)
        singular = model.singular_values_.astype(np.float32)

    with resources.operation("component_seed_stability"):
        model2 = TruncatedSVD(n_components=n_components, algorithm="randomized", n_iter=int(config["pca"]["n_iter"]), random_state=seed + int(config["pca"]["stability_seed_offset"]))
        embedding2 = model2.fit_transform(X).astype(np.float32)
        corr = np.corrcoef(np.vstack([loadings, model2.components_.astype(np.float32)]))[:n_components, n_components:]
        rows, cols = linear_sum_assignment(-np.abs(corr))
        stability = pd.DataFrame({"primary_component": rows + 1, "matched_component": cols + 1, "absolute_loading_correlation": np.abs(corr[rows, cols])})
        embedding_corr = [float(abs(np.corrcoef(embedding[:, i], embedding2[:, cols[np.where(rows == i)[0][0]]])[0, 1])) for i in rows]
        stability["absolute_score_correlation"] = embedding_corr
        del embedding2, model2, corr

    checkpoint = p["output"] / "representations" / "03a_pca_representations.npz"
    with resources.operation("save_pca_checkpoint"):
        atomic_npz(
            checkpoint, X_pca50=embedding, X_pca30=embedding[:, :30], X_pca20=embedding[:, :20],
            loadings=loadings, explained_variance=variance, explained_variance_ratio=variance_ratio,
            singular_values=singular, hvg_names=np.array(adata.var_names[hvg].astype(str).tolist(), dtype=str),
            obs_names=np.array(adata.obs_names.astype(str).tolist(), dtype=str),
        )
        atomic_csv(pd.DataFrame({"component": np.arange(1, n_components + 1), "explained_variance": variance, "explained_variance_ratio": variance_ratio, "cumulative_explained_variance_ratio": np.cumsum(variance_ratio)}), p["output"] / "figure_sources" / "03a_pca_explained_variance.csv")
        atomic_csv(stability, p["output"] / "audit" / "03a_pca_component_stability.csv")

    summary = {
        "status": "passed", "algorithm": config["pca"]["algorithm"],
        "dimensions": {"X_pca20": [adata.n_obs, 20], "X_pca30": [adata.n_obs, 30], "X_pca50": [adata.n_obs, 50]},
        "shared_basis": True, "independent_representations": False,
        "cumulative_explained_variance_ratio_20": float(variance_ratio[:20].sum()),
        "cumulative_explained_variance_ratio_30": float(variance_ratio[:30].sum()),
        "cumulative_explained_variance_ratio_50": float(variance_ratio.sum()),
        "median_component_loading_stability": float(stability["absolute_loading_correlation"].median()),
        "minimum_component_loading_stability": float(stability["absolute_loading_correlation"].min()),
        "checkpoint": str(checkpoint), "checkpoint_sha256": sha256(checkpoint),
        "runtime_seconds": time.perf_counter() - started,
        "source_sha256_before": source_sha_before, "source_sha256_after": sha256(p["compact"]),
        "raw_counts_modified": False, "inference_status": "descriptive_only",
    }
    atomic_json(summary, p["output"] / "audit" / "03a_pca_summary.json")
    resources.write(p["output"], reset_combined=True)
    write_stage_metadata(stage, config, p, resources, summary)
    del adata, X, embedding, loadings
    close_and_collect()
    return summary


def neighbor_jaccard(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.empty(a.shape[0], dtype=np.float32)
    for i in range(a.shape[0]):
        sa, sb = set(a[i].tolist()), set(b[i].tolist())
        out[i] = len(sa & sb) / max(1, len(sa | sb))
    return out


def run_diffusion() -> dict[str, Any]:
    import scanpy as sc
    from sklearn.neighbors import NearestNeighbors

    config, p = load_config(), configured_paths()
    validate_environment(config, p)
    ensure_tree(p["output"])
    stage = "03b_diffusion_representation"
    resources = ResourceLog(stage, config)
    source_sha_before = sha256(p["compact"])
    started = time.perf_counter()
    pca_path = p["output"] / "representations" / "03a_pca_representations.npz"
    with resources.operation("load_pca_checkpoint"):
        with np.load(pca_path) as data:
            pca50 = data["X_pca50"].astype(np.float32)
            pca30 = data["X_pca30"].astype(np.float32)
            names = data["obs_names"].astype(str)
    n, k = pca50.shape[0], int(config["diffusion"]["n_neighbors"])
    estimated_graph_bytes = int(2 * (n * k * (4 + 4) + (n + 1) * 4))
    if estimated_graph_bytes / 1024**3 >= float(config["warning_rss_gib"]):
        raise MemoryError("Estimated diffusion graph storage exceeds warning threshold")

    with resources.operation("sparse_neighbor_graph"):
        work = ad.AnnData(X=pca50)
        work.obs_names = names
        sc.pp.neighbors(work, n_neighbors=k, metric=config["diffusion"]["metric"], method="umap", random_state=int(config["random_seed"]))
        connectivities = work.obsp["connectivities"].tocsr().astype(np.float32)
        distances = work.obsp["distances"].tocsr().astype(np.float32)
        n_components_graph, component_labels = connected_components(connectivities, directed=False)

    with resources.operation("diffusion_components"):
        sc.tl.diffmap(work, n_comps=int(config["diffusion"]["n_components"]), random_state=int(config["random_seed"]))
        diffmap = work.obsm["X_diffmap"].astype(np.float32)

    with resources.operation("neighbor_graph_stability"):
        nn50 = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=1).fit(pca50).kneighbors(return_distance=False)[:, 1:]
        nn30 = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=1).fit(pca30).kneighbors(return_distance=False)[:, 1:]
        jaccard = neighbor_jaccard(nn50, nn30)

    checkpoint = p["output"] / "representations" / "03b_diffusion_representation.npz"
    with resources.operation("save_diffusion_checkpoint"):
        atomic_npz(checkpoint, X_diffmap=diffmap, obs_names=names, graph_component=component_labels.astype(np.int32))
        if bool(config["diffusion"]["retain_sparse_graphs"]):
            atomic_sparse_npz(connectivities, p["output"] / "representations" / "03b_connectivities.npz")
            atomic_sparse_npz(distances, p["output"] / "representations" / "03b_distances.npz")

    summary = {
        "status": "passed", "dimensions": [int(diffmap.shape[0]), int(diffmap.shape[1])],
        "n_neighbors": k, "estimated_graph_bytes": estimated_graph_bytes,
        "connectivities_nnz": int(connectivities.nnz), "distances_nnz": int(distances.nnz),
        "disconnected_components": int(n_components_graph),
        "largest_component_fraction": float(np.bincount(component_labels).max() / len(component_labels)),
        "pca30_vs_pca50_neighbor_jaccard_median": float(np.median(jaccard)),
        "pca30_vs_pca50_neighbor_jaccard_q05": float(np.quantile(jaccard, 0.05)),
        "checkpoint": str(checkpoint), "checkpoint_sha256": sha256(checkpoint),
        "runtime_seconds": time.perf_counter() - started,
        "source_sha256_before": source_sha_before, "source_sha256_after": sha256(p["compact"]),
        "inference_status": "descriptive_only",
    }
    atomic_json(summary, p["output"] / "audit" / "03b_diffusion_summary.json")
    resources.write(p["output"])
    write_stage_metadata(stage, config, p, resources, summary)
    del work, pca50, pca30, diffmap, connectivities, distances
    close_and_collect()
    return summary


def history_frame(history: dict[str, Any]) -> pd.DataFrame:
    pieces = []
    for metric, values in history.items():
        frame = pd.DataFrame(values).reset_index().rename(columns={"index": "epoch"})
        value_cols = [c for c in frame.columns if c != "epoch"]
        for col in value_cols:
            pieces.append(pd.DataFrame({"epoch": frame["epoch"], "metric": metric if len(value_cols) == 1 else f"{metric}:{col}", "value": frame[col]}))
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=["epoch", "metric", "value"])


def run_scvi() -> dict[str, Any]:
    import scvi
    import torch

    config, p = load_config(), configured_paths()
    validate_environment(config, p)
    ensure_tree(p["output"])
    stage = "03c_scvi_representation"
    resources = ResourceLog(stage, config)
    source_sha_before = sha256(p["compact"])
    started = time.perf_counter()
    seed = int(config["random_seed"])
    np.random.seed(seed)
    torch.manual_seed(seed)
    scvi.settings.seed = seed
    sc = config["scvi"]

    with resources.operation("load_raw_hvg_counts"):
        source, X, hvg = load_compact_hvg(p, int(config["pca"]["n_hvg"]))
        work = ad.AnnData(X=X, obs=pd.DataFrame(index=source.obs_names.copy()), var=pd.DataFrame(index=source.var_names[hvg].copy()))

    with resources.operation("configure_scvi_no_batch"):
        scvi.model.SCVI.setup_anndata(work)
        model = scvi.model.SCVI(
            work, n_hidden=int(sc["n_hidden"]), n_latent=int(sc["n_latent"]), n_layers=int(sc["n_layers"]),
            dropout_rate=float(sc["dropout_rate"]), dispersion=sc["dispersion"], gene_likelihood=sc["gene_likelihood"],
        )
        use_gpu = bool(torch.cuda.is_available())
        accelerator = "gpu" if use_gpu else "cpu"
        devices: Any = 1

    with resources.operation("train_scvi"):
        model.train(
            max_epochs=int(sc["max_epochs"]), accelerator=accelerator, devices=devices,
            train_size=float(sc["train_size"]), batch_size=int(sc["batch_size"]),
            early_stopping=bool(sc["early_stopping"]),
            datasplitter_kwargs={"num_workers": int(sc["num_workers"])},
            early_stopping_patience=int(sc["early_stopping_patience"]),
            enable_progress_bar=False,
        )

    with resources.operation("extract_scvi_latent_and_history"):
        latent = model.get_latent_representation(batch_size=int(sc["batch_size"])).astype(np.float32)
        history = history_frame(model.history)
        trained_epochs = int(history["epoch"].max() + 1) if not history.empty else None

    latent_path = p["output"] / "representations" / "03c_scvi_representation.npz"
    model_path = p["output"] / "models" / "03c_scvi_state_dict.pt"
    with resources.operation("save_scvi_checkpoint"):
        atomic_npz(latent_path, X_scvi=latent, obs_names=np.array(source.obs_names.astype(str).tolist(), dtype=str))
        tmp_model = model_path.with_name(f".{model_path.name}.tmp.{os.getpid()}")
        torch.save(model.module.state_dict(), tmp_model)
        os.replace(tmp_model, model_path)
        atomic_csv(history, p["output"] / "figure_sources" / "03c_scvi_training_history.csv")
        atomic_json({
            "model_class": "scvi.model.SCVI", "n_hidden": int(sc["n_hidden"]), "n_latent": int(sc["n_latent"]),
            "n_layers": int(sc["n_layers"]), "dropout_rate": float(sc["dropout_rate"]),
            "dispersion": sc["dispersion"], "gene_likelihood": sc["gene_likelihood"],
            "batch_covariate": None, "batch_size": int(sc["batch_size"]), "num_workers": int(sc["num_workers"]),
            "seed": seed, "accelerator": accelerator, "trained_epochs": trained_epochs,
            "state_dict": str(model_path), "state_dict_sha256": sha256(model_path),
        }, p["output"] / "models" / "03c_scvi_model_parameters.json")

    summary = {
        "status": "passed", "dimensions": [int(latent.shape[0]), int(latent.shape[1])],
        "accelerator": accelerator, "gpu_available": use_gpu, "batch_covariate": None,
        "batch_size": int(sc["batch_size"]), "num_workers": int(sc["num_workers"]),
        "trained_epochs": trained_epochs, "early_stopping_requested": bool(sc["early_stopping"]),
        "checkpoint": str(latent_path), "checkpoint_sha256": sha256(latent_path),
        "model_state_sha256": sha256(model_path), "runtime_seconds": time.perf_counter() - started,
        "source_sha256_before": source_sha_before, "source_sha256_after": sha256(p["compact"]),
        "dense_normalized_expression_stored": False, "inference_status": "descriptive_only",
    }
    atomic_json(summary, p["output"] / "audit" / "03c_scvi_summary.json")
    resources.write(p["output"])
    write_stage_metadata(stage, config, p, resources, summary, ["scvi-tools", "torch", "lightning"])
    del source, work, X, model, latent
    close_and_collect()
    return summary


def knn(embedding: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.neighbors import NearestNeighbors

    distances, indices = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=1).fit(embedding).kneighbors()
    return distances[:, 1:].astype(np.float32), indices[:, 1:].astype(np.int32)


def same_neighbor_fraction(values: np.ndarray, indices: np.ndarray, valid: np.ndarray | None = None) -> float:
    if valid is None:
        valid = np.ones(len(values), dtype=bool)
    scores = []
    for i in np.flatnonzero(valid):
        neigh = indices[i]
        usable = valid[neigh]
        if usable.any():
            scores.append(float(np.mean(values[neigh[usable]] == values[i])))
    return float(np.mean(scores)) if scores else float("nan")


def local_distortion(reference: np.ndarray, target: np.ndarray, reference_indices: np.ndarray, sample_size: int, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    cells = np.repeat(np.arange(reference.shape[0]), reference_indices.shape[1])
    neigh = reference_indices.ravel()
    if len(cells) > sample_size:
        take = rng.choice(len(cells), size=sample_size, replace=False)
        cells, neigh = cells[take], neigh[take]
    d_ref = np.linalg.norm(reference[cells] - reference[neigh], axis=1)
    d_target = np.linalg.norm(target[cells] - target[neigh], axis=1)
    valid = (d_ref > 0) & np.isfinite(d_ref) & np.isfinite(d_target)
    d_ref, d_target = d_ref[valid], d_target[valid]
    scale = np.median(d_target) / max(np.median(d_ref), 1e-12)
    rel = np.abs(d_target - scale * d_ref) / np.maximum(scale * d_ref, 1e-12)
    corr = spearmanr(d_ref, d_target).statistic
    return {"median_relative_local_distance_distortion": float(np.median(rel)), "local_distance_spearman": float(corr), "distance_scale_vs_pca50": float(scale)}


def run_quality() -> dict[str, Any]:
    from sklearn.metrics import silhouette_score
    import umap

    config, p = load_config(), configured_paths()
    validate_environment(config, p)
    ensure_tree(p["output"])
    stage = "03d_representation_quality"
    resources = ResourceLog(stage, config)
    source_sha_before = sha256(p["compact"])
    started = time.perf_counter()
    seed = int(config["random_seed"])
    quality = config["quality"]

    with resources.operation("load_representation_checkpoints"):
        source = ad.read_h5ad(p["compact"], backed="r")
        obs = source.obs.copy()
        source_names = source.obs_names.astype(str).to_numpy()
        source.file.close()
        with np.load(p["output"] / "representations" / "03a_pca_representations.npz") as data:
            reps = {"X_pca20": data["X_pca20"].astype(np.float32), "X_pca30": data["X_pca30"].astype(np.float32), "X_pca50": data["X_pca50"].astype(np.float32)}
            if not np.array_equal(data["obs_names"].astype(str), source_names):
                raise RuntimeError("PCA checkpoint cell order mismatch")
        with np.load(p["output"] / "representations" / "03b_diffusion_representation.npz") as data:
            reps["X_diffmap"] = data["X_diffmap"].astype(np.float32)
            if not np.array_equal(data["obs_names"].astype(str), source_names):
                raise RuntimeError("Diffusion checkpoint cell order mismatch")
        with np.load(p["output"] / "representations" / "03c_scvi_representation.npz") as data:
            reps["X_scvi"] = data["X_scvi"].astype(np.float32)
            if not np.array_equal(data["obs_names"].astype(str), source_names):
                raise RuntimeError("scVI checkpoint cell order mismatch")

    conditions = obs["condition"].astype(str).to_numpy()
    labels = obs["marker_inferred_label"].astype(str).to_numpy()
    clades = obs["souporcell_clade"].astype("string")
    clade_valid = clades.notna().to_numpy()
    clade_values = clades.fillna("unmapped").astype(str).to_numpy()
    k = int(quality["n_neighbors"])
    metrics, rare_rows, knn_indices = [], [], {}

    with resources.operation("quality_metrics_and_knn"):
        for name, embedding in reps.items():
            finite = np.isfinite(embedding).all(axis=1)
            distances, indices = knn(embedding, k)
            knn_indices[name] = indices
            rows = np.arange(len(indices))
            graph = sparse.csr_matrix((np.ones(indices.size, dtype=np.uint8), (np.repeat(rows, k), indices.ravel())), shape=(len(rows), len(rows)))
            graph = graph.maximum(graph.T)
            n_graph_components, graph_labels = connected_components(graph, directed=False)
            label_sil = silhouette_score(embedding, labels, sample_size=min(int(quality["silhouette_sample_size"]), len(labels)), random_state=seed)
            condition_sil = silhouette_score(embedding, conditions, sample_size=min(int(quality["silhouette_sample_size"]), len(labels)), random_state=seed)
            metrics.append({
                "representation": name, "n_cells": embedding.shape[0], "n_dimensions": embedding.shape[1],
                "n_nonfinite_cells": int((~finite).sum()),
                "condition_mixing_fraction": float(1 - same_neighbor_fraction(conditions, indices)),
                "marker_label_neighbor_preservation": same_neighbor_fraction(labels, indices),
                "souporcell_clade_neighbor_preservation_mapped": same_neighbor_fraction(clade_values, indices, clade_valid),
                "silhouette_marker_label": float(label_sil), "silhouette_condition": float(condition_sil),
                "graph_components": int(n_graph_components), "largest_graph_component_fraction": float(np.bincount(graph_labels).max() / len(graph_labels)),
                "neighborhood_overlap_vs_pca50": 1.0 if name == "X_pca50" else float(np.mean(neighbor_jaccard(indices, knn_indices.get("X_pca50", indices)))) if "X_pca50" in knn_indices else float("nan"),
            })
            counts = pd.Series(labels).value_counts()
            for label, count in counts.items():
                if count / len(labels) <= float(quality["rare_state_max_fraction"]):
                    mask = labels == label
                    purity = np.mean([np.mean(labels[indices[i]] == label) for i in np.flatnonzero(mask)])
                    rare_rows.append({"representation": name, "label": label, "n_cells": int(count), "fraction": float(count / len(labels)), "same_label_neighbor_fraction": float(purity)})

        base_indices = knn_indices["X_pca50"]
        for row in metrics:
            name = row["representation"]
            if name != "X_pca50":
                row["neighborhood_overlap_vs_pca50"] = float(np.mean(neighbor_jaccard(knn_indices[name], base_indices)))
            row.update(local_distortion(reps["X_pca50"], reps[name], base_indices, int(quality["distance_pair_sample_size"]), seed))

    metric_frame = pd.DataFrame(metrics)
    rare_frame = pd.DataFrame(rare_rows)
    atomic_csv(metric_frame, p["output"] / "figure_sources" / "03d_representation_quality_metrics.csv")
    atomic_csv(rare_frame, p["output"] / "figure_sources" / "03d_rare_state_sensitivity.csv")

    umaps = {}
    with resources.operation("display_only_umaps_after_all_representations"):
        for name in ["X_pca50", "X_diffmap", "X_scvi"]:
            reducer = umap.UMAP(
                n_neighbors=int(quality["umap_n_neighbors"]), min_dist=float(quality["umap_min_dist"]),
                n_components=2, metric=quality["umap_metric"], random_state=seed, transform_seed=seed, low_memory=True, n_jobs=1,
            )
            umaps[f"X_umap_{name.removeprefix('X_')}"] = reducer.fit_transform(reps[name]).astype(np.float32)

    with resources.operation("diagnostic_figures_and_sources"):
        figure_source = pd.DataFrame({"obs_name": source_names, "condition": conditions, "marker_inferred_label": labels, "souporcell_clade": clade_values})
        for name, coords in umaps.items():
            figure_source[f"{name}_1"] = coords[:, 0]
            figure_source[f"{name}_2"] = coords[:, 1]
        atomic_csv(figure_source, p["output"] / "figure_sources" / "03d_display_umap_coordinates.csv.gz", compression="gzip")
        make_quality_figures(figure_source, metric_frame, p["output"], seed)

    with resources.operation("assemble_lightweight_representation_object"):
        representation = ad.AnnData(obs=obs.copy())
        representation.obs_names = pd.Index(source_names)
        for name, embedding in reps.items():
            representation.obsm[name] = embedding
        for name, coords in umaps.items():
            representation.obsm[name] = coords
        representation.uns["dataset_b_inference_status"] = "descriptive_only"
        representation.uns["source_compact_h5ad"] = str(p["compact"])
        representation.uns["source_compact_sha256"] = source_sha_before
        representation.uns["raw_expression_included"] = False
        representation.uns["umap_role"] = "display_only"
        representation.uns["pca_views_share_one_basis"] = True
        representation.uns["condition_is_biological_not_batch"] = True
        representation.uns["forbidden_methods_run"] = []
        p["representation"].parent.mkdir(parents=True, exist_ok=True)
        tmp = p["representation"].with_name(f".{p['representation'].name}.tmp.{os.getpid()}.h5ad")
        representation.write_h5ad(tmp, compression="gzip")
        backed = ad.read_h5ad(tmp, backed="r")
        try:
            required = {*reps, *umaps}
            if backed.n_obs != len(obs) or not required.issubset(backed.obsm.keys()) or backed.n_vars != 0:
                raise RuntimeError("Representation H5AD backed validation failed")
        finally:
            backed.file.close()
        os.replace(tmp, p["representation"])

    object_sha = sha256(p["representation"])
    summary = {
        "status": "passed", "representation_object": str(p["representation"]),
        "representation_object_shape": [int(representation.n_obs), int(representation.n_vars)],
        "representation_object_sha256": object_sha, "raw_expression_included": False,
        "representations": {name: [int(v) for v in arr.shape] for name, arr in {**reps, **umaps}.items()},
        "runtime_seconds": time.perf_counter() - started,
        "source_sha256_before": source_sha_before, "source_sha256_after": sha256(p["compact"]),
        "inference_status": "descriptive_only", "scgeo_run": False, "comparators_run": False,
    }
    atomic_json(summary, p["output"] / "audit" / "03d_quality_and_object_summary.json")
    resources.write(p["output"])
    write_stage_metadata(stage, config, p, resources, summary)
    del representation, reps, umaps
    close_and_collect()
    return summary


def make_quality_figures(source: pd.DataFrame, metrics: pd.DataFrame, output: Path, seed: int) -> None:
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(seed)
    n = min(12000, len(source))
    take = np.sort(rng.choice(len(source), size=n, replace=False))
    conditions = sorted(source["condition"].unique())
    labels = sorted(source["marker_inferred_label"].unique())
    condition_colors = dict(zip(conditions, plt.cm.Set1(np.linspace(0, 1, len(conditions)))))
    label_colors = dict(zip(labels, plt.cm.tab10(np.linspace(0, 1, len(labels)))))
    bases = ["pca50", "diffmap", "scvi"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for col, base in enumerate(bases):
        x, y = source[f"X_umap_{base}_1"].to_numpy()[take], source[f"X_umap_{base}_2"].to_numpy()[take]
        for value in conditions:
            mask = source["condition"].to_numpy()[take] == value
            axes[0, col].scatter(x[mask], y[mask], s=1, alpha=.45, color=condition_colors[value], rasterized=True, label=value)
        for value in labels:
            mask = source["marker_inferred_label"].to_numpy()[take] == value
            axes[1, col].scatter(x[mask], y[mask], s=1, alpha=.45, color=label_colors[value], rasterized=True, label=value)
        axes[0, col].set_title(f"{base}: condition (descriptive)")
        axes[1, col].set_title(f"{base}: conservative marker label")
        for ax in axes[:, col]:
            ax.set_xticks([]); ax.set_yticks([])
    axes[0, 0].legend(markerscale=5, frameon=False)
    axes[1, 0].legend(markerscale=5, frameon=False, fontsize=7)
    fig.suptitle("GSE249479 display-only UMAP diagnostics; no batch correction")
    fig.tight_layout()
    for ext in ["png", "svg"]:
        fig.savefig(output / "figures" / f"03d_display_umap_diagnostics.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    panels = [
        ("marker_label_neighbor_preservation", "Marker-label neighbour preservation"),
        ("souporcell_clade_neighbor_preservation_mapped", "Mapped-clade neighbour preservation"),
        ("neighborhood_overlap_vs_pca50", "Neighbour overlap vs PCA50"),
    ]
    for ax, (column, title) in zip(axes, panels):
        ax.bar(metrics["representation"], metrics[column], color="#4c78a8")
        ax.set_title(title); ax.tick_params(axis="x", rotation=45); ax.set_ylim(0, 1)
    fig.tight_layout()
    for ext in ["png", "svg"]:
        fig.savefig(output / "figures" / f"03d_representation_quality.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)

    (output / "alt_text" / "03d_display_umap_diagnostics.txt").write_text(
        "Six display-only UMAP panels compare PCA50, diffusion-map and scVI coordinates. The top row colours cells by PBS, TNF or LPS; the bottom row colours conservative marker-based labels. UMAP geometry is not used for inference and no batch correction was applied.\n",
        encoding="utf-8",
    )
    (output / "alt_text" / "03d_representation_quality.txt").write_text(
        "Three bar charts compare five representations on conservative marker-label neighbour preservation, mapped SouporCell-clade neighbour preservation, and neighbourhood overlap relative to PCA50. PCA20, PCA30 and PCA50 are nested views of one shared basis.\n",
        encoding="utf-8",
    )
