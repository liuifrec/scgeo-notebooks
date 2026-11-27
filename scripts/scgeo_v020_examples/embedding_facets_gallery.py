"""Generate self-contained or public-data examples of ScGeo embedding faceting.

Examples never fit UMAP, RNA velocity, or kinetics. Synthetic vectors are
analytic; public examples use the original supplied X_umap/velocity_umap.
Public files are opt-in and verified against checksums from the ScGeo real-data
plotting gallery. Nothing is downloaded implicitly.

    python examples/embedding_facets_gallery.py --dataset synthetic --output-dir results/facets
    python examples/embedding_facets_gallery.py --dataset pbmc --input-h5ad pbmc3k_processed.h5ad
    python examples/embedding_facets_gallery.py --dataset pancreas --input-h5ad cellrank_pancreas_preprocessed.h5ad
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

import scgeo as sg


PUBLIC = {
    "pbmc": {"sha256": "0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38",
             "n_obs": 2638, "label": "louvain", "size_bytes": 24653425,
             "url": "https://exampledata.scverse.org/scanpy/pbmc3k.h5ad",
             "filename": "pbmc3k_processed.h5ad",
             "source": "https://scanpy.scverse.org/en/stable/generated/scanpy.datasets.pbmc3k_processed.html"},
    "pancreas": {"sha256": "e8d672d74f5cbb866fe9b5a7acba828a5a0fa9816331b1ce5a1c1e630de01564",
                 "n_obs": 2531, "label": "clusters", "size_bytes": 146972595,
                 "url": "https://ndownloader.figshare.com/files/25030028",
                 "filename": "cellrank_pancreas_preprocessed.h5ad",
                 "source": "https://cellrank.readthedocs.io/en/stable/api/_autosummary/datasets/cellrank.datasets.pancreas.html"},
}


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def synthetic(seed=31):
    """Return a constructed display example, NOT a fitted UMAP or trajectory."""
    rng = np.random.default_rng(seed)
    centers = np.array([[-3.2, -0.5], [0.6, 2.2], [3.0, -0.7]], dtype=float)
    names, coords, velocity, batches, conditions = [], [], [], [], []
    shifts = [(0.0, 0.0), (0.32, -0.12), (-0.3, 0.3)]
    phase = [0.15, -0.55, 0.85]
    for condition, delta, angle in zip(("Baseline", "Perturbation", "Recovery"), shifts, phase):
        for state, center in zip(("State A", "State B", "State C"), centers):
            for batch in ("Batch 1", "Batch 2"):
                count = 45
                point = rng.normal(size=(count, 2)) * np.array([0.70, 0.48]) + center + delta
                direction = np.column_stack([np.cos(angle + point[:, 1] * 0.11),
                                             np.sin(angle + point[:, 0] * 0.06)])
                direction *= 0.20
                coords.extend(point)
                velocity.extend(direction + rng.normal(0, 0.018, (count, 2)))
                names.extend([state] * count)
                batches.extend([batch] * count)
                conditions.extend([condition] * count)
    obs = pd.DataFrame({
        "condition": pd.Categorical(conditions, categories=["Baseline", "Perturbation", "Recovery"]),
        "state": pd.Categorical(names, categories=["State A", "State B", "State C"]),
        "batch": batches,
    }, index=[f"simulated_{j:04d}" for j in range(len(names))])
    result = AnnData(obs=obs)
    result.obsm["X_umap"] = np.asarray(coords)
    result.obsm["velocity_umap"] = np.asarray(velocity)
    result.obs["illustrative_score"] = rng.normal(size=result.n_obs)
    result.uns["state_colors"] = ["#168C91", "#E79B5F", "#897DB5"]
    return result


def _download_public(dataset, dest_dir):
    """Explicit, checksum-pinned public download; never fetch implicitly."""
    spec = PUBLIC[dataset]
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / spec["filename"]
    if path.exists():
        if path.stat().st_size == spec["size_bytes"] and _sha256(path) == spec["sha256"]:
            return path
        raise ValueError(f"Existing public input does not match the pinned source: {path}")
    partial = path.with_suffix(path.suffix + ".partial")
    try:
        request = Request(
            spec["url"],
            headers={
                "User-Agent": "Mozilla/5.0 (compatible; ScGeo-public-gallery/0.1; +https://github.com/liuifrec/scgeo)",
                "Accept": "application/octet-stream, */*",
            },
        )
        with urlopen(request, timeout=80) as response, partial.open("wb") as dest:
            while block := response.read(1024 * 1024):
                dest.write(block)
                if dest.tell() > spec["size_bytes"]:
                    raise ValueError("Remote file exceeds pinned source byte count.")
        if partial.stat().st_size != spec["size_bytes"] or _sha256(partial) != spec["sha256"]:
            raise ValueError("Downloaded public dataset failed pinned size/SHA-256 verification.")
        partial.replace(path)
    finally:
        partial.unlink(missing_ok=True)
    return path


def _prepare_public(path, dataset):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"The input dataset does not exist: {path}")
    expected = PUBLIC[dataset]
    digest = _sha256(path)
    if digest != expected["sha256"]:
        raise ValueError(
            f"SHA-256 mismatch for {dataset}: {digest}. Expected {expected['sha256']}; "
            "use the pinned public source, not another dataset under the same filename."
        )
    raw = read_h5ad(path, backed="r")
    try:
        if raw.n_obs != expected["n_obs"]:
            raise ValueError("Unexpected observation count for pinned public dataset.")
        name = expected["label"]
        if name not in raw.obs or "X_umap" not in raw.obsm:
            raise KeyError(f"Missing supplied public annotation {name} or X_umap.")
        selected = pd.Series(np.ones(raw.n_obs, dtype=bool), index=raw.obs_names)
        if dataset == "pancreas":
            # Choose the four largest *observed* annotations for a legible
            # focused panel; this is a display subset, not a lineage claim.
            keep = raw.obs[name].value_counts().head(4).index
            selected = raw.obs[name].isin(keep)
        obs = raw.obs.loc[selected].copy(deep=True)
        view = AnnData(obs=obs)
        view.obsm["X_umap"] = np.asarray(raw.obsm["X_umap"])[selected.to_numpy()].copy()
        if dataset == "pancreas":
            if "velocity_umap" not in raw.obsm:
                raise KeyError("Public pancreas file lacks upstream projected velocity_umap.")
            view.obsm["velocity_umap"] = np.asarray(raw.obsm["velocity_umap"])[selected.to_numpy()].copy()
        if f"{name}_colors" in raw.uns:
            view.uns[f"{name}_colors"] = np.asarray(raw.uns[f"{name}_colors"]).copy()
        return view, {"sha256": digest, "path": str(path), "n_full": raw.n_obs,
                      "n_selected": view.n_obs, "annotation": name, "source": expected["source"],
                      "subset_policy": "top four observed annotation groups by cell count" if dataset == "pancreas" else "all cells"}
    finally:
        if getattr(raw, "file", None) is not None:
            raw.file.close()


def run(output_dir, dataset="synthetic", input_h5ad=None, download_public=False,
        formats=("png", "svg", "pdf")):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    if dataset == "synthetic":
        if input_h5ad is not None or download_public:
            raise ValueError("Synthetic gallery does not use public-data input/download options.")
        adata = synthetic()
        options = [
            ("synthetic_one_way_grid", {"facet": "condition", "color": "state", "velocity": "grid", 
                                       "velocity_options": {"grid_size": 22, "n_neighbors": 24, "min_neighbors": 5},
                                       "context": "all"}),
            ("synthetic_one_way_stream", {"facet": "condition", "color": "state", "velocity": "stream",
                                         "velocity_options": {"grid_size": 26, "n_neighbors": 24, "min_neighbors": 5},
                                         "context": "all"}),
            ("synthetic_two_way", {"facet_row": "batch", "facet_col": "condition", "color": "state", "velocity": "grid",
                                   "velocity_options": {"grid_size": 19, "n_neighbors": 18, "min_neighbors": 4}}),
            ("synthetic_numeric", {"facet": "condition", "color": "illustrative_score"}),
            ("synthetic_no_velocity", {"facet": "condition", "color": "state"}),
        ]
        provenance = {"dataset": "fully_synthetic", "seed": 31,
                      "synthetic_geometry": "constructed 2D display, not fitted UMAP",
                      "synthetic_dynamics": "analytic projected vectors; no RNA velocity inference"}
    else:
        if dataset not in PUBLIC:
            raise ValueError("dataset must be 'synthetic', 'pbmc' or 'pancreas'.")
        if download_public:
            if input_h5ad is not None:
                raise ValueError("Choose --download-public or --input-h5ad, not both.")
            input_h5ad = _download_public(dataset, output_dir / "_public_inputs")
        if input_h5ad is None:
            raise ValueError("Public datasets require --input-h5ad or explicit --download-public.")
        adata, provenance = _prepare_public(input_h5ad, dataset)
        name = PUBLIC[dataset]["label"]
        if dataset == "pbmc":
            options = [("public_pbmc_umap_facets", {
                "facet": name, "color": name, "col_wrap": 4,
                "max_points": 1800, "context": "all",
            })]
        else:
            options = [("public_pancreas_velocity_facets", {
                "facet": name, "color": name, "col_wrap": 4,
                "velocity": "stream", "context": "all",
                "velocity_options": {"grid_size": 25, "n_neighbors": 30,
                                     "min_neighbors": 5, "min_coherence": 0.2},
            })]
    for stem, kwargs in options:
        result = sg.pl.embedding_facets(
            adata, return_result=True, show=False,
            title=("Synthetic example" if dataset == "synthetic" else f"Public {dataset} dataset")
            + " · shared embedding and scales", **kwargs)
        result.metadata["gallery_provenance"] = provenance
        files = result.save(output_dir / stem, formats=formats, dpi=200, bbox_inches="tight")
        plt.close(result.fig)
        results.append({"name": stem, "cells": adata.n_obs,
                        "panels": result.metadata["n_panels"],
                        "supported": pd.to_numeric(result.data["facets"].n_grid_supported,
                                                    errors="coerce").fillna(0).astype(int).tolist(),
                        "figures": {k: str(v) for k, v in files["figures"].items()}})
    summary = {"dataset": dataset, "provenance": provenance, "figures": results}
    (output_dir / "faceting_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=("synthetic", "pbmc", "pancreas"), default="synthetic")
    parser.add_argument("--input-h5ad", type=Path)
    parser.add_argument("--download-public", action="store_true",
                        help="Explicitly download and checksum-verify a public input; never used by default")
    parser.add_argument("--output-dir", type=Path, default=Path("results/faceting_gallery"))
    parser.add_argument("--formats", nargs="+", choices=("png", "svg", "pdf"), default=["png", "svg", "pdf"])
    args = parser.parse_args()
    result = run(args.output_dir, dataset=args.dataset,
                 input_h5ad=args.input_h5ad, download_public=args.download_public,
                 formats=args.formats)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
