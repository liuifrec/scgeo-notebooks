"""General-purpose ScGeo atlas gallery; public inputs are strictly opt-in.

Without file arguments, generates an entirely synthetic *linear* display
projection. This is not a fitted UMAP or an observed RNA-velocity trajectory.

    python examples/geometry_dynamics_atlas_gallery.py --output-dir results/atlas
    python examples/geometry_dynamics_atlas_gallery.py --pbmc /data/pbmc3k_processed.h5ad
    python examples/geometry_dynamics_atlas_gallery.py --pancreas /data/cellrank_pancreas_preprocessed.h5ad

No in-house samples, archived manuscripts, or local study-specific configuration
are referenced, and the script performs no downloads.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad
from scipy import sparse

import scgeo as sg
from plotting_gallery import make_gallery_data

PUBLIC = {
    "pbmc": {
        "sha256": "0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38",
        "label": "louvain",
        "keys": ["X_umap"],
        "description": "Scanpy processed public PBMC3k; supplied UMAP only; no velocity or contrast inferred",
    },
    "pancreas": {
        "sha256": "e8d672d74f5cbb866fe9b5a7acba828a5a0fa9816331b1ce5a1c1e630de01564",
        "label": "clusters",
        "keys": ["X_umap", "velocity_umap"],
        "description": "CellRank E15.5 public pancreas; supplied UMAP and projected velocities, no refit",
    },
}


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_public(dataset, path):
    source = PUBLIC[dataset]
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = _sha256(path)
    if digest != source["sha256"]:
        raise ValueError(f"Public input hash mismatch for {dataset}: {digest}")
    raw = read_h5ad(path, backed="r")
    try:
        adata = AnnData(obs=raw.obs.copy(deep=True))
        for key in source["keys"]:
            adata.obsm[key] = np.asarray(raw.obsm[key]).copy()
        if source["label"] not in adata.obs:
            raise ValueError(f"Expected public annotation {source['label']!r} is absent.")
        return adata, digest
    finally:
        if getattr(raw, "file", None) is not None:
            raw.file.close()


def run_gallery(output_dir, formats=("png", "svg", "pdf"), *, pbmc=None, pancreas=None):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    exports = {}
    adata = make_gallery_data(seed=17)
    states = ["State A", "State B", "State C"]
    adata.obs["state"] = pd.Categorical(adata.obs["state"], categories=states)
    # Demonstrate reuse of Scanpy's standard categorical color convention.
    adata.uns["state_colors"] = ["#3474A6", "#D7834B", "#419C91"]
    adata.uns["paga"] = {
        "groups": "state",
        "connectivities": sparse.csr_matrix([
            [0, .62, .29], [.62, 0, .72], [.29, .72, 0],
        ]),
    }
    sg.tl.robust_shift(
        adata, rep="X_toy", condition_key="condition", group0="control",
        group1="treated", by="state", sample_key="sample",
        center="geometric_median", n_boot=35, seed=17,
    )
    report = sg.get.state_report(adata, node_key="state")
    result = sg.pl.geometry_dynamics_atlas(
        adata, basis="display", color="state", node_key="state",
        condition_key="condition", group0="control", group1="treated",
        sample_key="sample", report=report,
        panels=("states", "velocity", "paga", "evidence"),
        velocity_mode="stream",
        velocity_options={"grid_size": 26, "n_neighbors": 28,
                          "min_neighbors": 6, "min_coherence": .12},
        paga_options={"graph_node_order": states, "min_cells": 10},
        highlight_states=["State A"],
        title="ScGeo geometry–dynamics atlas · simulated input",
        return_result=True, show=False,
    )
    result.metadata["data_provenance"] = {
        "kind": "fully synthetic",
        "seed": 17,
        "embedding": "first two coordinates of generated X_toy, NOT a fitted UMAP",
        "dynamics": "constructed vectors, NOT fitted RNA kinetics",
        "paga": "constructed symmetric graph, NOT inferred topology",
        "effect": "ScGeo robust_shift run on synthetic three-dimensional representation",
        "biological_claim": "none",
    }
    exports["synthetic_full"] = result.save(out / "synthetic_full", formats=formats, dpi=180)
    plt.close(result.fig)

    # Standard AnnData does not need velocity, a PAGA graph, or an analysis run.
    simple = AnnData(obs=adata.obs[["state"]].copy(deep=True))
    simple.obsm["X_display"] = adata.obsm["X_display"].copy()
    result = sg.pl.geometry_dynamics_atlas(
        simple, basis="display", color="state", panels="auto", show=False,
        title="Cell-state atlas from ordinary AnnData",
    )
    assert result.metadata["panels"] == ["states"]
    result.metadata["data_provenance"] = {"kind": "synthetic", "seed": 17,
                                          "available_panels": ["states"]}
    exports["synthetic_minimal"] = result.save(
        out / "synthetic_minimal", formats=formats, dpi=180)
    plt.close(result.fig)

    for name, path in (("pbmc", pbmc), ("pancreas", pancreas)):
        if path is None:
            continue
        real, digest = _read_public(name, path)
        label = PUBLIC[name]["label"]
        result = sg.pl.geometry_dynamics_atlas(
            real, basis="umap", color=label,
            panels="auto", max_points=5000, random_state=17,
            velocity_options={"grid_size": 28, "n_neighbors": 30,
                              "min_neighbors": 6, "min_coherence": .2},
            show=False, title=f"Public {name.upper()} · supplied inputs",
        )
        result.metadata["data_provenance"] = {
            "kind": "public checksum-pinned",
            "source_sha256": digest,
            "source_description": PUBLIC[name]["description"],
            "n_cells": real.n_obs,
            "velocity_projection_refitted": False,
        }
        exports[name] = result.save(out / name, formats=formats, dpi=180)
        plt.close(result.fig)
    summary = {
        "renderer": "scgeo.pl.geometry_dynamics_atlas",
        "no_private_data": True,
        "no_downloads": True,
        "executed": list(exports),
        "output_manifests": {
            k: str(v["metadata"]) for k, v in exports.items()
        },
    }
    (out / "gallery_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/geometry_dynamics_atlas"))
    parser.add_argument("--formats", nargs="+", choices=("png", "svg", "pdf"),
                        default=("png", "svg", "pdf"))
    parser.add_argument("--pbmc", type=Path, help="Existing checksum-matched public Scanpy PBMC3k H5AD")
    parser.add_argument("--pancreas", type=Path, help="Existing checksum-matched public CellRank pancreas H5AD")
    args = parser.parse_args()
    print(json.dumps(run_gallery(args.output_dir, args.formats,
                                 pbmc=args.pbmc, pancreas=args.pancreas), indent=2))


if __name__ == "__main__":
    main()
