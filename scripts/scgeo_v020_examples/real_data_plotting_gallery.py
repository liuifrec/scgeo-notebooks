"""Validate plotting on checksum-pinned, public PBMC and pancreas data.

Preparation and rendering are separate stages. No biological conditions,
sample identifiers, velocities, or time labels are invented. PBMC mapping is
an internal same-dataset check in the supplied PCA representation, not an
independent classification benchmark. Pancreas vectors are supplied upstream.

    python examples/real_data_plotting_gallery.py --dataset all --stage prepare --download
    python examples/real_data_plotting_gallery.py --dataset all --stage render
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import platform
import warnings
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

import scgeo as sg

SOURCES = {
    "pbmc": {
        "name": "Scanpy processed PBMC3k",
        "url": "https://exampledata.scverse.org/scanpy/pbmc3k.h5ad",
        "documentation": "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.datasets.pbmc3k_processed.html",
        "sha256": "0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38",
        "size_bytes": 24653425,
        "filename": "pbmc3k_processed.h5ad",
    },
    "pancreas": {
        "name": "CellRank preprocessed pancreas E15.5",
        "url": "https://ndownloader.figshare.com/files/25030028",
        "documentation": "https://cellrank.readthedocs.io/en/stable/api/_autosummary/datasets/cellrank.datasets.pancreas.html",
        "sha256": "e8d672d74f5cbb866fe9b5a7acba828a5a0fa9816331b1ce5a1c1e630de01564",
        "size_bytes": 146972595,
        "filename": "cellrank_pancreas_preprocessed.h5ad",
    },
}
FLOW_COLUMNS = ["Supplied annotation", "Held-out kNN label"]
DENSITY_GROUPS = ["CD4 T cells", "CD14+ Monocytes", "B cells", "Dendritic cells"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _json(path, payload):
    Path(path).write_text(
        json.dumps(_plain(payload), indent=2) + "\n", encoding="utf-8"
    )


def require_source(dataset, path, download=False):
    """Download only on request; never overwrite an existing source file."""
    source = SOURCES[dataset]
    path = Path(path).expanduser().resolve()
    if not path.exists():
        if not download:
            raise FileNotFoundError(
                f"Missing {source['name']}: {path}. Supply its local H5AD path or use --download."
            )
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name + ".download")
        try:
            with (
                urlopen(source["url"], timeout=30) as response,
                partial.open("wb") as handle,
            ):
                received = 0
                for block in iter(lambda: response.read(1024 * 1024), b""):
                    received += len(block)
                    if received > source["size_bytes"]:
                        raise ValueError(
                            "Public response exceeds the pinned dataset size."
                        )
                    handle.write(block)
            if (
                partial.stat().st_size != source["size_bytes"]
                or sha256(partial) != source["sha256"]
            ):
                raise ValueError(
                    "Downloaded data do not match the pinned size and SHA-256."
                )
            partial.replace(path)
        finally:
            if partial.exists():
                partial.unlink()
    if path.stat().st_size != source["size_bytes"] or sha256(path) != source["sha256"]:
        raise ValueError(
            f"{path} does not match the pinned public source; it was not modified."
        )
    return path


def _read_coordinates(dataset, path):
    """Copy observations and supplied coordinates, without rewriting raw data.

    The derived container has no expression matrix or genes. This also avoids
    carrying legacy H5AD layer placeholders into the new plotting cache.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        original = read_h5ad(path, backed="r")
    label_key = "louvain" if dataset == "pbmc" else "clusters"
    keys = ["X_pca", "X_umap"] if dataset == "pbmc" else ["X_umap", "velocity_umap"]
    try:
        result = AnnData(obs=original.obs.copy(deep=True))
        for key in keys:
            result.obsm[key] = np.asarray(original.obsm[key]).copy()
            if not np.isfinite(result.obsm[key]).all():
                raise ValueError(f"Unexpected nonfinite supplied values in {key}.")
        categories = original.obs[label_key].cat.categories.astype(str).tolist()
        colors = [str(value) for value in original.uns[f"{label_key}_colors"]]
        if len(categories) != len(colors):
            raise ValueError("Supplied category palette does not match its labels.")
        result.uns[f"{label_key}_colors"] = colors
        source_schema = {
            "shape": list(original.shape),
            "obs_columns": list(original.obs.columns),
            "obsm": {key: list(value.shape) for key, value in original.obsm.items()},
            "layers": [None if key is None else str(key) for key in original.layers],
            "sample_identity_fields": [],
            "read_warnings": sorted({str(warning.message) for warning in caught}),
        }
        if dataset == "pancreas":
            source_schema["upstream_velocity_params"] = _plain(
                original.uns.get("velocity_params", {})
            )
            source_schema["upstream_recover_dynamics"] = _plain(
                original.uns.get("recover_dynamics", {})
            )
            source_schema["observed_day_counts"] = {
                str(key): int(value)
                for key, value in original.obs["day"].value_counts().items()
            }
    finally:
        original.file.close()
    return result, source_schema, dict(zip(categories, colors))


def prepare_dataset(
    dataset, path, prepared_dir, *, download=False, seed=17, query_fraction=0.25, k=25
):
    """Create a derived cache and exact observation/split audit tables."""
    path = require_source(dataset, path, download)
    prepared_dir = Path(prepared_dir)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    adata, schema, palette = _read_coordinates(dataset, path)
    metadata = {
        "dataset": dataset,
        "source": {**SOURCES[dataset], "local_path": str(path)},
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "source_schema": schema,
        "palette": palette,
        "derived_container": "Coordinates and observations only; no expression matrix, refitted embedding, synthetic condition, or sample identities.",
        "versions": {
            name: importlib.metadata.version(name)
            for name in [
                "scgeo",
                "anndata",
                "numpy",
                "pandas",
                "scikit-learn",
                "matplotlib",
            ]
        },
    }
    objects = {"full": adata}
    if dataset == "pbmc":
        if not 0 < query_fraction < 1:
            raise ValueError("query_fraction must lie strictly between zero and one.")
        rng = np.random.default_rng(seed)
        query_mask = np.zeros(adata.n_obs, dtype=bool)
        for category in adata.obs.louvain.cat.categories:
            positions = np.flatnonzero(adata.obs.louvain.eq(category).to_numpy())
            n_query = max(
                1,
                min(len(positions) - 1, int(np.ceil(len(positions) * query_fraction))),
            )
            query_mask[rng.permutation(positions)[:n_query]] = True
        adata.obs["reference_query_role"] = pd.Categorical(
            np.where(query_mask, "query", "reference"),
            categories=["reference", "query"],
        )
        reference, query = adata[~query_mask].copy(), adata[query_mask].copy()
        sg.tl.map_knn(
            reference,
            query,
            label_key="louvain",
            rep="X_pca",
            k=k,
            out_label_key="mapped_label",
            out_conf_key="knn_vote_fraction",
            out_ent_key="knn_vote_entropy",
            out_ood_key="knn_distance_outlier",
        )
        categories = adata.obs.louvain.cat.categories
        query.obs[FLOW_COLUMNS[0]] = query.obs.louvain.copy()
        query.obs[FLOW_COLUMNS[1]] = pd.Categorical(
            query.obs.mapped_label, categories=categories
        )
        query.obs["annotation_agreement"] = query.obs.louvain.astype(str).eq(
            query.obs.mapped_label.astype(str)
        )
        objects.update({"reference": reference, "query": query})
        metadata["mapping"] = {
            "seed": int(seed),
            "query_fraction": float(query_fraction),
            "split_rule": "Within supplied annotation categories, ceil(n*query_fraction) query cells, bounded to keep at least one reference and one query cell; one seeded RNG in supplied category order.",
            "reference_cells": reference.n_obs,
            "query_cells": query.n_obs,
            "supplied_representation": "All 50 supplied PCA coordinates; embedding fitted upstream on the full dataset, including query cells.",
            "scope": "Held-out labels within the same dataset and a fixed shared representation; descriptive software validation, not independent classification accuracy or biological ground truth.",
            "map_knn": _plain(query.uns["scgeo"]["map_knn"]),
            "annotation_agreement_count": int(query.obs.annotation_agreement.sum()),
            "annotation_agreement_fraction": float(
                query.obs.annotation_agreement.mean()
            ),
            "vote_fraction_meaning": "Fraction of k reference neighbors voting for the assigned label; not a calibrated probability or biological confidence.",
        }
        pd.crosstab(
            query.obs[FLOW_COLUMNS[0]], query.obs[FLOW_COLUMNS[1]], dropna=False
        ).to_csv(prepared_dir / "pbmc_query_cross_tab.csv")
        adata.obs.groupby(
            ["louvain", "reference_query_role"], observed=True
        ).size().rename("n_cells").to_csv(prepared_dir / "pbmc_split_counts.csv")
    else:
        metadata["velocity_scope"] = (
            "All supplied pancreas cells at the recorded single embryonic time point; original upstream X_umap and velocity_umap, no refit, new projection, inferred pseudotime, or fabricated condition contrast."
        )
    metadata["files"] = {}
    for name, obj in objects.items():
        file = prepared_dir / f"{dataset}_{name}.h5ad"
        obj.write_h5ad(file, compression="gzip")
        metadata["files"][name] = {
            "filename": file.name,
            "sha256": sha256(file),
            "n_cells": obj.n_obs,
        }
        obj.obs.to_csv(
            prepared_dir / f"{dataset}_{name}_observations.csv", index_label="cell"
        )
    _json(prepared_dir / f"{dataset}_manifest.json", metadata)
    print(f"Prepared {dataset}: {adata.n_obs} observed cells; raw SHA-256 unchanged.")
    return metadata


def load_prepared(dataset, prepared_dir):
    prepared_dir = Path(prepared_dir)
    metadata = json.loads((prepared_dir / f"{dataset}_manifest.json").read_text())
    if metadata["source"]["sha256"] != SOURCES[dataset]["sha256"]:
        raise ValueError(
            "Prepared data refer to a different public source fingerprint."
        )
    objects = {}
    for name, details in metadata["files"].items():
        path = prepared_dir / details["filename"]
        if sha256(path) != details["sha256"]:
            raise ValueError(f"Prepared file changed after preparation: {path}")
        objects[name] = read_h5ad(path)
    return objects, metadata


def _combine(fig, results, axes, plot):
    return sg.pl.PlotResult(
        fig=fig,
        axes=axes,
        data={
            f"{name}_{key}": table
            for name, result in results.items()
            for key, table in result.data.items()
        },
        metadata={
            "plot": plot,
            "panels": {name: result.metadata for name, result in results.items()},
        },
        warnings=[
            warning for result in results.values() for warning in result.warnings
        ],
        artists={name: result.artists for name, result in results.items()},
    )


def _save(result, output_dir, stem, metadata, observations, caption, formats):
    result.metadata["public_data_provenance"] = metadata
    result.metadata["gallery_caption"] = caption
    if isinstance(result.data, pd.DataFrame):
        result.data = {"cells": result.data}
    result.data["source_observations"] = observations.copy()
    result.fig.text(
        0.055,
        0.036,
        "PUBLIC DATA · DESCRIPTIVE SOFTWARE VALIDATION",
        fontsize=8,
        weight="bold",
        color="#627584",
    )
    result.fig.text(0.055, 0.012, caption, fontsize=7, color="#627584")
    exported = result.save(
        Path(output_dir) / stem, formats=formats, dpi=220, bbox_inches="tight"
    )
    plt.close(result.fig)
    return exported


def render_pbmc(prepared_dir, output_dir, formats):
    objects, metadata = load_prepared("pbmc", prepared_dir)
    full, query = objects["full"], objects["query"]
    palette = metadata["palette"]
    order = {
        column: [
            label
            for label in full.obs.louvain.cat.categories
            if query.obs[column].eq(label).any()
        ]
        for column in FLOW_COLUMNS
    }
    exports = {}
    flow_kwargs = {
        "columns": FLOW_COLUMNS,
        "min_count": 1,
        "drop_na": False,
        "normalize": False,
        "category_order": order,
        "order": "input",
        "palette": palette,
        "color_by": "target",
        "return_result": True,
        "show": False,
    }
    for kind, plotter, figsize in [
        ("alluvial", sg.pl.state_flow_alluvial, (11, 7.5)),
        ("chord", sg.pl.state_flow_chord, (10, 8.5)),
    ]:
        fig, ax = plt.subplots(figsize=figsize)
        fig.subplots_adjust(left=0.12, right=0.88, bottom=0.15, top=0.87)
        result = plotter(
            query,
            **flow_kwargs,
            ax=ax,
            label_fontsize=8.5,
            title=f"PBMC3k: supplied annotation and held-out mapping ({kind})",
        )
        exports[kind] = _save(
            result,
            output_dir,
            f"real_pbmc_{kind}",
            metadata,
            query.obs,
            f"Same {query.n_obs} query cells on both sides; ribbons show annotation coassignments, not movement or time.",
            formats,
        )

    xy = np.asarray(full.obsm["X_umap"])
    low, high = xy.min(0) - 0.4, xy.max(0) + 0.4
    extent = (float(low[0]), float(high[0]), float(low[1]), float(high[1]))
    fig = plt.figure(figsize=(13.5, 9.4))
    fig.suptitle(
        "PBMC3k: common density intensity and relative shape",
        x=0.48,
        y=0.97,
        fontsize=14,
    )
    rows = {}
    axes_map = {}
    for row, (name, normalization, y) in enumerate(
        [("common", "global", 0.53), ("relative", "per_group", 0.14)]
    ):
        axes = {
            group: fig.add_axes([0.055 + i * 0.215, y, 0.195, 0.31])
            for i, group in enumerate(DENSITY_GROUPS)
        }
        cax = fig.add_axes([0.93, y + 0.02, 0.014, 0.265])
        fig.text(
            0.055,
            y + 0.355,
            "Common scale: transformed cell counts"
            if row == 0
            else "Relative scale: each annotation peak = 1",
            fontsize=10,
            weight="bold",
            color="#344C5C",
        )
        result = sg.pl.embedding_density(
            full,
            groupby="louvain",
            basis="umap",
            groups=DENSITY_GROUPS,
            gridsize=65,
            grid_extent=extent,
            xlim=extent[:2],
            ylim=extent[2:],
            normalize=normalization,
            contour=True,
            smooth_k=3,
            log1p=True,
            cmap="magma",
            transparent_background=False,
            background="white",
            axes=axes,
            cax=cax,
            return_result=True,
            show=False,
        )
        rows[name] = result
        axes_map.update({f"{name}:{group}": ax for group, ax in axes.items()})
        axes_map[f"{name}:colorbar"] = cax
    density = _combine(fig, rows, axes_map, "pbmc_density_scale_comparison")
    exports["density"] = _save(
        density,
        output_dir,
        "real_pbmc_density",
        metadata,
        full.obs,
        "Four supplied annotations; identical UMAP grid and transforms. Relative intensity does not compare population abundance.",
        formats,
    )

    fig, ax = plt.subplots(figsize=(9, 8))
    fig.subplots_adjust(left=0.10, right=0.84, bottom=0.16, top=0.88)
    ax.set_aspect("equal", adjustable="box")
    confidence = sg.pl.score_embedding(
        query,
        score_key="knn_vote_fraction",
        basis="umap",
        norm=Normalize(0, 1),
        ax=ax,
        return_result=True,
        show=False,
        title="PBMC3k: held-out kNN vote fraction",
    )
    confidence.artists["scores"].colorbar.set_label("Neighbor vote fraction")
    exports["vote_fraction"] = _save(
        confidence,
        output_dir,
        "real_pbmc_vote_fraction",
        metadata,
        query.obs,
        "Neighbor vote fraction on supplied UMAP; query labels held out, but the PCA was fitted upstream on the full dataset.",
        formats,
    )
    return exports


def render_pancreas(prepared_dir, output_dir, formats):
    objects, metadata = load_prepared("pancreas", prepared_dir)
    adata = objects["full"]
    xy = np.asarray(adata.obsm["X_umap"])
    low, high = xy.min(0) - 0.4, xy.max(0) + 0.4
    extent = (float(low[0]), float(high[0]), float(low[1]), float(high[1]))
    fig = plt.figure(figsize=(13.5, 7.8))
    fig.suptitle(
        "Pancreas E15.5: supplied projected velocity", x=0.51, y=0.96, fontsize=15
    )
    axes = [
        fig.add_axes([0.055, 0.26, 0.415, 0.61]),
        fig.add_axes([0.535, 0.26, 0.415, 0.61]),
    ]
    panels = {}
    common = {
        "basis": "umap",
        "velocity_basis": "umap",
        "background_key": "clusters",
        "palette": metadata["palette"],
        "grid_extent": extent,
        "grid_size": 32,
        "n_neighbors": 30,
        "min_neighbors": 5,
        "min_coherence": 0.2,
        "xlim": (extent[0] - 0.8, extent[1] + 0.8),
        "ylim": (extent[2] - 0.8, extent[3] + 0.8),
        "background_size": 5,
        "background_alpha": 0.25,
        "color": "#16837D",
        "scale": 7.0,
        "arrow_linewidth": 0.85,
        "arrow_head_size": 2.8,
        "legend": False,
        "return_result": True,
        "show": False,
    }
    for i, mode in enumerate(["grid", "stream"]):
        panels[mode] = sg.pl.velocity_field(
            adata,
            **common,
            mode=mode,
            ax=axes[i],
            title="Supported local grid"
            if mode == "grid"
            else "Supported streamline field",
        )
    axes[0].set_title(
        f"Grid · display multiplier ×{panels['grid'].metadata['scale_resolved']:.2g}"
    )
    handles = [
        Line2D(
            [], [], marker="o", linestyle="none", color=color, markersize=5, label=label
        )
        for label, color in metadata["palette"].items()
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.09),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    result = _combine(
        fig,
        panels,
        {"grid": axes[0], "stream": axes[1]},
        "pancreas_supplied_velocity_comparison",
    )
    return _save(
        result,
        output_dir,
        "real_pancreas_velocity",
        metadata,
        adata.obs,
        "Existing upstream fitted vectors and UMAP; no kinetics or pseudotime refit. Field paths do not establish lineage or causality.",
        formats,
    )


def validation_summary(exports, prepared_dir, output_dir):
    """Audit the actual exported tables and write a portable compact report."""

    def require(condition, message):
        if not condition:
            raise ValueError(f"Real-data export audit failed: {message}")

    def table(export, name):
        return pd.read_csv(export["tables"][name])

    def plot_metadata(export):
        return json.loads(Path(export["metadata"]).read_text())["plot"]

    def file_record(path):
        path = Path(path)
        return {
            "file": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }

    def source_record(metadata):
        source = dict(metadata["source"])
        raw_path = Path(source.pop("local_path"))
        verified = raw_path.is_file()
        if verified:
            require(
                sha256(raw_path) == source["sha256"],
                f"original {metadata['dataset']} file changed",
            )
        schema = {
            key: value
            for key, value in metadata["source_schema"].items()
            if key != "read_warnings"
        }
        return {
            "source": source,
            "raw_hash_reverified_at_render": verified,
            "schema": schema,
            "prepared_files": metadata["files"],
        }

    package_dir = Path(sg.__file__).resolve().parent
    module_names = [
        "_state_flow.py",
        "_state_flow_alluvial.py",
        "_state_flow_chord.py",
        "_density.py",
        "_score.py",
        "_velocity_field.py",
        "_velocity_grid.py",
        "_result.py",
        "_style.py",
        "_utils.py",
    ]
    code_hashes = {
        f"scgeo/pl/{name}": sha256(package_dir / "pl" / name) for name in module_names
    }
    mapping_module = Path(inspect.getsourcefile(sg.tl.map_knn)).resolve()
    code_hashes[str(mapping_module.relative_to(package_dir.parent))] = sha256(
        mapping_module
    )
    code_hashes["examples/real_data_plotting_gallery.py"] = sha256(__file__)
    summary = {
        "schema_version": "scgeo.real_data_plotting_validation.v1",
        "rendered_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Public-data descriptive software validation; same-cell coassignments and supplied embeddings/vectors, not independent annotation accuracy, lineage, causality, or biological inference.",
        "software_versions": {
            "python": platform.python_version(),
            **{
                name: importlib.metadata.version(name)
                for name in [
                    "scgeo",
                    "anndata",
                    "numpy",
                    "pandas",
                    "scipy",
                    "scikit-learn",
                    "matplotlib",
                ]
            },
        },
        "source_code_sha256": code_hashes,
        "datasets": {},
        "exports": {},
    }
    for dataset, figures in exports.items():
        summary["exports"][dataset] = {
            name: {
                "figures": {
                    fmt: file_record(path) for fmt, path in export["figures"].items()
                },
                "tables": {
                    key: file_record(path) for key, path in export["tables"].items()
                },
                "metadata": file_record(export["metadata"]),
            }
            for name, export in figures.items()
        }

    if "pbmc" in exports:
        objects, metadata = load_prepared("pbmc", prepared_dir)
        full, reference, query = (
            objects[key] for key in ["full", "reference", "query"]
        )
        require(
            set(reference.obs_names).isdisjoint(query.obs_names)
            and set(reference.obs_names) | set(query.obs_names) == set(full.obs_names),
            "reference/query identities do not form a disjoint complete partition",
        )
        pbmc = source_record(metadata)
        pbmc["observed_annotation_counts"] = {
            str(label): int(full.obs.louvain.eq(label).sum())
            for label in full.obs.louvain.cat.categories
        }
        pbmc["mapping"] = metadata["mapping"]
        pbmc["mapping"]["annotation_disagreement_count"] = int(
            (~query.obs.annotation_agreement).sum()
        )
        pbmc["split_counts"] = (
            full.obs.groupby(["louvain", "reference_query_role"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .to_dict("records")
        )
        expected_links = {
            tuple(map(str, labels)): int(count)
            for labels, count in query.obs[FLOW_COLUMNS].value_counts().items()
            if count > 0
        }
        flow = {}
        for kind in ["alluvial", "chord"]:
            export = exports["pbmc"][kind]
            links, nodes, cells = (
                table(export, name) for name in ["links", "nodes", "cells"]
            )
            actual_links = {
                (row.source_label, row.target_label): int(row.count)
                for row in links.itertuples()
            }
            require(
                actual_links == expected_links,
                f"{kind} cross-tab differs from query labels",
            )
            require(
                cells.cell.tolist() == query.obs_names.tolist(),
                f"{kind} cell identities differ",
            )
            require(cells.included.all(), f"{kind} unexpectedly omitted query cells")
            stage_mass = nodes.groupby("stage_index")["count"].sum()
            require(
                (stage_mass == query.n_obs).all(), f"{kind} node mass is not conserved"
            )
            require(
                int(links["count"].sum()) == query.n_obs,
                f"{kind} link mass is not conserved",
            )
            require(
                np.isclose(links.fraction.sum(), 1),
                f"{kind} fractions do not sum to one",
            )
            flow[kind] = {
                "query_cells": len(cells),
                "positive_links": len(links),
                "displayed_links": int(links.displayed.sum()),
                "link_count_sum": int(links["count"].sum()),
                "link_fraction_sum": float(links.fraction.sum()),
                "node_count_by_stage": {
                    str(key): int(value) for key, value in stage_mass.items()
                },
                "matches_observed_query_cross_tab": True,
                "all_cell_identities_preserved": True,
            }
        pbmc["flow_mass_reconciliation"] = flow
        pbmc["observed_query_cross_tab"] = [
            {
                "supplied_annotation": labels[0],
                "mapped_label": labels[1],
                "count": count,
            }
            for labels, count in sorted(expected_links.items())
        ]

        density_export = exports["pbmc"]["density"]
        density_meta = plot_metadata(density_export)["panels"]
        grids = {
            name: table(density_export, f"{name}_grid")
            for name in ["common", "relative"]
        }
        require(
            grids["common"][["group", "ix", "iy", "x", "y", "raw_count"]].equals(
                grids["relative"][["group", "ix", "iy", "x", "y", "raw_count"]]
            ),
            "density panels do not share identical raw counts and bin coordinates",
        )
        density = {}
        for name, grid in grids.items():
            observed_counts = {
                label: int(full.obs.louvain.eq(label).sum()) for label in DENSITY_GROUPS
            }
            raw_counts = {
                str(key): int(value)
                for key, value in grid.groupby("group").raw_count.sum().items()
            }
            require(
                raw_counts == observed_counts,
                f"{name} density bin counts do not match source cells",
            )
            require(
                np.allclose(grid.transformed_value, np.log1p(grid.smoothed_count)),
                f"{name} density transform differs from one log1p",
            )
            if name == "common":
                require(
                    np.allclose(grid.plotted_value, grid.transformed_value),
                    "common intensity was renormalized",
                )
            else:
                peaks = grid.groupby("group").transformed_value.transform("max")
                require(
                    np.allclose(grid.plotted_value, grid.transformed_value / peaks),
                    "relative intensity does not match group peaks",
                )
            density[name] = {
                "grid_rows": len(grid),
                "raw_count_sum": int(grid.raw_count.sum()),
                "raw_count_by_annotation": raw_counts,
                "plotted_peak_by_annotation": grid.groupby("group")
                .plotted_value.max()
                .to_dict(),
                "color_limits": density_meta[name]["color_limits"],
                "contour_levels": density_meta[name]["contour_levels"],
                "normalization_rule": density_meta[name]["normalization_rule"],
            }
        density["same_bin_coordinates_and_raw_counts"] = True
        density["grid_size"] = density_meta["common"]["gridsize"]
        density["grid_extent"] = density_meta["common"]["grid_extent"]
        density["smooth_k"] = density_meta["common"]["smooth_k"]
        density["log1p"] = density_meta["common"]["log1p"]
        pbmc["density"] = density
        votes = table(exports["pbmc"]["vote_fraction"], "cells").set_index("cell")
        require(
            votes.index.tolist() == query.obs_names.tolist(),
            "vote plot cell identities differ",
        )
        require(
            np.allclose(votes.value, query.obs.knn_vote_fraction),
            "vote colors differ from stored mapping votes",
        )
        pbmc["vote_fraction"] = {
            "cells": len(votes),
            "minimum": float(votes.value.min()),
            "maximum": float(votes.value.max()),
            "mean": float(votes.value.mean()),
            "color_limits": [0, 1],
            "matches_stored_mapping_values": True,
            "calibrated_probability": False,
        }
        summary["datasets"]["pbmc"] = pbmc

    if "pancreas" in exports:
        objects, metadata = load_prepared("pancreas", prepared_dir)
        full = objects["full"]
        pancreas = source_record(metadata)
        export = exports["pancreas"]["velocity"]
        panels = plot_metadata(export)["panels"]
        grids = {name: table(export, f"{name}_grid") for name in ["grid", "stream"]}
        require(
            grids["grid"]
            .drop(columns="displayed")
            .equals(grids["stream"].drop(columns="displayed")),
            "grid and stream modes use different estimated fields",
        )
        metrics = {}
        for name in ["grid", "stream"]:
            cells, paths, vectors = (
                table(export, f"{name}_{key}") for key in ["cells", "paths", "vectors"]
            )
            require(
                cells.cell.tolist() == full.obs_names.tolist(),
                f"{name} velocity source identities differ",
            )
            require(
                np.allclose(cells[["x", "y"]], full.obsm["X_umap"]),
                f"{name} source UMAP differs",
            )
            require(
                np.allclose(cells[["vx", "vy"]], full.obsm["velocity_umap"]),
                f"{name} supplied velocity differs",
            )
            require(
                cells.eligible.all(),
                f"{name} unexpectedly excludes supplied finite vectors",
            )
            grid = grids[name]
            details = panels[name]
            metrics[name] = {
                key: details[key]
                for key in [
                    "grid_size",
                    "grid_extent",
                    "n_neighbors",
                    "min_neighbors",
                    "radius",
                    "radius_source",
                    "min_coherence",
                    "scale_requested",
                    "scale_resolved",
                    "stream_segments_checked",
                    "stream_segments_removed",
                    "stream_paths_rendered",
                    "grid_arrows_rendered",
                ]
            }
            metrics[name].update(
                {
                    "eligible_cells": int(cells.eligible.sum()),
                    "grid_nodes": len(grid),
                    "supported_nodes": int(grid.supported.sum()),
                    "support_reason_counts": grid.reason.value_counts().to_dict(),
                    "path_vertices": len(paths),
                    "exported_paths": int(paths.path_id.nunique()),
                    "displayed_vectors": int(vectors.displayed.sum()),
                    "supplied_coordinates_and_vectors_preserved": True,
                }
            )
        pancreas["observed_cluster_counts"] = full.obs.clusters.value_counts().to_dict()
        pancreas["same_estimated_field_in_both_modes"] = True
        pancreas["velocity_scope"] = metadata["velocity_scope"]
        pancreas["velocity"] = metrics
        summary["datasets"]["pancreas"] = pancreas

    summary["all_export_audits_passed"] = True
    _json(Path(output_dir) / "validation_summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dataset", choices=["pbmc", "pancreas", "all"], default="all")
    parser.add_argument("--stage", choices=["prepare", "render", "all"], default="all")
    cache = Path.home() / ".cache" / "scgeo"
    parser.add_argument(
        "--pbmc-h5ad", type=Path, default=cache / SOURCES["pbmc"]["filename"]
    )
    parser.add_argument(
        "--pancreas-h5ad", type=Path, default=cache / SOURCES["pancreas"]["filename"]
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Explicitly download missing checksum-pinned public inputs.",
    )
    parser.add_argument(
        "--prepared-dir", type=Path, default=Path("results/real_data_plotting/prepared")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/real_data_plotting")
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=["png", "svg", "pdf"],
        default=["png", "svg", "pdf"],
    )
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--query-fraction", type=float, default=0.25)
    parser.add_argument("--k", type=int, default=25)
    args = parser.parse_args()
    datasets = ["pbmc", "pancreas"] if args.dataset == "all" else [args.dataset]
    if args.stage in {"prepare", "all"}:
        for dataset in datasets:
            prepare_dataset(
                dataset,
                getattr(args, f"{dataset}_h5ad"),
                args.prepared_dir,
                download=args.download,
                seed=args.seed,
                query_fraction=args.query_fraction,
                k=args.k,
            )
    if args.stage in {"render", "all"}:
        if "pbmc" in datasets and (
            not hasattr(sg.pl, "state_flow_chord")
            or "return_result"
            not in inspect.signature(sg.pl.embedding_density).parameters
        ):
            raise RuntimeError(
                "Rendering requires the unreleased flow/density API in the matching source checkout."
            )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        exports = {}
        with sg.pl.style_context("paper"):
            if "pbmc" in datasets:
                exports["pbmc"] = render_pbmc(
                    args.prepared_dir, args.output_dir, tuple(args.formats)
                )
            if "pancreas" in datasets:
                exports["pancreas"] = {
                    "velocity": render_pancreas(
                        args.prepared_dir, args.output_dir, tuple(args.formats)
                    )
                }
        validation_summary(exports, args.prepared_dir, args.output_dir)
        print(
            f"Rendered public-data figures and exact sources to {args.output_dir.resolve()}"
        )


if __name__ == "__main__":
    main()
