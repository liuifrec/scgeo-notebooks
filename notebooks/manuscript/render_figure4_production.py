"""Render the accepted manuscript Figure 4 with a layout-only legend fix.

This script intentionally does not rerun differential-expression or enrichment
analysis. It reads the accepted analysis outputs produced by
``notebooks/manuscript/test_driver_genes.ipynb`` and regenerates only the
Figure 4 graphics. The scientific inputs, genes, classifications, thresholds,
colors, scales, and panel meanings are unchanged.

The production correction is limited to the middle-column Scanpy dotplots
(panels B, E, and H): the source canvas and composite column are widened and
the size/color legends are given more room with smaller legend typography so
that their labels do not overlap.

Run from the repository root, for example::

    python notebooks/manuscript/render_figure4_production.py

Optional::

    python notebooks/manuscript/render_figure4_production.py \
        --data-dir data --output-dir production/figure4
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.gridspec import GridSpec


DISCORDANT_SET = [
    "stress-associated / ambiguous",
    "stressed erythroid-like",
    "stressed myeloid-like",
    "macrophage-like / remodeling myeloid",
    "divergent erythroid-like",
    "activated APC-like",
    "transitional erythroid",
]

CANONICAL_SET = [
    "HSC / primitive progenitor-like",
    "granulocytic / emergency myeloid",
    "B-lineage / pre-B-like",
    "T/NK-like lymphoid",
]

ERYTHROID_DISCORDANT = [
    "stress-associated / ambiguous",
    "stressed erythroid-like",
    "divergent erythroid-like",
    "transitional erythroid",
]

ERYTHROID_CANONICAL = [
    "erythroid progenitor-like",
    "stable erythroid",
]

MYELOID_DISCORDANT = [
    "stressed myeloid-like",
    "macrophage-like / remodeling myeloid",
    "activated APC-like",
]

MYELOID_CANONICAL = [
    "granulocytic / emergency myeloid",
    "stable APC-like",
]


PALETTE_GLOBAL = {
    "canonical_set": "#4C78A8",
    "discordant_set": "#F58518",
    "other": "#D9D9D9",
}

PALETTE_ERY = {
    "ery_canonical": "#4C78A8",
    "ery_discordant": "#E45756",
    "other": "#D9D9D9",
}

PALETTE_MY = {
    "my_canonical": "#2E8B57",
    "my_discordant": "#B279A2",
    "other": "#D9D9D9",
}


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=repo_root / "data",
        help="Directory containing the accepted GSE280305 analysis outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "production" / "figure4",
        help="Directory for production Figure 4 outputs.",
    )
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(
            f"Required accepted-analysis file not found: {path}\n"
            "Do not regenerate upstream analysis for this production fix. "
            "Use the workstation copy that produced the accepted Figure 4."
        )
    return path


def filter_deg(df: pd.DataFrame, lfc: float = 0.25, padj: float = 0.05) -> pd.DataFrame:
    return df.query("pvals_adj < @padj and logfoldchanges > @lfc").copy()


def top_genes(df: pd.DataFrame, n: int = 10) -> list[str]:
    return df["names"].dropna().head(n).tolist()


def shorten_term(term: object, max_len: int = 24) -> str:
    term = str(term)
    return term if len(term) <= max_len else term[: max_len - 3] + "..."


def assign_accepted_state_sets(adata) -> None:
    labels = adata.obs["cluster_label_manual"]

    adata.obs["state_set"] = "other"
    adata.obs.loc[labels.isin(DISCORDANT_SET), "state_set"] = "discordant_set"
    adata.obs.loc[labels.isin(CANONICAL_SET), "state_set"] = "canonical_set"

    adata.obs["ery_state_set"] = "other"
    adata.obs.loc[labels.isin(ERYTHROID_DISCORDANT), "ery_state_set"] = "ery_discordant"
    adata.obs.loc[labels.isin(ERYTHROID_CANONICAL), "ery_state_set"] = "ery_canonical"

    adata.obs["my_state_set"] = "other"
    adata.obs.loc[labels.isin(MYELOID_DISCORDANT), "my_state_set"] = "my_discordant"
    adata.obs.loc[labels.isin(MYELOID_CANONICAL), "my_state_set"] = "my_canonical"


def plot_set_umap_manual(
    adata,
    key: str,
    ax,
    colors: dict[str, str],
    basis: str = "umap_scanorama",
    title: str | None = None,
    other_label: str = "other",
    size: float = 8,
    alpha_other: float = 0.08,
    alpha_main: float = 0.90,
) -> None:
    key_basis = basis if basis.startswith("X_") else f"X_{basis}"
    xy = adata.obsm[key_basis]
    s = adata.obs[key].astype(str)

    if other_label in set(s):
        mask_other = (s == other_label).to_numpy()
        ax.scatter(
            xy[mask_other, 0],
            xy[mask_other, 1],
            s=size,
            c=colors.get(other_label, "#D9D9D9"),
            alpha=alpha_other,
            linewidths=0,
            rasterized=True,
        )

    for cat, col in colors.items():
        if cat == other_label:
            continue
        mask = (s == cat).to_numpy()
        if mask.sum() == 0:
            continue
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=size,
            c=col,
            alpha=alpha_main,
            linewidths=0,
            rasterized=True,
            label=cat,
        )

    ax.set_title(title, fontsize=14, loc="left", pad=4)
    ax.set_xlabel("UMAP_SCANORAMA1", fontsize=11)
    ax.set_ylabel("UMAP_SCANORAMA2", fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="center left", bbox_to_anchor=(1.00, 0.5))


def plot_enrichment_bar(
    enr_df: pd.DataFrame,
    ax,
    title: str,
    color: str = "#4C78A8",
    top_n: int = 5,
    term_max_len: int = 30,
) -> None:
    df = enr_df.sort_values("Adjusted P-value").head(top_n).copy()
    df = df[df["Adjusted P-value"] > 0].copy()
    df["score"] = -np.log10(df["Adjusted P-value"])
    df["Term_short"] = df["Term"].map(lambda x: shorten_term(x, max_len=term_max_len))
    df = df.sort_values("score", ascending=True)

    ax.barh(df["Term_short"], df["score"], color=color)
    ax.set_xlabel("-log10(adj p)", fontsize=15)
    ax.set_title(title, fontsize=18, loc="left", pad=6)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="x", labelsize=14)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def plot_deg_volcano(
    df: pd.DataFrame,
    ax,
    title: str,
    gene_col: str = "names",
    lfc_col: str = "logfoldchanges",
    p_col: str = "pvals_adj",
    top_n: int = 5,
    lfc_thr: float = 0.25,
    p_thr: float = 0.05,
    up_color: str = "#D55E00",
    ns_color: str = "#BDBDBD",
    x_clip: tuple[float, float] = (-0.5, 6),
    ylim: tuple[float, float] | None = (0, 320),
) -> None:
    d = df.copy()
    d[lfc_col] = pd.to_numeric(d[lfc_col], errors="coerce")
    d[p_col] = pd.to_numeric(d[p_col], errors="coerce")
    d["neglog10p"] = -np.log10(np.clip(d[p_col], 1e-300, 1.0))
    d = d[d[lfc_col] > x_clip[0]].copy()
    d["lfc_plot"] = d[lfc_col].clip(lower=x_clip[0], upper=x_clip[1])
    sig = (d[p_col] < p_thr) & (d[lfc_col] > lfc_thr)

    ax.scatter(
        d.loc[~sig, "lfc_plot"],
        d.loc[~sig, "neglog10p"],
        s=8,
        c=ns_color,
        alpha=0.45,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        d.loc[sig, "lfc_plot"],
        d.loc[sig, "neglog10p"],
        s=10,
        c=up_color,
        alpha=0.80,
        linewidths=0,
        rasterized=True,
    )
    ax.axvline(lfc_thr, linestyle="--", linewidth=1, color="black", alpha=0.6)
    ax.axhline(-np.log10(p_thr), linestyle="--", linewidth=1, color="black", alpha=0.6)

    top = d.loc[sig].sort_values([p_col, lfc_col], ascending=[True, False]).head(top_n)
    offsets = [(4, 8), (4, -10), (6, 14), (6, -16), (8, 20), (8, -22)]
    for i, (_, row) in enumerate(top.iterrows()):
        dx, dy = offsets[i % len(offsets)]
        ax.annotate(
            str(row[gene_col]),
            xy=(row["lfc_plot"], row["neglog10p"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=11,
            ha="left",
            va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="black", alpha=0.5),
        )

    ax.set_xlabel("logFC (clipped for display)", fontsize=12)
    ax.set_ylabel("-log10(adj p)", fontsize=12)
    ax.set_title(title, fontsize=14, loc="left", pad=4)
    ax.set_xlim(*x_clip)
    if ylim is not None:
        ax.set_ylim(*ylim)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def crop_white_rgba(
    img,
    white_thresh: int = 250,
    pad_left: int = 45,
    pad_right: int = 8,
    pad_top: int = 8,
    pad_bottom: int = 18,
):
    rgb = img[..., :3]
    mask = np.any(rgb < white_thresh, axis=2)
    if not mask.any():
        return img

    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    r0 = max(rows[0] - pad_top, 0)
    r1 = min(rows[-1] + pad_bottom + 1, img.shape[0])
    c0 = max(cols[0] - pad_left, 0)
    c1 = min(cols[-1] + pad_right + 1, img.shape[1])
    return img[r0:r1, c0:c1, :]


def draw_scanpy_dotplot_on_ax(
    adata,
    mask,
    var_names: list[str],
    groupby: str,
    ax,
    layer: str = "matrix",
    title: str | None = None,
    figsize: tuple[float, float] = (5.2, 8.5),
) -> None:
    """Draw the accepted dotplot with production-only legend spacing fixes."""
    dp = sc.pl.dotplot(
        adata[mask],
        var_names=var_names,
        groupby=groupby,
        layer=layer,
        swap_axes=True,
        show=False,
        return_fig=True,
        figsize=figsize,
    )

    # Same legend words as the accepted plot, but with explicit line breaks and
    # a wider legend column to prevent the size and color legends from colliding.
    try:
        dp.legend(
            show=True,
            show_size_legend=True,
            size_title="Fraction of cells\nin group (%)",
            show_colorbar=True,
            colorbar_title="Mean expression\nin group",
            width=2.6,
        )
    except TypeError:
        # Compatibility fallback for older Scanpy versions.
        try:
            dp.legend(width=2.6)
        except TypeError:
            pass

    dp.make_figure()

    # Preserve the accepted large main-panel typography, while keeping only the
    # legend typography compact enough to remain non-overlapping.
    for txt in dp.fig.findobj(match=plt.Text):
        txt.set_fontsize(15)

    ax_dict = getattr(dp, "ax_dict", {}) or {}
    for legend_key in ("size_legend_ax", "color_legend_ax"):
        legend_ax = ax_dict.get(legend_key)
        if legend_ax is None:
            continue
        for txt in legend_ax.findobj(match=plt.Text):
            txt.set_fontsize(9)

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        dp.fig.savefig(
            tmp_path,
            dpi=320,
            bbox_inches="tight",
            pad_inches=0.08,
            facecolor="white",
        )
        img = mpimg.imread(tmp_path)
    finally:
        plt.close(dp.fig)
        try:
            os.unlink(tmp_path)
        except FileNotFoundError:
            pass

    img = crop_white_rgba(
        img,
        white_thresh=250,
        pad_left=6,
        pad_right=12,
        pad_top=8,
        pad_bottom=12,
    )

    ax.imshow(img, aspect="auto")
    ax.axis("off")
    ax.set_title(title, fontsize=18, loc="left", pad=6)


def load_accepted_inputs(data_dir: Path):
    adata_path = require_file(data_dir / "scgeo_gse280305_final_integrated.h5ad")
    deg_global_path = require_file(data_dir / "deg_state_set_global.csv")
    deg_ery_path = require_file(data_dir / "deg_state_set_ery.csv")
    deg_my_path = require_file(data_dir / "deg_state_set_my.csv")
    enrichment_path = require_file(data_dir / "enrichment_hallmark_state_sets.csv")

    adata = sc.read_h5ad(adata_path)
    deg_global = pd.read_csv(deg_global_path)
    deg_ery = pd.read_csv(deg_ery_path)
    deg_my = pd.read_csv(deg_my_path)
    enrichment = pd.read_csv(enrichment_path)

    required_obs = {"cluster_label_manual"}
    missing_obs = required_obs.difference(adata.obs.columns)
    if missing_obs:
        raise KeyError(f"Missing required AnnData obs columns: {sorted(missing_obs)}")
    if "X_umap_scanorama" not in adata.obsm:
        raise KeyError("AnnData is missing X_umap_scanorama used by the accepted Figure 4.")
    if "matrix" not in adata.layers:
        raise KeyError("AnnData is missing the accepted 'matrix' layer used by the dotplots.")

    assign_accepted_state_sets(adata)
    return adata, deg_global, deg_ery, deg_my, enrichment


def select_enrichment(enrichment: pd.DataFrame, contrast: str) -> pd.DataFrame:
    out = enrichment.loc[enrichment["contrast"] == contrast].copy()
    if out.empty:
        raise ValueError(f"No enrichment rows found for accepted contrast: {contrast}")
    return out


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    adata, deg_global, deg_ery, deg_my, enrichment = load_accepted_inputs(data_dir)

    deg_global_sig = filter_deg(deg_global)
    deg_ery_sig = filter_deg(deg_ery)
    deg_my_sig = filter_deg(deg_my)
    top_global = top_genes(deg_global_sig, n=10)
    top_ery = top_genes(deg_ery_sig, n=10)
    top_my = top_genes(deg_my_sig, n=10)

    if not (len(top_global) == len(top_ery) == len(top_my) == 10):
        raise ValueError(
            "Accepted Figure 4 expected 10 genes in each dotplot. "
            f"Observed global={len(top_global)}, erythroid={len(top_ery)}, myeloid={len(top_my)}."
        )

    enr_global = select_enrichment(enrichment, "global_discordant_vs_canonical")
    enr_ery = select_enrichment(enrichment, "ery_discordant_vs_canonical")
    enr_my = select_enrichment(enrichment, "my_discordant_vs_canonical")

    print("Accepted top genes retained")
    print("Global:", top_global)
    print("Erythroid:", top_ery)
    print("Myeloid/APC:", top_my)

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 18,
            "axes.labelsize": 14,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
        }
    )

    # Production-only layout correction: slightly wider overall figure and
    # middle column. No scientific inputs or plotting scales are changed.
    fig = plt.figure(figsize=(22, 18), constrained_layout=True)
    gs = GridSpec(
        nrows=3,
        ncols=3,
        figure=fig,
        height_ratios=[1.0, 1.0, 1.0],
        width_ratios=[1.05, 0.90, 1.15],
    )

    ax_a = fig.add_subplot(gs[0, 0])
    plot_deg_volcano(
        deg_global,
        ax=ax_a,
        title="A  Global discordant vs canonical",
        top_n=5,
        x_clip=(-0.5, 6),
        ylim=(0, 320),
    )

    ax_b = fig.add_subplot(gs[0, 1])
    mask_global = adata.obs["state_set"].isin(["canonical_set", "discordant_set"]).to_numpy()
    draw_scanpy_dotplot_on_ax(
        adata,
        mask=mask_global,
        var_names=top_global,
        groupby="state_set",
        ax=ax_b,
        layer="matrix",
        title="B  Top discordant genes",
    )

    ax_c = fig.add_subplot(gs[0, 2])
    plot_enrichment_bar(
        enr_global,
        ax=ax_c,
        title="C  Global enrichment",
        color="#4C78A8",
        top_n=6,
    )

    ax_d = fig.add_subplot(gs[1, 0])
    plot_set_umap_manual(
        adata,
        key="ery_state_set",
        ax=ax_d,
        colors=PALETTE_ERY,
        basis="umap_scanorama",
        title="D  Erythroid discordant vs canonical",
    )

    ax_e = fig.add_subplot(gs[1, 1])
    mask_ery = adata.obs["ery_state_set"].isin(["ery_canonical", "ery_discordant"]).to_numpy()
    draw_scanpy_dotplot_on_ax(
        adata,
        mask=mask_ery,
        var_names=top_ery,
        groupby="ery_state_set",
        ax=ax_e,
        layer="matrix",
        title="E  Erythroid discordant genes",
    )

    ax_f = fig.add_subplot(gs[1, 2])
    plot_enrichment_bar(
        enr_ery,
        ax=ax_f,
        title="F  Erythroid enrichment",
        color="#E45756",
        top_n=6,
    )

    ax_g = fig.add_subplot(gs[2, 0])
    plot_set_umap_manual(
        adata,
        key="my_state_set",
        ax=ax_g,
        colors=PALETTE_MY,
        basis="umap_scanorama",
        title="G  Myeloid/APC discordant vs canonical",
    )

    ax_h = fig.add_subplot(gs[2, 1])
    mask_my = adata.obs["my_state_set"].isin(["my_canonical", "my_discordant"]).to_numpy()
    draw_scanpy_dotplot_on_ax(
        adata,
        mask=mask_my,
        var_names=top_my,
        groupby="my_state_set",
        ax=ax_h,
        layer="matrix",
        title="H  Myeloid/APC discordant genes",
    )

    ax_i = fig.add_subplot(gs[2, 2])
    plot_enrichment_bar(
        enr_my,
        ax=ax_i,
        title="I  Myeloid/APC enrichment",
        color="#B279A2",
        top_n=6,
    )

    tif_path = output_dir / "Figure4_production.tiff"
    pdf_path = output_dir / "Figure4_production.pdf"
    png_path = output_dir / "Figure4_production_preview.png"

    fig.savefig(tif_path, dpi=args.dpi, format="tiff", bbox_inches="tight")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, dpi=180, format="png", bbox_inches="tight")
    plt.close(fig)

    print("Production-only Figure 4 render complete")
    print("TIFF:", tif_path)
    print("PDF:", pdf_path)
    print("Preview:", png_path)
    print("No differential-expression or enrichment analysis was recomputed.")


if __name__ == "__main__":
    main()
