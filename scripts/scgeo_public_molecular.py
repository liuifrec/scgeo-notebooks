"""Public molecular shift accounting, separate from frozen evidence rendering.

The original HVG normalization differs between B and C. Never use arbitrary
loadings with a frozen representation. Recover C's centered PCA loadings from
the exact score matrix, validate the complete affine map, and fail closed.
An explicit fallback may fit a new, separately labeled linear representation.
"""
from __future__ import annotations

import gc
import hashlib
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.utils.sparsefuncs import inplace_row_scale
import scgeo as sg

from scgeo_v020_applications import ROOT, PRIMARY, OTHER, require, report_for, sha256

RTOL, ATOL = 1e-5, 1e-6
SEED = 1729
REP = "X_molecular_linear50"
FOCUS_B = ("HSC/quiescent", "Activated HSC", "MPP/progenitor")
FOCUS_C = {"major": ("Fibroblast/stromal", "Myeloid"),
           "fibroblast": ("Col13a1-like matrix fibroblast",
                          "Col14a1-like matrix fibroblast", "Myofibroblast")}


class MappingError(ValueError):
    pass


def normalized_hvgs(path, dataset, obs_names, n_hvg, target=10000.0):
    """Read row blocks of original raw counts; never densify the expression."""
    source = ad.read_h5ad(path, backed="r")
    try:
        require(source.obs_names.equals(obs_names), "Compact cell order differs")
        require(source.var_names.is_unique, "Duplicate original feature identifiers")
        selected = np.flatnonzero(source.var["highly_variable"].to_numpy(bool))
        require(len(selected) == n_hvg if dataset == "B" else len(selected) >= n_hvg,
                "Unexpected frozen HVG mask")
        selected = selected[:n_hvg]
        var = source.var.iloc[selected].copy()
        blocks = []
        for start in range(0, source.n_obs, 2048):
            stop = min(start + 2048, source.n_obs)
            raw = source.X[start:stop]
            require(sparse.issparse(raw), "Expected sparse raw counts")
            x = raw[:, selected].tocsr().astype(np.float32, copy=True)
            require(np.isfinite(x.data).all() and (x.data >= 0).all(),
                    "Nonfinite or negative count values")
            if dataset == "B":
                totals = source.obs["total_counts"].iloc[start:stop].to_numpy(np.float32)
                scale = np.divide(np.float32(target), totals, out=np.zeros_like(totals),
                                  where=totals > 0)
                require(np.isfinite(totals).all() and (totals >= 0).all(), "Invalid counts")
                inplace_row_scale(x, scale)
                x.data = np.log1p(x.data).astype(np.float32, copy=False)
            else:
                totals = np.asarray(x.sum(axis=1)).ravel().astype(np.float32)
                scale = np.divide(target, totals, out=np.zeros_like(totals), where=totals > 0)
                x = (sparse.diags(scale.astype(np.float32)) @ x).tocsr()
                np.log1p(x.data, out=x.data)
            blocks.append(x)
        return sparse.vstack(blocks, format="csr"), var
    finally:
        source.file.close()


def verify_map(x, z, loadings, *, centered):
    """Audit every coordinate, with the same maximum-error criterion as ScGeo."""
    z = np.asarray(z, dtype=np.float64)
    L = np.asarray(loadings, dtype=np.float64)
    if (L.shape != (x.shape[1], z.shape[1]) or z.shape[0] != x.shape[0]
            or not np.isfinite(L).all() or not np.isfinite(z).all()):
        raise MappingError("Invalid loading or score dimensions")
    gram_error = float(np.max(np.abs(L.T @ L - np.eye(L.shape[1]))))
    if gram_error > 5e-4:
        raise MappingError(f"Nonorthonormal recovered loadings: {gram_error:g}")
    intercept = np.asarray(z.mean(axis=0) - np.asarray(x.mean(axis=0, dtype=np.float64)).ravel() @ L)
    scale = float(np.max(np.abs(z - z.mean(axis=0))))
    tolerance = ATOL + RTOL * scale
    residual = 0.0
    for start in range(0, len(z), 4096):
        sl = slice(start, min(start + 4096, len(z)))
        err = z[sl] - np.asarray(x[sl] @ L) - intercept
        residual = max(residual, float(np.max(np.abs(err))))
    if not np.isfinite(residual) or residual > tolerance:
        raise MappingError(f"Full frozen affine residual {residual:g} exceeds {tolerance:g}")
    if not centered and np.max(np.abs(intercept)) > tolerance:
        raise MappingError("Uncentered SVD requires a near-zero intercept")
    return {"max_abs_residual": residual, "tolerance": tolerance,
            "max_loading_orthogonality_error": gram_error,
            "intercept": intercept.tolist(), "n_verified_cells": len(z),
            "mapping_rtol": RTOL, "mapping_atol": ATOL}


def recover_centered_pca(x, z):
    """For true centered PCA, Xc.T @ Z = L @ (Z.T @ Z).

    Accumulate in double precision, then validate all original scores.
    This is a recovery attempt, not a new supervised regression or PCA fit.
    """
    z = np.asarray(z, dtype=np.float64)
    zc = z - z.mean(axis=0)
    gram = zc.T @ zc
    require(np.linalg.matrix_rank(gram) == gram.shape[0], "Rank-deficient PCA scores")
    cross = np.zeros((x.shape[1], z.shape[1]), dtype=np.float64)
    for start in range(0, len(z), 4096):
        sl = slice(start, min(start + 4096, len(z)))
        cross += x[sl].astype(np.float64).T @ zc[sl]
    # Xc is implicit because Zc sums to zero.
    return np.linalg.solve(gram, cross.T).T


def refit_linear(x, dataset):
    """Explicit new representation, never claimed to explain frozen scores."""
    if dataset == "B":
        model = TruncatedSVD(n_components=50, algorithm="randomized", n_iter=7,
                             random_state=SEED)
    else:
        model = PCA(n_components=50, svd_solver="arpack", random_state=SEED)
    model.fit(x)
    L = model.components_.T.astype(np.float64)
    intercept = -np.asarray(model.mean_, dtype=float) @ L if dataset == "C" else 0.0
    z = (x @ L + intercept).astype(np.float32)
    return z, L


def load_molecular(app, rep):
    dataset = app.dataset
    config = json.loads((ROOT / "configs" / (
        "gse249479_representations_v1.json" if dataset == "B"
        else "gse211713_representations_v1.json")).read_text())
    n_hvg = config["pca"]["n_hvg"]
    x, var = normalized_hvgs(app.paths["compact"], dataset, rep.obs_names,
                             n_hvg, config["pca"]["normalization_target_sum"])
    original = np.asarray(rep.obsm["X_pca50"], dtype=np.float32)
    hvg_names = var.index.astype(str).to_numpy()
    hvg_hash = hashlib.sha256("\n".join(hvg_names).encode()).hexdigest()
    provenance = {"dataset": dataset, "hvg_count": n_hvg, "hvg_order_sha256": hvg_hash,
                  "normalization": "log1p(10000 * raw HVG counts / " +
                  ("original total_counts)" if dataset == "B" else "HVG row sum)"),
                  "original_basis": "X_pca50"}
    try:
        if dataset == "B":
            base = app.frozen / "gse249479_dataset_b"
            checkpoint = base / "representations/03a_pca_representations.npz"
            report = base / "audit/03a_pca_summary.json"
            if not checkpoint.is_file() or not report.is_file():
                raise MappingError("Original B SVD checkpoint or its audit is absent")
            app.record(report)
            audit = json.loads(report.read_text())
            require(audit["status"] == "passed", "Original SVD stage did not pass")
            expected = app.records[str(app.paths["compact"].resolve())]["sha256"]
            require(audit["source_sha256_before"] == expected == audit["source_sha256_after"],
                    "SVD checkpoint refers to a different input")
            app.record(checkpoint, audit["checkpoint_sha256"])
            with np.load(checkpoint, allow_pickle=False) as saved:
                require(np.array_equal(saved["hvg_names"].astype(str), hvg_names)
                        and np.array_equal(saved["obs_names"].astype(str), rep.obs_names.astype(str))
                        and np.array_equal(saved["X_pca50"].astype(np.float32), original),
                        "SVD checkpoint feature/cell/score identity mismatch")
                L = saved["loadings"].T.astype(np.float64)
            provenance["method"] = "original_uncentered_svd_checkpoint"
        else:
            stage = app.frozen / "gse211713_dataset_c/c6_v2/metadata/pca_stage.json"
            if stage.is_file():
                app.record(stage)
                audit = json.loads(stage.read_text())
                require(audit["status"] == "passed" and audit["hvg_sha256"] == hvg_hash,
                        "Frozen C PCA feature order differs")
            provenance["hvg_stage_audit_present"] = stage.is_file()
            L = recover_centered_pca(x, original)
            provenance["method"] = "recovered_centered_pca_loadings"
        provenance.update(verify_map(x, original, L, centered=dataset == "C"))
        z = original
        provenance["basis_scope"] = "verified_frozen_pca50_coordinates"
    except MappingError as exc:
        if os.environ.get("SCGEO_MOLECULAR_ALLOW_REFIT") != "1":
            raise MappingError(f"{exc}. No attribution produced. To explicitly analyze a NEW "
                               "linear basis instead, set SCGEO_MOLECULAR_ALLOW_REFIT=1.") from exc
        z, L = refit_linear(x, dataset)
        provenance.update(verify_map(x, z, L, centered=dataset == "C"))
        provenance.update(method="new_linear_fit", basis_scope="new_representation_not_frozen",
                          frozen_mapping_failure=str(exc))
    # Keep stable original IDs; make readable labels with a reversible feature map.
    symbols = (var["gene_symbol"].fillna("").astype(str).to_numpy()
               if dataset == "C" else hvg_names.copy())
    symbols = np.array([s or i for s, i in zip(symbols, hvg_names)], dtype=str)
    counts = pd.Series(symbols).value_counts()
    labels = [s if counts[s] == 1 else f"{s} [{i}]" for s, i in zip(symbols, hvg_names)]
    if len(set(labels)) != len(labels):
        labels = hvg_names.tolist()
    var["original_feature_id"] = hvg_names
    var["display_gene_symbol"] = symbols
    var.index = pd.Index(labels)
    work = ad.AnnData(X=x, obs=rep.obs.copy(), var=var)
    work.obsm[REP] = z
    work.varm["molecular_loadings"] = L
    display = "X_umap_pca50" if dataset == "B" else "X_umap_display_only"
    work.obsm[display] = np.array(rep.obsm[display], copy=True)
    work.var[["original_feature_id", "display_gene_symbol"]].to_csv(
        app.out / "molecular_feature_map.csv", index_label="feature")
    np.savez_compressed(app.out / "verified_linear_map.npz",
                       loadings=L, original_feature_ids=hvg_names,
                       feature_labels=np.asarray(labels, dtype=str),
                       intercept=np.asarray(provenance["intercept"]))
    (app.out / "mapping_audit.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Mapping:", provenance["basis_scope"], provenance["method"],
          "residual", provenance["max_abs_residual"], "tolerance", provenance["tolerance"],
          flush=True)
    return work, provenance, display


def programs_for(app, work):
    path = ROOT / "configs" / ("gse249479_signatures_v1.json" if app.dataset == "B"
                              else "gse211713_annotation_markers_v1.json")
    app.record(path)
    config = json.loads(path.read_text())
    specs = config["signatures"] if app.dataset == "B" else {
        **config["major_modules"], **config["refined_modules"]}
    lookup = {}
    for label, symbol in zip(work.var_names, work.var["display_gene_symbol"]):
        lookup.setdefault(symbol, []).append(label)
    modules, rows = {}, []
    for program, spec in specs.items():
        members = []
        for symbol in dict.fromkeys(spec["genes"]):
            hits = lookup.get(symbol, [])
            resolved = hits[0] if len(hits) == 1 else None
            if resolved is not None:
                members.append(resolved)
            rows.append({"program": program, "requested_symbol": symbol,
                         "feature": resolved, "status": "matched" if resolved else
                         "ambiguous" if hits else "absent_from_hvg_basis"})
        # At least two unambiguous measured members, not the complete original program.
        if len(members) >= 2:
            modules[program] = members
    audit = pd.DataFrame(rows)
    audit["program_used"] = audit.program.isin(modules)
    audit.to_csv(app.out / "program_membership_audit.csv", index=False)
    print("HVG-restricted programs:", len(modules), "of", len(specs), flush=True)
    return modules


def choose_features(genes):
    valid = genes.loc[np.isfinite(genes.contribution)]
    positive = valid.loc[valid.contribution.gt(0)].nlargest(1, "contribution")
    negative = valid.loc[valid.contribution.lt(0)].nsmallest(1, "contribution")
    chosen = list(pd.concat([positive, negative]).gene)
    return chosen or valid.assign(magnitude=valid.contribution.abs()).nlargest(2, "magnitude").gene.tolist()


def qualified_mask(work, row, state_key, condition, group0, group1, sample=None):
    base = work.obs[state_key].astype(str).eq(str(row["state"]))
    if sample is None:
        mask = base & work.obs[condition].isin([group0, group1])
    else:
        mask = base & (
            (work.obs[condition].eq(group0) &
             work.obs[sample].astype(str).isin(json.loads(row["eligible_mouse_ids_group0"]))) |
            (work.obs[condition].eq(group1) &
             work.obs[sample].astype(str).isin(json.loads(row["eligible_mouse_ids_group1"]))))
    for group, suffix in ((group0, "0"), (group1, "1")):
        require(int((mask & work.obs[condition].eq(group)).sum()) == int(row[f"n_cells{suffix}"]),
                "Molecular cell selection does not match frozen evidence")
    return mask.to_numpy()


def analyze_state(app, work, row, state_key, condition, group0, group1,
                  modules, provenance, display, name, plot, sample=None):
    state = str(row["state"])
    identity = {"target_state": state, "contrast": str(row["contrast"]),
                "hierarchy": str(row["hierarchy"]), "basis_scope": provenance["basis_scope"]}
    if row["effect_status"] != "assessed":
        return {**identity, "status": "skipped_frozen_coverage", "n_cells0": row["n_cells0"],
                "n_cells1": row["n_cells1"]}
    mask = qualified_mask(work, row, state_key, condition, group0, group1, sample)
    a = work[mask].copy()
    shift = sg.tl.robust_shift(a, rep=REP, condition_key=condition,
        group0=group0, group1=group1, sample_key=sample, center="geometric_median",
        n_boot=0, seed=SEED, store_key="molecular_shift")
    observed = float(shift["global"]["delta_norm"])
    original = float(row["delta_norm"])
    frozen = provenance["basis_scope"] == "verified_frozen_pca50_coordinates"
    if frozen:
        require(np.isclose(observed, original, rtol=1e-4, atol=1e-5),
                f"New robust point estimate differs from frozen displacement: {name}: "
                f"{observed:g} vs {original:g}")
    result = sg.tl.explain_shift(a, shift_key="molecular_shift",
        loadings_key="molecular_loadings", modules=modules or None,
        min_cells=20, min_samples=3 if sample else 2, n_boot=200 if sample else 0,
        seed=SEED, mapping_rtol=RTOL, mapping_atol=ATOL, store_key="molecular_explanation")
    summary = result["summary"].iloc[0].to_dict()
    if sample is not None:
        require(summary["n_samples0"] == int(row["n_samples0"])
                and summary["n_samples1"] == int(row["n_samples1"]), "Qualified mouse counts differ")
    summary = {**summary, **identity, "frozen_delta_norm": original,
               "new_delta_norm": observed, "frozen_magnitude_checked": frozen,
               "mapping_max_abs_residual": result["mapping"]["max_abs_residual"]}
    for key in ("genes", "programs", "samples"):
        table = result[key].copy()
        for col, value in identity.items():
            table[col] = value
        if key == "genes":
            table["original_feature_id"] = table.gene.map(a.var["original_feature_id"])
            table["gene_symbol"] = table.gene.map(a.var["display_gene_symbol"])
        table.to_csv(app.out / f"{name}_{key}.csv.gz", index=False)
    # No p-values: B pooled cells; C conditional sample stability only.
    scope = ("Pooled descriptive attribution; no biological replicates" if sample is None else
             "Qualified cross-sectional mice; conditional sample bootstrap and sample-range curves")
    scope += "; " + provenance["basis_scope"] + "; geometric accounting, not causal regulation"
    if plot and summary["status"] in ("ok", "insufficient_samples"):
        genes = choose_features(result["genes"])
        if genes:
            app.save(sg.pl.shift_curves(a, var_names=genes, state="all",
                key="molecular_explanation", basis=display,
                uncertainty="sample_range" if sample else "none",
                show_samples=sample is not None, max_points=2500, random_state=SEED,
                title=f"{state} | {group0} vs {group1} | Signed gene contributions",
                return_result=True, show=False), f"{name}_gene_curves", scope)
        p = result["programs"]
        finite = p.loc[np.isfinite(p.contribution)]
        if not finite.empty:
            selected = finite.assign(magnitude=finite.contribution.abs()).nlargest(1, "magnitude")
            app.save(sg.pl.shift_curves(a, program_names=selected.program.tolist(), state="all",
                key="molecular_explanation", basis=display,
                uncertainty="sample_range" if sample else "none",
                show_samples=sample is not None, max_points=2500, random_state=SEED,
                title=f"{state} | {group0} vs {group1} | HVG-restricted program",
                return_result=True, show=False), f"{name}_program_curves",
                scope + "; overlapping program sums are not percentages")
        app.save(sg.pl.shift_embedding(a, states=["all"], key="molecular_explanation",
            basis=display, max_points=5000, random_state=SEED, figsize=(7, 6),
            title=f"{state} | {group0} vs {group1} | Schematic center connector",
            return_result=True, show=False), f"{name}_shift_embedding",
            scope + "; connector is not a lineage or RNA velocity")
    print(name, summary["status"], "norm", observed, flush=True)
    return summary


def run_dataset(app, work, provenance, display, modules):
    rows = []
    if app.dataset == "B":
        mapped = work.obs["marker_inferred_label"].astype(str).map(app.config["detailed_state_mapping"])
        require(mapped.notna().all(), "Unmapped HSPC annotation")
        work.obs["application_state"] = mapped
        for treated in ("TNF", "LPS"):
            report = report_for(app.b_evidence(work, treated))
            for i, (_, row) in enumerate(report.iterrows()):
                rows.append(analyze_state(app, work, row, "application_state", "condition",
                    "PBS", treated, modules, provenance, display, f"B_{treated}_{i:02d}",
                    row["state"] in FOCUS_B))
    else:
        dose = pd.to_numeric(work.obs["dose_gy"]).astype(float)
        month = pd.to_numeric(work.obs["month_post_irradiation"]).astype(float)
        work.obs["application_group"] = np.select(
            [dose.eq(0), dose.eq(17) & month.isin([1, 2]), dose.eq(17) & month.isin([4, 5])],
            ["control", "17Gy_early", "17Gy_late"], default="other_dose_time")
        work = work[work.obs["application_group"].ne("other_dose_time")].copy()
        roster = work.obs.groupby("application_group", observed=True)["mouse_id"].nunique()
        require(roster.reindex(["control", "17Gy_early", "17Gy_late"]).tolist() == [5,4,4],
                "Unexpected primary mouse roster")
        for hierarchy, column in (("major", "state_major"), ("fibroblast", "fibroblast_state")):
            approved = app.config["states"][hierarchy]
            values = work.obs[column].astype("string")
            work.obs["application_state"] = values.where(values.isin(approved), OTHER).fillna(OTHER)
            for contrast, group0, group1 in PRIMARY:
                frame, _ = app.c_evidence(work, contrast, group0, group1, hierarchy)
                for i, (_, row) in enumerate(report_for(frame).iterrows()):
                    rows.append(analyze_state(app, work, row, "application_state",
                        "application_group", group0, group1, modules, provenance, display,
                        f"C_{hierarchy}_{contrast}_{i:02d}", row["state"] in FOCUS_C[hierarchy],
                        sample="mouse_id"))
                gc.collect()
    summary = pd.DataFrame(rows)
    require(summary.status.isin(["ok", "insufficient_samples"]).any(), "No supported molecular results")
    summary.to_csv(app.out / "molecular_state_summary.csv", index=False)
    # Hash new numerical outputs and mapping artifacts, independently of plot exports.
    artifacts = {p.name: sha256(p) for p in sorted(app.out.iterdir())
                 if p.name in ("mapping_audit.json", "molecular_feature_map.csv",
                               "verified_linear_map.npz", "program_membership_audit.csv",
                               "molecular_state_summary.csv") or p.name.endswith(".csv.gz")}
    (app.out / "molecular_manifest.json").write_text(json.dumps(
        {"basis": provenance, "scope": app.scope, "state_rows": len(summary),
         "artifact_sha256": artifacts}, indent=2) + "\n")
    return summary
