"""Reproduce the compact HSPC result snapshot from the exact reviewed archive.

No archive member is extracted to disk. Reject other archives, incomplete
rosters, inconsistent accounting, or fabricated sample inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

import numpy as np
import pandas as pd

ARCHIVE_SHA = "0a3ccd640df283558bcba766f735e8304b9a316190955329819c819ac249b227"
STATES = ["HSC/quiescent", "Activated HSC", "MPP/progenitor",
          "Megakaryocyte/erythroid", "Myeloid", "Lymphoid", "Ambiguous HSPC"]
CONTRASTS = ["PBS_vs_TNF", "PBS_vs_LPS"]
GENES = ["AC105402.3", "HLA-DQA1", "STAT1", "SIPA1L1", "XIST", "FOS"]


def r(value):
    return round(float(value), 9)


def require(test, message):
    if not test:
        raise ValueError(message)


def extract(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    require(digest == ARCHIVE_SHA, "Archive is not the reviewed result bundle")
    with tarfile.open(path, "r:gz") as t:
        names = t.getnames()
        require(len(names) == len(set(names)) == 46, "Unexpected archive members")

        def read(name):
            with t.extractfile(name) as f:
                return pd.read_csv(f, compression="gzip" if name.endswith(".gz") else None)

        s = read("molecular_state_summary.csv")
        m = read("program_membership_audit.csv")
        fm = read("molecular_feature_map.csv")
        audit = json.load(t.extractfile("mapping_audit.json"))
        require(len(s) == 14 and set(s.target_state) == set(STATES)
                and set(s.contrast) == set(CONTRASTS)
                and not s.duplicated(["target_state", "contrast"]).any(),
                "Unexpected state/contrast roster")
        require(s.status.eq("insufficient_samples").all()
                and s.n_samples0.eq(0).all() and s.n_samples1.eq(0).all()
                and s.n_boot.eq(0).all() and s.interval_scope.eq("unavailable").all(),
                "Pooled descriptive contract violated")
        require(s.frozen_magnitude_checked.all() and s.frozen_delta_norm.eq(s.new_delta_norm).all(),
                "Frozen magnitudes differ")
        require(audit["basis_scope"] == "verified_frozen_pca50_coordinates"
                and audit["method"] == "original_uncentered_svd_checkpoint"
                and audit["n_verified_cells"] == 34432
                and audit["max_abs_residual"] <= audit["tolerance"],
                "Original SVD map not verified")
        require(len(fm) == 3000 and fm.original_feature_id.is_unique and fm.feature.is_unique
                and hashlib.sha256("\n".join(fm.original_feature_id.astype(str)).encode()).hexdigest()
                == audit["hvg_order_sha256"], "Frozen feature order differs")
        gs, ps = {}, {}
        for name in names:
            if name.endswith("_genes.csv.gz"):
                f = read(name)
                key = (f.target_state.iloc[0], f.contrast.iloc[0])
                require(key not in gs and len(f) == 3000 and f.original_feature_id.is_unique
                        and set(f.original_feature_id) == set(fm.original_feature_id),
                        "Gene table roster or feature set differs")
                require(f.status.eq("insufficient_samples").all()
                        and f[["ci_low", "ci_high", "sign_stability"]].isna().all().all(),
                        "Unexpected gene uncertainty for pooled cells")
                require((f.expression_delta * f.axis_loading - f.contribution).abs().max() < 1e-12,
                        "Gene product conservation failed")
                gs[key] = f.set_index("original_feature_id").sort_index()
            elif name.endswith("_programs.csv.gz"):
                f = read(name)
                key = (f.target_state.iloc[0], f.contrast.iloc[0])
                require(key not in ps and f.program.is_unique
                        and f[["ci_low", "ci_high"]].isna().all().all(),
                        "Unexpected program roster or uncertainty")
                ps[key] = f.set_index("program")
            elif name.endswith("_samples.csv.gz"):
                f = read(name)
                require(len(f) == 2 and f["sample"].eq("__pooled__").all(),
                        "Pooled-center rows must never be interpreted as sample replicates")
        expected = {(state, contrast) for state in STATES for contrast in CONTRASTS}
        require(set(gs) == set(ps) == expected, "Incomplete gene/program contrast roster")
        programs = list(ps[STATES[0], CONTRASTS[0]].index)
        require(len(programs) == 18, "Unexpected program count")
        rows, directions, gene, prog = [], [], [], []
        for i, st in enumerate(STATES):
            x, y = [gs[st, k].axis_loading for k in CONTRASTS]
            require(x.index.equals(y.index), "Feature identities differ between axes")
            directions.append([r(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))),
                               r(np.linalg.norm(x)), r(np.linalg.norm(y))])
            for j, k in enumerate(CONTRASTS):
                q = s[(s.target_state == st) & (s.contrast == k)].iloc[0]
                f, pp = gs[st, k], ps[st, k]
                require(abs(f.contribution.sum() - q.delta_norm) < 1e-6,
                        "Full gene sum does not reconstruct the observed shift")
                require(list(pp.index) == programs, "Program ordering differs")
                for program in programs:
                    matched = m.loc[(m.program == program) & (m.status == "matched"), "feature"]
                    require(len(matched) == pp.loc[program, "n_genes"]
                            and len(set(matched)) == len(matched)
                            and set(matched).issubset(f.index)
                            and abs(pp.loc[program, "contribution"]
                                    - f.loc[matched, "contribution"].sum()) < 1e-12,
                            "Program contribution does not match its measured members")
                rows.append([i, j, int(q.n_cells0), int(q.n_cells1),
                             r(q.delta_norm), r(q.contribution_sum), r(q.mapping_max_abs_residual)])
                prog.append([i, j, *[r(pp.loc[program, "contribution"]) for program in programs]])
                if i < 3:
                    for g, symbol in enumerate(GENES):
                        hit = f.loc[symbol]
                        require(hit.gene_symbol == symbol, "Selected feature identity differs")
                        gene.append([i, j, g, r(hit.expression_delta),
                                     r(hit.axis_loading), r(hit.contribution)])
    data = {
        "archive_sha256": digest, "states": STATES, "contrasts": CONTRASTS,
        "summary_columns": ["state_idx", "contrast_idx", "n_cells0", "n_cells1",
                            "delta_norm", "contribution_sum", "mapping_residual"],
        "summary": rows,
        "direction_columns": ["cosine", "norm_tnf", "norm_lps"], "directions": directions,
        "gene_names": GENES,
        "gene_columns": ["state_idx", "contrast_idx", "gene_idx", "expression_delta",
                         "axis_loading", "contribution"],
        "genes": gene, "program_names": programs,
        "program_sizes": [int(ps[STATES[0], CONTRASTS[0]].loc[q, "n_genes"]) for q in programs],
        "program_requested": [int(m.program.eq(q).sum()) for q in programs],
        "programs": prog,
        "validation": {
            "state_rows": len(s), "gene_rows": sum(len(f) for f in gs.values()),
            "program_rows": sum(len(f) for f in ps.values()),
            "biological_samples": 0, "bootstrap_draws": 0, "status": "insufficient_samples",
            "max_frozen_error": float(abs(s.frozen_delta_norm - s.new_delta_norm).max()),
            "max_sum_error": float(abs(s.delta_norm - s.contribution_sum).max()),
            "max_center_reconstruction_error": float(s.reconstruction_error.max()),
            "matched_memberships": int(m.status.eq("matched").sum()),
            "requested_memberships": len(m),
            "usable_programs": int(m.groupby("program").program_used.first().sum()),
            "total_programs": m.program.nunique(),
        },
        "mapping_audit": {k: v for k, v in audit.items() if k != "intercept"},
    }
    return "{\n" + ",\n".join("  " + json.dumps(k) + ": " + json.dumps(v, separators=(",", ":"))
                            for k, v in data.items()) + "\n}\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--check", required=True, type=Path)
    args = parser.parse_args()
    if extract(args.archive) != args.check.read_text():
        raise SystemExit("FAIL: reviewed HSPC snapshot differs")
    print("PASS: exact HSPC snapshot reproduced")
