"""Re-extract the compact public GSE211713 result snapshot from the review archive.

The exact archive is not distributed with the gallery. No archive members are
extracted to disk; only the curated numerical fields are published.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

import numpy as np
import pandas as pd

ARCHIVE_SHA = "5560cc696876e2ab04c6424871dbe8fbf1dc9ac3bbebb5d5bfe22dd46fd92393"
CONTRASTS = ["control_vs_17gy_early", "control_vs_17gy_late", "17gy_early_vs_late"]
GENES = ["Cdkn1a", "Ifi27l2a", "Col3a1", "Spp1", "Igkc", "Sftpc"]
PROGRAMS = ["fibroblast/stromal", "Col13a1_like_matrix_fibroblast",
            "Col14a1_like_matrix_fibroblast", "myofibroblast", "epithelial", "B_cell", "macrophage"]


def r(value):
    return round(float(value), 9)


def extract(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != ARCHIVE_SHA:
        raise ValueError("Archive is not the reviewed result bundle")
    with tarfile.open(path, "r:gz") as t:
        def read(name):
            with t.extractfile(name) as f:
                return pd.read_csv(f, compression="gzip" if name.endswith(".gz") else None)

        s = read("molecular_state_summary.csv")
        states = list(s.target_state.unique())
        assert len(s) == 27 and len(states) == 9 and s.status.eq("ok").all()
        assert s.frozen_magnitude_checked.all() and (s.frozen_delta_norm == s.new_delta_norm).all()
        summary, genes, programs, mice, axes = [], [], [], [], {}
        for name in sorted(n.name for n in t if n.name.endswith("_genes.csv.gz")):
            g = read(name)
            first = g.iloc[0]
            si, ci = states.index(first.target_state), CONTRASTS.index(first.contrast)
            assert len(g) == 3500 and g.status.eq("ok").all()
            assert g.original_feature_id.is_unique
            a = g.set_index("original_feature_id").axis_loading.sort_index()
            axes[si, ci] = a
            row = s[(s.target_state == first.target_state) & (s.contrast == first.contrast)].iloc[0]
            assert abs(g.contribution.sum() - row.delta_norm) < 2e-5
            summary.append([si, ci, *[int(row[k]) for k in ("n_cells0", "n_cells1", "n_samples0", "n_samples1")],
                            *[r(row[k]) for k in ("delta_norm", "contribution_sum", "mapping_max_abs_residual", "mapping_default_state_tolerance")]])
            if first.hierarchy != "fibroblast":
                continue
            for j, symbol in enumerate(GENES):
                if not ((ci == 0 and j < 3) or (ci == 2 and j >= 3)):
                    continue
                hit = g[g.gene_symbol == symbol]
                assert len(hit) == 1
                q = hit.iloc[0]
                genes.append([si, ci, j, *[r(q[k]) for k in ("expression_delta", "axis_loading", "contribution", "ci_low", "ci_high", "sign_stability")]])
            m = read(name.replace("_genes", "_programs")).set_index("program")
            for j, program in enumerate(PROGRAMS):
                q = m.loc[program]
                programs.append([si, ci, j, int(q.n_genes), r(q.contribution)])
            if ci == 2:
                for q in read(name.replace("_genes", "_samples")).itertuples():
                    assert str(q.sample).startswith("GSM") and str(q.sample)[3:].isdigit()
                    mice.append([si, q.condition, str(q.sample), int(q.n_cells), r(q.shift_score)])
        directions = []
        for i in range(len(states)):
            a, b = axes[i, 0], axes[i, 1]
            assert a.index.equals(b.index)
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            directions.append([i, r(np.dot(a, b) / (na * nb)), r(na), r(nb)])
        audit = json.load(t.extractfile("mapping_audit.json"))
    out = {
        "archive_sha256": digest, "states": states, "contrasts": CONTRASTS,
        "gene_names": GENES, "program_names": PROGRAMS,
        "summary_columns": ["state_idx", "contrast_idx", "n_cells0", "n_cells1", "n_samples0", "n_samples1", "delta_norm", "contribution_sum", "mapping_residual", "default_tolerance"],
        "summary": sorted(summary),
        "direction_columns": ["state_idx", "cosine_early_late", "norm_early_loading_axis", "norm_late_loading_axis"],
        "directions": directions,
        "gene_columns": ["state_idx", "contrast_idx", "gene_idx", "expression_delta", "axis_loading", "contribution", "ci_low", "ci_high", "sign_stability"],
        "genes": sorted(genes),
        "program_columns": ["state_idx", "contrast_idx", "program_idx", "n_genes", "contribution"],
        "programs": sorted(programs),
        "mouse_columns": ["state_idx", "condition", "public_accession", "n_cells", "shift_score"],
        "mice": sorted(mice),
        "mapping_audit": {k: v for k, v in audit.items() if k != "intercept"},
    }
    return "{\n" + ",\n".join("  " + json.dumps(k) + ": " + json.dumps(v, separators=(",", ":")) for k, v in out.items()) + "\n}\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--check", type=Path, required=True, help="Compare with the committed snapshot; never overwrite it")
    args = parser.parse_args()
    result = extract(args.archive)
    if result != args.check.read_text():
        raise SystemExit("FAIL: curated snapshot differs")
    print("PASS: exact curated snapshot reproduced from the reviewed archive")
