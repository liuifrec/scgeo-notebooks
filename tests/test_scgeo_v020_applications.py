"""Synthetic compatibility checks only; no real-data files in CI."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
import scgeo as sg

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from scgeo_v020_applications import (
    REPS, adapt_b, adapt_c, adapt_centers, adjacency, plot_evidence, audit_exports,
)


def fixture():
    rng = np.random.default_rng(81)
    obs = pd.DataFrame([{"application_state": state, "application_group": group,
                         "condition": "PBS" if group == "control" else "TNF",
                         "mouse_id": f"{group}_{i}"}
                        for group in ("control", "17Gy_early") for i in range(3)
                        for state in ("A", "B", "Other retained cells") for _ in range(32)])
    obs.index = pd.Index([f"c{i}" for i in range(len(obs))])
    a = ad.AnnData(obs=obs)
    xy = rng.normal(size=(a.n_obs, 2))
    xy[:, 0] += obs["application_state"].map({"A": -3., "B": 3., "Other retained cells": 0.})
    xy[:, 1] += obs["application_group"].eq("17Gy_early") * .5
    a.obsm["X_display"] = xy
    graph = sparse.diags([np.ones(a.n_obs-1), np.ones(a.n_obs-1)], [-1,1], format="csr")
    return a, graph


class Recorder:
    def __init__(self, out):
        self.out = out
        self.saved = []

    def save(self, result, name, note):
        try:
            result.metadata["fixture_note"] = note
            result.save(self.out / name, formats=("png","svg","pdf"), dpi=70)
            self.saved.append(name)
        finally:
            plt.close(result.fig)


class PublicApplicationTests(unittest.TestCase):
    def test_b_counts_and_renderers(self):
        a, graph = fixture()
        states = ["A","B"]
        raw = pd.DataFrame([dict(state=state, contrast="PBS_vs_TNF", hierarchy="detailed",
                 representation=rep, inference_status="descriptive_only", status="usable",
                 n_cells_PBS=96, n_cells_TNF=96, robust_displacement_magnitude=.5,
                 normalized_effect=.4) for state in states for rep in REPS])
        frame = adapt_b(raw, a, "TNF", states)
        self.assertTrue(frame["descriptive_only"].all())
        self.assertNotIn("n_samples0", frame)
        self.assertNotIn("normalized_magnitude_ci95_low", frame)
        broken = raw.copy()
        broken.loc[0, "n_cells_PBS"] = 95
        with self.assertRaises(ValueError):
            adapt_b(broken, a, "TNF", states)
        names, matrix = adjacency(a, graph)
        np.testing.assert_allclose(matrix, matrix.T)
        self.assertEqual(set(names), {"A","B","Other retained cells"})
        with tempfile.TemporaryDirectory() as temp:
            app = Recorder(Path(temp))
            plot_evidence(app, a, frame, "PBS", "TNF", "condition", "display", "b", names)
            self.assertEqual(len(app.saved), 4)
            self.assertEqual(audit_exports(app.out)["figures"], 12)

    def test_c_mouse_qualification_and_renderers(self):
        a, graph = fixture()
        states = ["A","B"]
        rows, centers = [], []
        for state in states:
            for rep in REPS:
                row = dict(state=state, contrast="control_vs_17gy_early", hierarchy="major",
                     representation=rep, inference_status="replicate_aware_primary",
                     coverage_status="qualified_per_mouse_before_scgeo",
                     delta_norm=.6, normalized_delta_norm=.4, normalized_ci95_low=.2,
                     normalized_ci95_high=.8, raw_ci95_low=.3, raw_ci95_high=1.2)
                for suffix, group in (("0","control"),("1","17Gy_early")):
                    ids = [f"{group}_{i}" for i in range(3)]
                    row[f"eligible_mouse_ids_group{suffix}"] = json.dumps(ids)
                    row[f"excluded_mouse_ids_group{suffix}"] = "[]"
                    row[f"cells_per_mouse_group{suffix}"] = json.dumps(dict.fromkeys(ids,32))
                    row[f"n_mice_group{suffix}"] = 3
                    if rep == "X_pca50":
                        for i, mouse in enumerate(ids):
                            centers.append(dict(state=state, contrast=row["contrast"], hierarchy="major",
                                representation=rep, group=group, mouse_id=mouse, state_cells=32,
                                dimensions=50, coordinate_1=float(i)))
                rows.append(row)
        raw = pd.DataFrame(rows)
        frame = adapt_c(raw, a, "control_vs_17gy_early", "control", "17Gy_early", "major",
                        states, 30, 3)
        self.assertTrue(frame["n_samples0"].eq(3).all())
        self.assertTrue(frame["n_cells0"].eq(96).all())
        tidy = adapt_centers(pd.DataFrame(centers), frame, "control_vs_17gy_early", "major",
                             "control", "17Gy_early")
        self.assertEqual(len(tidy), 12)
        broken = raw.copy()
        broken.loc[0,"excluded_mouse_ids_group0"] = '["unknown"]'
        with self.assertRaises(ValueError):
            adapt_c(broken, a, "control_vs_17gy_early", "control", "17Gy_early",
                    "major", states, 30, 3)
        names, _ = adjacency(a, graph)
        with tempfile.TemporaryDirectory() as temp:
            app = Recorder(Path(temp))
            plot_evidence(app, a, frame, "control", "17Gy_early", "application_group",
                          "display", "c", names, sample="mouse_id")
            app.save(sg.pl.replicate_effect(tidy, coordinate=0, summary=None,
                return_result=True, show=False), "centers", "fixture")
            self.assertEqual(audit_exports(app.out)["figures"], 15)


if __name__ == "__main__":
    unittest.main()
