"""Synthetic tests of the molecular application contracts, not real-data validation."""
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
from sklearn.decomposition import PCA, TruncatedSVD
import scgeo as sg

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from scgeo_public_molecular import (
    REP, normalized_hvgs, recover_centered_pca, verify_map, MappingError,
    qualified_mask, analyze_state,
)


class Recorder:
    def __init__(self, path):
        self.out = path
        self.names = []

    def save(self, result, name, note):
        try:
            result.save(self.out / name, formats=("png", "svg", "pdf"), dpi=70)
            self.names.append(name)
        finally:
            plt.close(result.fig)


class MolecularContracts(unittest.TestCase):
    def test_distinct_normalization_denominators(self):
        a = ad.AnnData(sparse.csr_matrix([[2., 3., 5.], [1., 1., 8.]]),
                       obs=pd.DataFrame({"total_counts": [10., 10.]}, index=["a", "b"]),
                       var=pd.DataFrame({"highly_variable": [True, True, False]},
                                        index=["g1", "g2", "g3"]))
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "input.h5ad"
            a.write_h5ad(path)
            b, bv = normalized_hvgs(path, "B", a.obs_names, 2)
            c, cv = normalized_hvgs(path, "C", a.obs_names, 2)
            np.testing.assert_allclose(b.toarray(), np.log1p([[2000, 3000], [1000, 1000]]), rtol=1e-6)
            np.testing.assert_allclose(c.toarray(), np.log1p([[4000, 6000], [5000, 5000]]), rtol=1e-6)
            self.assertEqual(list(bv.index), list(cv.index))

    def test_pca_recovery_and_invalid_projection(self):
        rng = np.random.default_rng(17)
        x = sparse.csr_matrix(rng.poisson(2, size=(180, 15)).astype(np.float32))
        model = PCA(n_components=5, svd_solver="arpack", random_state=17)
        z = model.fit_transform(x).astype(np.float32)
        L = recover_centered_pca(x, z)
        audit = verify_map(x, z, L, centered=True)
        self.assertEqual(audit["n_verified_cells"], 180)
        np.testing.assert_allclose(L, model.components_.T, rtol=2e-4, atol=2e-5)
        wrong = z.copy()
        wrong[0, 0] += 0.5
        with self.assertRaises(MappingError):
            verify_map(x, wrong, L, centered=True)
        svd = TruncatedSVD(n_components=5, random_state=17).fit(x)
        verify_map(x, svd.transform(x), svd.components_.T, centered=False)

    def test_real_api_signed_attribution_and_qualified_samples(self):
        rng = np.random.default_rng(21)
        obs = pd.DataFrame([dict(state="A", condition=g, sample=f"{g}_{j}")
                            for g in ("control", "treated") for j in range(3) for _ in range(25)])
        obs.index = pd.Index([f"c{i}" for i in range(len(obs))])
        x = rng.poisson(3, size=(len(obs), 12)).astype(np.float32)
        x[obs.condition.eq("treated"), :3] += 3.
        a = ad.AnnData(sparse.csr_matrix(x), obs=obs)
        a.var_names = [f"g{i}" for i in range(a.n_vars)]
        a.var["original_feature_id"] = a.var_names
        a.var["display_gene_symbol"] = a.var_names
        model = PCA(n_components=5, svd_solver="arpack", random_state=17)
        a.obsm[REP] = model.fit_transform(a.X)
        a.obsm["X_display"] = a.obsm[REP][:,:2].copy()
        a.varm["molecular_loadings"] = recover_centered_pca(a.X, a.obsm[REP])
        verify_map(a.X, a.obsm[REP], a.varm["molecular_loadings"], centered=True)
        row = dict(state="A", contrast="control_vs_treated", hierarchy="major",
                   effect_status="assessed", n_cells0=75, n_cells1=75, n_samples0=3, n_samples1=3,
                   eligible_mouse_ids_group0=json.dumps([f"control_{i}" for i in range(3)]),
                   eligible_mouse_ids_group1=json.dumps([f"treated_{i}" for i in range(3)]))
        frozen = sg.tl.robust_shift(a, rep=REP, condition_key="condition",
            group0="control", group1="treated", sample_key="sample", n_boot=0)
        row["delta_norm"] = frozen["global"]["delta_norm"]
        provenance = {"basis_scope": "verified_frozen_pca50_coordinates"}
        with tempfile.TemporaryDirectory() as d:
            app = Recorder(Path(d))
            result = analyze_state(app, a, row, "state", "condition", "control", "treated",
                {"planted": ["g0", "g1", "g2"]}, provenance, "display", "fixture",
                True, sample="sample")
            self.assertEqual(result["status"], "ok")
            self.assertEqual(result["n_boot"], 200)
            self.assertAlmostEqual(result["contribution_sum"], result["new_delta_norm"], places=4)
            self.assertEqual(len(app.names), 3)
            genes = pd.read_csv(app.out / "fixture_genes.csv.gz")
            self.assertEqual(len(genes), a.n_vars)
            self.assertTrue(genes.ci_low.notna().all())
            self.assertEqual(set(genes.basis_scope), {"verified_frozen_pca50_coordinates"})
            bad = dict(row, n_cells0=74)
            with self.assertRaises(ValueError):
                qualified_mask(a, bad, "state", "condition", "control", "treated", "sample")
            # No biological sample ID is supplied for this second explanation.
            pooled = sg.tl.robust_shift(a, rep=REP, condition_key="condition",
                group0="control", group1="treated", n_boot=0)
            row["delta_norm"] = pooled["global"]["delta_norm"]
            b = analyze_state(app, a, row, "state", "condition", "control", "treated",
                {}, provenance, "display", "pooled", False)
            self.assertEqual(b["status"], "insufficient_samples")
            self.assertEqual(b["n_boot"], 0)
            self.assertTrue(pd.read_csv(app.out / "pooled_genes.csv.gz").ci_low.isna().all())


if __name__ == "__main__":
    unittest.main()
