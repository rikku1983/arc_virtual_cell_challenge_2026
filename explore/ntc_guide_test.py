"""Are the 46 non-targeting guides interchangeable?

For each guide, compare its 400-cell pseudobulk with the pooled other cells, and
compare that distance with random 400-cell groups (the sampling-noise null).
Run: .venv/bin/python explore/ntc_guide_test.py
"""
import gc
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
N_NULL = 200

rows = []
for ctx in "ABC":
    a = ad.read_h5ad(DATA / f"context_{ctx}.h5ad")
    X = sp.csr_matrix(a.X)
    umis = np.asarray(X.sum(1)).ravel()
    L = sp.diags(1e4 / umis) @ X
    L.data = np.log1p(L.data)
    L = L.tocsr()
    expressed = np.asarray(X.mean(0)).ravel() > 0.5  # skip near-zero genes
    L = L[:, expressed]
    total = np.asarray(L.sum(0)).ravel()
    n = L.shape[0]
    groups = a.obs["ntc_id"].astype(str).values

    def dist(idx):
        """RMS log-expression gap between a cell group and all remaining cells."""
        s = np.asarray(L[idx].sum(0)).ravel()
        m_in = s / len(idx)
        m_out = (total - s) / (n - len(idx))
        return float(np.sqrt(np.mean((m_in - m_out) ** 2)))

    null = np.array([dist(rng.choice(n, 400, replace=False)) for _ in range(N_NULL)])
    for g in sorted(set(groups)):
        idx = np.where(groups == g)[0]
        d = dist(idx)
        rows.append({
            "context": ctx, "ntc_id": g, "dist": d,
            "z_vs_random": (d - null.mean()) / null.std(),
            "frac_null_ge": float((null >= d).mean()),
            "median_umis": float(np.median(umis[idx])),
        })
    df = pd.DataFrame([r for r in rows if r["context"] == ctx])
    print(f"\n{ctx}: {expressed.sum()} genes | random-group dist {null.mean():.4f} ± {null.std():.4f} | "
          f"guides {df['dist'].median():.4f} (median), ratio {df['dist'].median() / null.mean():.2f}")
    print(df.sort_values("z_vs_random", ascending=False).head(5)
            [["ntc_id", "dist", "z_vs_random", "median_umis"]].round(4).to_string(index=False))
    del a, X, L
    gc.collect()

pd.DataFrame(rows).to_csv(OUT / "ntc_guide_test.csv", index=False)
