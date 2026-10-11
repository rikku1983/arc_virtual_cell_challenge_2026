"""First look at the 2026 validation control data (contexts A, B, C).

Loads one context at a time (memory is tight) and writes per-context summaries
to explore/out/. Run: .venv/bin/python explore/explore_controls.py
"""
import gc
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

genes = pd.read_csv(DATA / "gene_names.csv")
perts = pd.read_csv(DATA / "pert_counts.csv")
print("gene_names.csv:", genes.shape, list(genes.columns))
print("pert_counts.csv:", perts.shape, list(perts.columns))
gene_list = genes.iloc[:, 0].astype(str).tolist()
pert_list = perts["target_gene"].astype(str).tolist()
print("perts missing from gene list:", sorted(set(pert_list) - set(gene_list)))

summary = {}
mean_expr = {}  # per-context mean CPM per gene, for cross-context comparison
guide_frac = {}  # fraction of cells expressing each target gene

for ctx in "ABC":
    a = ad.read_h5ad(DATA / f"context_{ctx}.h5ad")
    X = a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)
    X = X.tocsr()
    assert list(a.var_names) == gene_list, f"{ctx}: var order differs from gene_names.csv"

    umis = np.asarray(X.sum(1)).ravel()
    ngenes = np.diff(X.indptr)
    is_int = bool(np.all(np.mod(X.data, 1) == 0))

    cpm = sp.diags(1e4 / umis) @ X  # counts per 10k
    mean_expr[ctx] = np.asarray(cpm.mean(0)).ravel()
    gene_det = np.asarray((X > 0).mean(0)).ravel()
    guide_frac[ctx] = pd.Series(gene_det, index=gene_list)

    # guide-to-guide variation: pseudobulk per ntc_id, correlation between guides
    log_cpm = cpm.copy()
    log_cpm.data = np.log1p(log_cpm.data)
    groups = a.obs["ntc_id"].astype(str).values
    ids = sorted(set(groups))
    pb = np.vstack([np.asarray(log_cpm[groups == g].mean(0)).ravel() for g in ids])
    expressed = mean_expr[ctx] > 1.0  # >1 per 10k
    corr = np.corrcoef(pb[:, expressed])
    off = corr[~np.eye(len(ids), dtype=bool)]
    guide_umis = pd.Series(umis).groupby(groups).median()

    summary[ctx] = {
        "shape": list(X.shape),
        "dtype": str(X.dtype),
        "integer_counts": is_int,
        "obs_columns": list(a.obs.columns),
        "n_ntc_ids": len(ids),
        "cells_per_ntc_min_max": [int(pd.Series(groups).value_counts().min()),
                                  int(pd.Series(groups).value_counts().max())],
        "umis_median": float(np.median(umis)),
        "umis_p5_p95": [float(np.percentile(umis, 5)), float(np.percentile(umis, 95))],
        "genes_detected_median": float(np.median(ngenes)),
        "genes_detected_max": int(ngenes.max()),
        "genes_detected_in_any_cell": int((gene_det > 0).sum()),
        "genes_mean_gt_1_per_10k": int(expressed.sum()),
        "guide_pseudobulk_corr_min_median": [float(off.min()), float(np.median(off))],
        "guide_median_umis_min_max": [float(guide_umis.min()), float(guide_umis.max())],
        "pert_targets_detected_in_<10%_cells": int((guide_frac[ctx][pert_list] < 0.1).sum()),
    }
    print(ctx, json.dumps(summary[ctx], indent=1))
    del a, X, cpm, log_cpm
    gc.collect()

# cross-context comparison on log mean expression
lm = pd.DataFrame({c: np.log1p(v) for c, v in mean_expr.items()}, index=gene_list)
print("cross-context correlation of log mean expression:\n", lm.corr().round(4))

# genes most specific to each context (log fold vs mean of the other two)
spec = {}
for c in "ABC":
    others = [o for o in "ABC" if o != c]
    d = lm[c] - lm[others].mean(1)
    spec[c] = {"up": d.nlargest(15).index.tolist(), "down": d.nsmallest(15).index.tolist()}
print("context-specific genes:", json.dumps(spec, indent=1))

tg = pd.DataFrame({c: guide_frac[c][pert_list] for c in "ABC"})
tg.to_csv(OUT / "target_gene_detection_frac.csv")
print("target genes: fraction of control cells with >0 counts, quantiles:\n",
      tg.describe().round(3))

json.dump({"summary": summary, "context_specific": spec}, open(OUT / "summary.json", "w"), indent=1)
lm.to_csv(OUT / "log_mean_cpm10k_by_context.csv")
