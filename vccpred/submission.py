"""Writing a prediction to .h5ad and packaging it as a .vcc with `vcc prep`.

The whole submission is built in memory (~17 GB of counts for the baseline,
~35 GB peak while stacking), so run this on a machine with enough RAM.
"""
import subprocess
from collections.abc import Iterator
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from vccpred.data import Bundle


def write_prediction(path: Path, bundle: Bundle, blocks: Iterator[tuple[str, str, sp.csr_matrix]]) -> None:
    """Collect (context, target_gene, counts) blocks and write them as one .h5ad."""
    mats, contexts, perts = [], [], []
    for ctx, pert, X in blocks:
        mats.append(sp.csr_matrix(X, dtype=np.float32))
        contexts += [ctx] * X.shape[0]
        perts += [pert] * X.shape[0]

    obs = pd.DataFrame({"target_gene": pd.Categorical(perts), "context": pd.Categorical(contexts)})
    obs.index = obs.index.astype(str)
    adata = ad.AnnData(X=sp.vstack(mats, format="csr"), obs=obs, var=pd.DataFrame(index=bundle.genes))
    adata.write_h5ad(path)


def vcc_prep(h5ad: Path, bundle: Bundle, out: Path) -> None:
    """Validate the prediction and package it for upload with the official CLI."""
    cmd = ["vcc", "prep", str(h5ad), "-g", str(bundle.genes_csv), "--perts", str(bundle.perts_csv),
           "-o", str(out), "--force"]
    subprocess.run(cmd, check=True)
