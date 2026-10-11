"""Baseline: perturbed cells look like control cells, with the target gene knocked down.

For every (context, perturbation), resample real control cells of that context and
binomially thin the target gene's counts: each count survives with p=residual
(CRISPRi knockdown is >80%, so the default keeps 15%). Writes a prediction that
meets SPEC.md.

  python models/control_resample/predict.py                  # full prediction
  python models/control_resample/predict.py --limit-perts 2  # smoke test
"""
import argparse
import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def knock_down(X: sp.csr_matrix, col: int, residual: float, rng: np.random.Generator) -> None:
    """Binomial thinning of one column in place."""
    hit = X.indices == col
    X.data[hit] = rng.binomial(X.data[hit].astype(np.int64), residual)
    X.eliminate_zeros()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=Path("data"))
    ap.add_argument("--out", type=Path, default=Path("predictions/control_resample"))
    ap.add_argument("--residual", type=float, default=0.15, help="fraction of target-gene counts kept")
    ap.add_argument("--exclude-guides", action="append", default=[], metavar="NTC_ID",
                    help="never sample cells of this non-targeting guide (repeatable)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit-perts", type=int, default=None, help="predict only the first N perturbations")
    args = ap.parse_args()

    manifest = json.loads((args.data_dir / "manifest.json").read_text())
    genes = pd.Index(pd.read_csv(args.data_dir / "gene_names.csv").iloc[:, 0].astype(str))
    perts = pd.read_csv(args.data_dir / "pert_counts.csv")["target_gene"].astype(str).tolist()
    perts = perts[: args.limit_perts]
    n_cells = manifest["cells_per_pert"]
    gene_pos = {g: i for i, g in enumerate(genes)}

    mats, obs_ctx, obs_pert = [], [], []
    for ctx in manifest["contexts"]:
        print(f"predict context {ctx}", flush=True)
        a = ad.read_h5ad(args.data_dir / f"context_{ctx}.h5ad")
        if not a.var_names.equals(genes):
            a = a[:, genes].copy()
        counts = sp.csr_matrix(a.X)
        pool = np.flatnonzero(~a.obs["ntc_id"].astype(str).isin(args.exclude_guides).to_numpy())
        del a

        for i, pert in enumerate(perts):
            rng = np.random.default_rng([args.seed, ord(ctx), i])
            rows = rng.choice(pool, size=n_cells, replace=n_cells > len(pool))
            X = counts[rows].astype(np.float32)  # copy; counts stays untouched
            if args.residual < 1.0:
                knock_down(X, gene_pos[pert], args.residual, rng)
            mats.append(X)
            obs_ctx += [ctx] * n_cells
            obs_pert += [pert] * n_cells
        del counts

    obs = pd.DataFrame({"target_gene": pd.Categorical(obs_pert), "context": pd.Categorical(obs_ctx)})
    obs.index = obs.index.astype(str)
    pred = ad.AnnData(X=sp.vstack(mats, format="csr"), obs=obs, var=pd.DataFrame(index=genes))
    del mats

    args.out.mkdir(parents=True, exist_ok=True)
    pred.write_h5ad(args.out / "prediction.h5ad")
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    config = {"model": "control_resample", "residual": args.residual, "exclude_guides": args.exclude_guides,
              "seed": args.seed, "n_perts": len(perts), "panel_id": manifest["panel_id"], "git_commit": commit}
    (args.out / "config.json").write_text(json.dumps(config, indent=1))
    print(f"wrote {args.out / 'prediction.h5ad'}: {pred.n_obs} cells, {pred.X.nnz:,} stored counts")


if __name__ == "__main__":
    main()
