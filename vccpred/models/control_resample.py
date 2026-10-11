"""Baseline: perturbed cells look like control cells, with the target gene knocked down."""
import numpy as np
import scipy.sparse as sp

from vccpred.data import ContextData
from vccpred.models.base import PerturbationModel


class ControlResample(PerturbationModel):
    """Resample real control cells of the context and thin the target gene's counts.

    residual: fraction of target-gene counts kept (CRISPRi knockdown is >80%, so
    the default keeps 15%). residual=1.0 gives the pure control baseline.
    exclude_guides: ntc_id values whose cells are never sampled.
    """

    name = "control_resample"

    def __init__(self, residual: float = 0.15, exclude_guides: tuple[str, ...] = ()):
        self.residual = residual
        self.exclude_guides = set(exclude_guides)

    def predict(self, ctx, target_gene, n_cells, rng):
        pool = np.flatnonzero(~np.isin(ctx.ntc_id, list(self.exclude_guides)))
        rows = rng.choice(pool, size=n_cells, replace=n_cells > len(pool))
        X = ctx.counts[rows].astype(np.float32)  # copy; ctx.counts stays untouched
        if self.residual < 1.0 and target_gene in ctx.gene_pos:
            self._knock_down(X, ctx.gene_pos[target_gene], rng)
        return X

    def _knock_down(self, X: sp.csr_matrix, col: int, rng: np.random.Generator) -> None:
        """Binomial thinning of one column in place: each count survives with p=residual."""
        hit = X.indices == col
        X.data[hit] = rng.binomial(X.data[hit].astype(np.int64), self.residual)
        X.eliminate_zeros()
