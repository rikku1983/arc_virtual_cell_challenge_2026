"""The interface every model implements."""
from abc import ABC, abstractmethod
from collections.abc import Iterator

import numpy as np
import scipy.sparse as sp

from vccpred.data import ContextData


class PerturbationModel(ABC):
    """Predicts raw-count profiles of CRISPRi-perturbed cells from a context's controls.

    Training happens outside this package: a trained model loads its weights or
    other artifacts in `__init__` (paths come in as --param values). Subclasses
    implement `predict`; override `predict_context` when it is cheaper to
    predict all perturbations of a context at once.
    """

    name: str = ""

    @abstractmethod
    def predict(
        self, ctx: ContextData, target_gene: str, n_cells: int, rng: np.random.Generator
    ) -> sp.csr_matrix:
        """Return n_cells × n_genes non-negative integer counts, columns in ctx.genes order."""

    def predict_context(
        self, ctx: ContextData, perts: list[str], n_cells: int, seed: int
    ) -> Iterator[tuple[str, sp.csr_matrix]]:
        """Yield (target_gene, counts) for every perturbation, in the order given."""
        for i, pert in enumerate(perts):
            rng = np.random.default_rng([seed, ord(ctx.name), i])
            yield pert, self.predict(ctx, pert, n_cells, rng)
