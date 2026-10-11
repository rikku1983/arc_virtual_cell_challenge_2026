"""Loading the control bundle (contexts, gene order, perturbation list)."""
import json
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


@dataclass
class ContextData:
    """Non-targeting control cells of one cell context, as raw counts."""

    name: str
    counts: sp.csr_matrix  # cells × genes, columns in official gene order
    ntc_id: np.ndarray  # non-targeting guide of each cell
    genes: pd.Index

    @cached_property
    def gene_pos(self) -> dict[str, int]:
        return {g: i for i, g in enumerate(self.genes)}

    @property
    def n_cells(self) -> int:
        return self.counts.shape[0]


class Bundle:
    """The downloaded control bundle; contexts are loaded on demand (memory is tight)."""

    def __init__(self, data_dir: Path = DATA_DIR):
        self.dir = Path(data_dir)
        self.manifest = json.loads((self.dir / "manifest.json").read_text())
        self.contexts: list[str] = self.manifest["contexts"]
        self.cells_per_pert: int = self.manifest["cells_per_pert"]

    @property
    def genes_csv(self) -> Path:
        return self.dir / "gene_names.csv"

    @property
    def perts_csv(self) -> Path:
        return self.dir / "pert_counts.csv"

    @cached_property
    def genes(self) -> pd.Index:
        return pd.Index(pd.read_csv(self.genes_csv).iloc[:, 0].astype(str))

    @cached_property
    def perts(self) -> list[str]:
        return pd.read_csv(self.perts_csv)["target_gene"].astype(str).tolist()

    def load(self, ctx: str) -> ContextData:
        a = ad.read_h5ad(self.dir / f"context_{ctx}.h5ad")
        if not a.var_names.equals(self.genes):
            a = a[:, self.genes].copy()
        return ContextData(
            name=ctx,
            counts=sp.csr_matrix(a.X),
            ntc_id=a.obs["ntc_id"].astype(str).to_numpy(),
            genes=self.genes,
        )
