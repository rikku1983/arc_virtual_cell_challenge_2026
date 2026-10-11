# Arc Virtual Cell Challenge 2026

Zero-shot prediction of CRISPRi Perturb-seq responses in unseen cell contexts.
Challenge notes: [know_how/data.md](know_how/data.md). Data goes in `data/` (from `vcc datasets download`).

## Layout

| Path | Contents |
| --- | --- |
| `vccpred/data.py` | `Bundle` (manifest, gene order, perturbation list, lazy context loading) and `ContextData` |
| `vccpred/models/base.py` | `PerturbationModel`, the interface every model implements |
| `vccpred/models/` | one module per model, registered in `models/__init__.py` |
| `vccpred/submission.py` | builds the submission `.h5ad` and packages the `.vcc` with `vcc prep` |
| `vccpred/__main__.py` | command line: pick a model, predict, package |
| `explore/` | one-off data exploration scripts |
| `predictions/<run>/` | `config.json`, `prediction.h5ad`, `prediction.vcc` (gitignored) |

## Setup

```bash
uv venv --python 3.11 && uv pip install -r requirements.txt
uv tool install vcc-cli==0.2.2      # provides `vcc` (prep, submit)
vcc datasets download               # control bundle into data/
```

## Running

```bash
.venv/bin/python -m vccpred --model control_resample                  # full submission
.venv/bin/python -m vccpred --model control_resample --limit-perts 3  # smoke test
.venv/bin/python -m vccpred --model control_resample --param residual=0.3 --tag r30
vcc submit predictions/control_resample/prediction.vcc
```

## Adding a model

1. Create `vccpred/models/my_model.py` with a subclass of `PerturbationModel`:
   - set `name`
   - implement `predict(ctx, target_gene, n_cells, rng)`. It returns an `n_cells × 18,533` scipy CSR matrix of
     non-negative integer counts, with columns in `ctx.genes` order. `ctx.counts` holds the context's control cells.
   - optionally override `predict_context` to predict a whole context in one batch.
   - training is not part of this package. Train in your own workflow, then load the result in `__init__`
     (e.g. `--param checkpoint=path/to/model.pt`).
2. Add the class to the list in `vccpred/models/__init__.py`.
3. Run `.venv/bin/python -m vccpred --model my_model --param key=value ...`. Each `--param` is passed to `__init__`.

The full submission is built in memory before it is written (~17 GB of counts for the baseline, ~35 GB peak
while stacking), and `vcc prep` loads it again to validate and package it. Run full submissions on a machine
with at least 64 GB of RAM; `--limit-perts` smoke tests fit on a laptop.

## Models

| Name | Idea |
| --- | --- |
| `control_resample` | Resample 400 real control cells, then binomially thin the target gene's counts to `residual` (default 0.15) |
