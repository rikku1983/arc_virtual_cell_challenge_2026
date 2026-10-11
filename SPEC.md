# Data and prediction spec

A model can be written in any language, live in any conda env, and run anywhere. It only has to read the input
below and write a prediction file that meets the output spec. `./package.sh` then checks and packages that file.

## Input: `data/`

`data/` holds the unzipped control bundle (`vcc datasets download controls`, file `vcc_2026_controls.zip`):

| File | Contents |
| --- | --- |
| `context_A.h5ad`, `context_B.h5ad`, `context_C.h5ad` | 18,400 non-targeting control cells × 18,533 genes each, raw integer counts (sparse, float32) |
| `gene_names.csv` | header `gene_name`, then the 18,533 gene symbols, **in the order a prediction must use** |
| `pert_counts.csv` | column `target_gene`: the 300 genes knocked down (CRISPRi) to predict, the same list in every context |
| `manifest.json` | `contexts` (`["A","B","C"]`), `cells_per_pert` (400), `panel_id` (`vcc2026-val-1`), per-context counts |

Each control file's `obs` has `target_gene` (always `non-targeting`), `context`, and `ntc_id`. `ntc_id` is one of
46 non-targeting guides, with 400 cells each. Its `var` index is the gene symbols. Median depth is about 20k UMIs
per cell. More detail is in [know_how/data.md](know_how/data.md).

**Final round** (released 2026-10-22): same structure, but contexts D, E, F, three different cell lines, and its
own perturbation list. Read contexts and counts from `manifest.json` rather than hard-coding them.

## Output: `predictions/<run>/prediction.h5ad`

One AnnData `.h5ad` covering every context:

| Part | Requirement |
| --- | --- |
| `obs` | columns `target_gene` (a symbol from `pert_counts.csv`) and `context` (a label from `manifest.json`). Exactly `cells_per_pert` (400) cells for every (context, perturbation) pair, and every perturbation in every context. **No `non-targeting` rows**, and don't copy the control cells. Row order is free. |
| `var` | index = the 18,533 symbols, **in `gene_names.csv` order** |
| `X` | **raw counts**: non-negative, whole numbers, finite (float32 storage is fine). No cell totals more than 1,000,000 counts. |
| `X` storage | **sparse** (CSR). Stored entries are capped by the scoring hardware (`vcc prep --max-nnz`, about 4.75e9 at full size). Explicitly stored zeros count, so a dense matrix fails. |
| shape | 360,000 × 18,533 for the validation round (3 contexts × 300 perturbations × 400 cells) |

Anything that changes values is **the model's job**: turning log-normalized output into integer counts,
making the matrix sparse, putting genes in the official order, filling genes the model doesn't predict.

Also recommended: a `config.json` in the same run folder recording the model, parameters, seed and code version.

For scale: the baseline has about 2.1e9 stored counts, a 16 GB `.h5ad` and a 3 GB `.vcc`.
