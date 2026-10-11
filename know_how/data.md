# About the Data — Virtual Cell Challenge 2026

This page describes the data behind the 2026 Challenge — what you receive, how it is structured, and the format your predictions must take.

**Key dates**

- Validation set released: **August 20, 2026**
- Final test set released: **October 22, 2026**

---

## Challenge Datasets

*New challenge datasets experimentally generated at Arc Institute.*

For the 2026 Challenge, Arc performed Perturb-seq experiments on **six cell lines** originating from different tissues, using **CRISPRi** and **10x Flex scRNA-seq**.

The six lines are split into two panels of three:

- **Validation round:** three cell lines, released as anonymized contexts **A, B, C**
- **Final test:** three *different* cell lines, released as anonymized contexts **D, E, F**

Which cell line a context is will not be published. What you receive instead is that context's unperturbed expression profile, which is the only thing your model needs to condition on.

### Zero-shot task

There is **no Challenge training set** this year. Given profiles from cells expressing non-targeting guides, plus the gene identifiers for the CRISPRi knockdowns to be predicted, your model must predict the post-perturbation expression profiles that would be measured after CRISPRi knockdown followed by 10x Flex single-cell profiling.

You may train on the H1 hESC data released for the 2025 Challenge, and on any other public or proprietary data you have the right to use (see [Arc Virtual Cell Atlas](#arc-virtual-cell-atlas)).

### What participants receive

- **Validation set** (released August 20, 2026): non-targeting control profiles for contexts A, B and C — 18,400 control cells per context, from 46 non-targeting guides at 400 cells each, with a median of ~20,000 UMIs per cell — plus the list of 300 gene perturbations to predict in every context. Predictions against this panel drive the live ranking leaderboard.
- **Final test set** (released October 22, 2026): non-targeting control profiles for three different cell lines, as contexts D, E and F, with their own perturbation list.
- **Final rankings are based solely on performance on the final test set.**

### Shape of one round

| Item | Value |
| --- | --- |
| Cell contexts | 3 (A, B, C in validation — or D, E, F in the final phase) |
| Perturbations | 300, the same set in all three contexts of the round |
| Cells to predict per perturbation | 400 |
| Cells in a complete submission | 360,000 (300 × 400 × 3 contexts) |
| Non-targeting control cells provided | 18,400 per context — 46 non-targeting guides × 400 cells, labelled in `ntc_id` |
| Genes | 18,533 |

Perturbations are labelled by gene symbol (e.g. `ADNP`). The validation and test rounds are scored independently — with different cell lines and different perturbation panels — so a validation score and a final score are not directly comparable in an absolute sense.

---

## File Formats

Control data is downloaded after registering and logging in, either from the website or with the `vcc` command-line tool. The test-phase bundle (available October 22, 2026) has the same structure, for contexts D, E and F.

### Control data — `controls.zip` (~630 MB)

One control file per cell context, plus reference files:

| File | Contents |
| --- | --- |
| `context_A.h5ad` | 18,400 non-targeting control cells × 18,533 genes, for cell context A |
| `context_B.h5ad` | the same, for cell context B |
| `context_C.h5ad` | the same, for cell context C |
| `gene_names.csv` | the 18,533 gene symbols, in the order your matrix must use |
| `pert_counts.csv` | the perturbations to predict |
| `manifest.json` | season, partition, panel id, contexts and per-context cell counts |

Each control file is an AnnData H5AD file holding **raw counts**.

**`obs`**

| cell index | `target_gene` | `context` | `ntc_id` |
| --- | --- | --- | --- |
| `A_ctrl_000000` | non-targeting | A | non-targeting-11 |
| `A_ctrl_000001` | non-targeting | A | non-targeting-11 |

`ntc_id` names which of the 46 non-targeting guides a control cell carries (400 cells per guide). `target_gene` is `non-targeting` on every row, so anything keying off that column is unaffected; use `ntc_id` to model guide-to-guide variation within the controls.

**`var`** — index is the gene name (e.g. `SAMD11`, `NOC2L`, …).

> **Control cells are model-only inputs.**
> Every cell in these files has `target_gene = "non-targeting"` — they are the unperturbed profile of that cell context, and they are what identifies the context to your model. Do **not** predict expression values for control cells; a submission that does so will be rejected.

### Perturbation list — `pert_counts.csv` (1.9 KB)

300 rows, one per perturbation.

| Field | Description |
| --- | --- |
| `target_gene` | Gene symbol targeted for perturbation. Predict every one of them, in every cell context, with exactly 400 cells each. |

Example values: `ABCD1`, `ACLY`, `ADNP`, `AGO1`, …

---

## Submission — `prediction.vcc`

Build **one** AnnData H5AD holding predicted expression profiles for all three contexts; `vcc prep` packages it into a single `.vcc` file for upload. There is one H5AD inside, not one per context — the `context` column separates them, using the same labels as the control files.

**`obs`**

| `target_gene` | `context` |
| --- | --- |
| ABCD1 | A |
| ABCD1 | B |
| ACLY | A |

- Gene symbols, not construct ids — `ADNP`, never `ADNP-1`.
- No non-targeting rows: every row is one of the perturbations you were asked to predict.
- Do not copy the downloaded control cells into your prediction.

**`var`** — all 18,533 genes, in the order given by `gene_names.csv`.

**`X`**

- Raw counts: non-negative, whole numbers, finite.
- No cell totalling more than 1,000,000 counts.
- Shape 360,000 rows × 18,533 columns: exactly the perturbations in `pert_counts.csv`, 400 cells each, in every context.
- Stored **sparse**: at most 4,750,000,000 stored entries (~13,200 per cell at this shape). Explicitly stored zeros count, so a dense array exceeds the cap on its own.

### Packaging

```bash
vcc prep prediction.h5ad -g gene_names.csv --perts pert_counts.csv -o prediction.vcc
```

`prep` validates the gene set, context labels, perturbation labels, cell counts and the raw-counts requirement locally, so problems surface in seconds rather than after a multi-gigabyte upload. Add `--dry-run` to check without writing a file.

---

## Target Gene Selection

All genes targeted for perturbation prediction had **>80% on-target knockdown**, and were selected to provide a strong set of perturbations and responses for a robust Challenge.

---

## Arc Virtual Cell Atlas

You may use any public or proprietary data you have the right to use to train your model. The Arc Virtual Cell Atlas is a collection of high-quality, curated, open datasets assembled to accelerate the creation of virtual cell models. It includes both observational and perturbational data from over 600 million cells (and growing).

### scBaseCount (Arc Institute)

A continuously updated single-cell RNA-seq database that uses an AI-driven, hierarchical agent workflow to automate discovery, metadata extraction and standardized preprocessing of Sequence Read Archive (SRA) data.

- Largest public repository of single-cell data: over 502 million cells (and expanding), spanning 27 organisms and 75 tissues.
- By continually discovering, annotating and reprocessing raw scRNA-seq data, it offers a harmonized repository that can serve as a foundation for AI-driven modeling and integrative meta-analyses.

### Tahoe-100M (Tahoe)

The world's largest single-cell dataset generated and open-sourced by Tahoe: 100 million cells from ~60,000 drug perturbation experiments, mapping the response of 50 cancer models to 1,100+ drug treatments.

### 2025 Virtual Cell Challenge (Arc Institute)

The dataset generated for Arc's 2025 Virtual Cell Challenge: single-cell responses to perturbations in a human embryonic stem cell line (**H1 hESC**) with deep cell coverage, deep sequencing, high on-target knockdown and the **same single-cell chemistry (10x Flex)** as this year's competition. Includes the training, validation and held-out test perturbation datasets used throughout that Challenge.

---

## Public Perturbation Datasets

Beyond the Arc Virtual Cell Atlas, these publicly available perturbation datasets may be useful:

| Dataset | Cell lines / downloads | Description |
| --- | --- | --- |
| **Replogle et al., 2022** | K562 Genome-wide (61.3 GB), K562 (9.9 GB), RPE1 (8.1 GB) | Genome-scale Perturb-seq targeting all expressed genes with CRISPRi across >2.5 million human cells (K562 and RPE1). The K562 genome-wide dataset overlaps most of the genes used in the Arc VCC training and validation datasets. |
| **Nadig et al., 2025** | HepG2 (5.2 GB), Jurkat (8.7 GB) | Single-cell CRISPR screens of DepMap Common Essential Genes in Jurkat and HepG2 cells. |
| **Jiang et al., 2025** | A549, MCF7, HT29, HAP1, BxPC3, K562 | Perturb-seq in six cancer cell lines from different tissues of origin: A549 (lung), MCF7 (breast), HT29 (colon), HAP1 (bone marrow), BxPC3 (pancreas), K562 (bone marrow). |
| **Srivatsan et al., 2020** | 3 cancer cell lines | Introduces sci-Plex, which uses nuclear hashing to quantify global transcriptional responses to thousands of independent perturbations at single-cell resolution; screens three cancer cell lines exposed to 188 compounds. |
| **McFaline-Figueroa et al., 2024** | Glioblastoma cell lines | Introduces sci-Plex-Gene-by-Environment, a platform for combined single-cell genetic and chemical screening at scale, applied to combinations of chemical and genetic perturbations in glioblastoma cell lines. |
