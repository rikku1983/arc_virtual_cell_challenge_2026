# Arc Virtual Cell Challenge 2026

Zero-shot prediction of CRISPRi Perturb-seq responses in unseen cell contexts.
Challenge notes: [know_how/data.md](know_how/data.md).

This repo holds the data spec, the prediction spec, and the packaging and submission steps.
Models are standalone: any code, any conda env. Each one only has to write a `prediction.h5ad` that meets
[SPEC.md](SPEC.md).

## Layout

| Path | Contents |
| --- | --- |
| [SPEC.md](SPEC.md) | what `data/` contains, and what a `prediction.h5ad` must look like |
| `package.sh` | `prediction.h5ad` → `prediction.vcc` with `vcc prep`, settings from `data/manifest.json` |
| `models/<name>/` | optional home for model code; `control_resample/` is the baseline and an example |
| `explore/` | one-off data exploration scripts |
| `data/` | the unzipped control bundle (gitignored) |
| `predictions/<run>/` | `prediction.h5ad`, `config.json`, `prediction.vcc` (gitignored) |

## Setup

```bash
conda create -n vcc python=3.11 -y && conda activate vcc
pip install -r models/control_resample/requirements.txt   # only needed for the baseline
pip install vcc-cli==0.2.2                                # `vcc` (prep, submit)
vcc login --token-stdin --store-plaintext                 # or copy ~/.config/vcc/credentials.json
vcc datasets download controls -d data && python -m zipfile -e data/vcc_2026_controls.zip data/
```

## Running: three separate steps

```bash
# 1. Predict, with any model in its own env. The baseline as a Slurm job:
sbatch models/control_resample/predict.sbatch                 # log: slurm-control_resample-<jobid>.out

# 2. Package (login node): validates against SPEC.md and writes prediction.vcc next to the h5ad
./package.sh predictions/control_resample/prediction.h5ad
./package.sh predictions/control_resample/prediction.h5ad --dry-run   # check only

# 3. Submit (login node, by hand: submissions per day are limited)
vcc submit predictions/control_resample/prediction.vcc -m control_resample -d "what changed" --wait
vcc submit predictions/control_resample/prediction.vcc -m control_resample --resume   # if the upload was interrupted
```

`vcc prep` loads the whole prediction into memory: about 20–35 GB of RAM for a full one. If the login node
kills step 2, run it inside an interactive allocation (`srun --mem=64G --pty bash`).

For a new model, only step 1 changes. Write `predictions/<run>/prediction.h5ad`, then run steps 2–3 on it.

## Models

| Name | Idea |
| --- | --- |
| `control_resample` | Resample 400 real control cells per perturbation, then binomially thin the target gene's counts to `--residual` (default 0.15). Smoke test: `python models/control_resample/predict.py --limit-perts 2 --out predictions/smoke` |
