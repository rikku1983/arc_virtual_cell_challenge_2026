#!/bin/bash
# Package a prediction (see SPEC.md) into a .vcc next to it, with `vcc prep`:
#   ./package.sh predictions/<run>/prediction.h5ad [extra vcc prep args, e.g. --dry-run]
# Contexts and cells per perturbation come from $DATA_DIR/manifest.json (default data/).
# vcc prep loads the whole h5ad (~20-35 GB RAM for a full prediction). If the login
# node kills it, run it inside an interactive allocation: srun --mem=64G --pty bash
set -euo pipefail

h5ad=${1:?usage: ./package.sh PREDICTION.h5ad [vcc prep args...]}
shift
data=${DATA_DIR:-data}
read -r contexts cells_per_pert < <(python3 -c "
import json, sys
m = json.load(open(sys.argv[1]))
print(','.join(m['contexts']), m['cells_per_pert'])" "$data/manifest.json")

vcc prep "$h5ad" -g "$data/gene_names.csv" --perts "$data/pert_counts.csv" \
    --contexts "$contexts" --cells-per-pert "${CELLS_PER_PERT:-$cells_per_pert}" \
    -o "$(dirname "$h5ad")/prediction.vcc" --force "$@"
