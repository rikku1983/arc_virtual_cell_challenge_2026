#!/bin/bash
# Run the control_resample baseline. Submit from the repo root:
#   sbatch models/control_resample/predict.sbatch [predict.py args...]
#SBATCH --job-name=control_resample
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out
# TODO: set these for your cluster, then remove the extra '#'
##SBATCH --partition=TODO
##SBATCH --account=TODO
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${VCC_ENV:-vcc}"
echo "host=$(hostname) date=$(date -Is) commit=$(git rev-parse --short HEAD)"

python models/control_resample/predict.py "$@"
