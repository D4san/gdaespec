#!/usr/bin/env bash
set -euo pipefail

RETRIEVAL_DIR="/mnt/c/Proyectos/Astro/gdaespec/training-and-experiments/Earth_like_Atmosphere/Retrieval Tests"
export GDAE_CAMPAIGN_DIR="${GDAE_CAMPAIGN_DIR:-/home/dasan/gdae_campaign_5obs_trial}"
export GDAE_MPI_LAUNCHER="${GDAE_MPI_LAUNCHER:-/home/dasan/anaconda3/envs/POSEIDON/bin/mpirun}"

cd "$RETRIEVAL_DIR"
source /home/dasan/anaconda3/etc/profile.d/conda.sh
conda activate POSEIDON

mkdir -p "$GDAE_CAMPAIGN_DIR/logs"
python -u campaign_run_gdae_queue.py --nproc 12 --include-test01 --rerun --keep-going > "$GDAE_CAMPAIGN_DIR/logs/gdae_queue_master.log" 2>&1
