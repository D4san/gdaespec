#!/usr/bin/env bash
set -euo pipefail

cd "/mnt/c/Proyectos/Astro/gdaespec/training-and-experiments/Earth_like_Atmosphere/Retrieval Tests"
source /home/dasan/anaconda3/etc/profile.d/conda.sh
conda activate POSEIDON

stamp="$(date +%Y%m%d_%H%M%S)"
python campaign_run_contam_queue.py --nproc 12 --include-test01 --keep-going \
  > "campaign_5obs/logs/contam_queue_master_resume_${stamp}.log" 2>&1
