#!/bin/bash
# ESMFold on the 500-protein pilot, with the image and arguments nf-core/proteinfold 2.0.0
# uses (quay.io/nf-core/proteinfold_esmfold:2.0.0, --num-recycles 4, no --chunk-size).
#
# The pipeline runs one task per protein. Here one job folds all 500 so the GPU is not
# requested 500 times; --max-tokens-per-batch 0 makes esm-fold fold one protein at a time,
# as the pipeline does. esm-fold logs each protein's fold time, and the time from job
# start to the first fold is the model load the pipeline pays once per protein.
#
# Usage: sbatch esmfold.sbatch <pilot dir> <fasta>
#SBATCH --job-name=pilot-esmfold
#SBATCH --partition=gpu
#SBATCH -G 1
#SBATCH -C GPU_MEM:32GB|GPU_MEM:48GB|GPU_MEM:80GB
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
set -euo pipefail
PILOT=$1
FASTA=$2
OUT=$PILOT/esmfold
mkdir -p "$OUT/pdb"

nvidia-smi --query-gpu=name,memory.total --format=csv,noheader > "$OUT/gpu.txt"
date +%s > "$OUT/start_epoch.txt"

apptainer exec --nv \
    --env TRITON_CACHE_DIR=/tmp/triton_cache --env XDG_CACHE_HOME=/tmp \
    --bind "$PILOT" --bind "$(dirname "$FASTA")" \
    "$PILOT/containers/proteinfold_esmfold_2.0.0.sif" \
    esm-fold -i "$FASTA" -o "$OUT/pdb" -m "$PILOT/esmfold_params" \
        --num-recycles 4 --max-tokens-per-batch 0 \
    2>&1 | tee "$OUT/esm-fold.log"

date +%s > "$OUT/end_epoch.txt"
