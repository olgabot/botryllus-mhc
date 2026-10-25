#!/bin/bash
# Predict 3Di structure letters for each Botryllus protein from its sequence
# with ProstT5, and write a Foldseek query database. Runs on a GPU.
#   sbatch ... foldseek_prostt5_createdb.sbatch <query.fa> <prostt5_weights> <out_db_prefix>
#SBATCH --job-name=bs-sp-prostt5
#SBATCH --output=logs/%x.%j.out
set -euo pipefail

QUERY=$1
WEIGHTS=$2
OUTDB=$3

FOLDSEEK=${FOLDSEEK:-$GROUP_HOME/opt/foldseek-10/foldseek/bin/foldseek}
"$FOLDSEEK" version
hostname
nvidia-smi --query-gpu=name,memory.total --format=csv

mkdir -p "$(dirname "$OUTDB")"
"$FOLDSEEK" createdb "$QUERY" "$OUTDB" \
    --prostt5-model "$WEIGHTS" --gpu 1 \
    --threads "$SLURM_CPUS_PER_TASK"
ls -la "$(dirname "$OUTDB")"
