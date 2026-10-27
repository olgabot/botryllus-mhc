#!/bin/bash
#SBATCH --job-name=233-kmerseek-alphabets
#SBATCH --partition=hns
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=233-kmerseek-alphabets-%A_%a.out
#
# Notebook 233 on Sherlock: 19 kmerseek alphabets (233_arms.tsv), each at the k that
# carries at least the bits of hp_lehninger2 k = 24. For each one,
# 233_kmerseek_alphabet_arm.py indexes human GENCODE v49 canonical proteins, fits
# E-values, searches all 72,617 Botryllus proteins and keeps the matches with E <= 10 or
# E = inf.
#
# 4 array tasks; task t runs rows t, t+4, t+8, ... of 233_arms.tsv one after another, so
# the many-class alphabets (rows 14-18) and protein20 (row 0) are spread across tasks.
# One alphabet per task would be too short a job: hp_lehninger2 k = 24 took 1.5 min on
# the laptop, and Sherlock asks for at least 10 min of work per job.
#
# Submit from a checkout or copy of these scripts on Sherlock:
#   sbatch notebooks/233_kmerseek_alphabets.sbatch
#
# kmerseek is v0.4.0, the release notebook 232's run used: the build in
# $SCRATCH/301-kmerseek-v0.4.0 (or 233_build_kmerseek.sbatch, if that is gone).
#
# Sizing: the only measured arm is hp_lehninger2 k = 24 on the laptop (PR 10, 22 Sept 2026):
# index 10 s, search 1 min 16 s, peak memory not recorded. Other alphabets make more
# chance regions at the same bits (notebook 241's fit: protein20 k = 6 had 3.1 M regions
# against 37 thousand for hp_lehninger2 k = 24), so one search can take ~100 times as long
# and a task holds ~5 of them. 24 h and 64 GB cover that; read the real use afterwards with
#   sacct -j <array job id> -o JobID%20,State,Elapsed,MaxRSS
# and size the next run from it.
#
# The index and the full search CSV live on node-local disk ($L_SCRATCH_JOB) and are
# deleted after each alphabet. Copy back only the small files:
#   rsync -av sherlock-dtn:/scratch/users/olgabot/233-kmerseek-alphabets/{evalue_le_10,ka_survival,logs} \
#       /Users/olga/data/botryllus/kmerseek-0.4-alphabets/
# The FUN proteome run (notebook 234) uses the same script with other inputs and
# memory, and gives each task a balanced set of alphabets measured from this run:
#   gbmr7 (row 12) 9 h 6 min, gbmr4 (row 9) 2 h 28 min, the other 17 alphabets 155 min
#   in all; peak memory 64 GB (the cap) for the tasks holding gbmr7 and gbmr4, 26 GB or
#   less for the rest. See notebooks/234_submit_fun_proteome.sh.
set -euo pipefail

: "${SCRATCH:?}" "${SLURM_SUBMIT_DIR:?}" "${SLURM_ARRAY_TASK_ID:?}" "${L_SCRATCH_JOB:?}"
SCRIPTS=$SLURM_SUBMIT_DIR/notebooks
WORK=${WORK:-$SCRATCH/233-kmerseek-alphabets}
export KMERSEEK=${KMERSEEK:-$SCRATCH/301-kmerseek-v0.4.0/target/release/kmerseek}
export HUMAN=${HUMAN:-$SCRATCH/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa}
export QUERY=${QUERY:-$WORK/Bs_proteins.fa}
export OUT_DIR=${OUT_DIR:-$WORK}
export WORK_DIR=$L_SCRATCH_JOB/233
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK
N_TASKS=4
# ROW_GROUPS, if set, gives each array task its own rows of 233_arms.tsv: groups split by
# ';', rows by ','. Task t runs group t. Unset: task t runs rows t, t+4, t+8, ...
ROW_GROUPS=${ROW_GROUPS:-}

for f in "$KMERSEEK" "$HUMAN" "$QUERY" "$SCRIPTS/233_arms.tsv" "$SCRIPTS/233_kmerseek_alphabet_arm.py"; do
    [ -s "$f" ] || { echo "233: missing input $f" >&2; exit 2; }
done

# gcc for kmerseek's libstdc++, python for the arm script (standard library only).
module load gcc/14.2.0 python/3.12.1
"$KMERSEEK" --version
df -h "$L_SCRATCH"

n_rows=$(( $(wc -l < "$SCRIPTS/233_arms.tsv") - 1 ))
if [ -n "$ROW_GROUPS" ]; then
    IFS=';' read -ra groups <<< "$ROW_GROUPS"
    IFS=',' read -ra rows <<< "${groups[$SLURM_ARRAY_TASK_ID]:?no row group for task $SLURM_ARRAY_TASK_ID}"
else
    rows=()
    for (( row = SLURM_ARRAY_TASK_ID; row < n_rows; row += N_TASKS )); do rows+=("$row"); done
fi
failed=()
for row in "${rows[@]}"; do
    echo "233: row $row start $(date)"
    # one alphabet failing must not stop the others in this task
    python3 "$SCRIPTS/233_kmerseek_alphabet_arm.py" "$row" || failed+=("$row")
    rm -rf "$WORK_DIR"
    echo "233: row $row end $(date)"
done

echo "233: task $SLURM_ARRAY_TASK_ID done. Time and memory: sacct -j ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID -o JobID,Elapsed,MaxRSS"
if (( ${#failed[@]} )); then
    echo "233: rows that failed: ${failed[*]}" >&2
    exit 1
fi
