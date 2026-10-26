#!/bin/bash
#SBATCH --job-name=233-kmerseek-alphabets
#SBATCH --account=ayeletv
#SBATCH --partition=hns
#SBATCH --array=0-18
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=233-kmerseek-alphabets-%A_%a.out
#
# Notebook 233 on Sherlock: one array task per row of 233_arms.tsv (19 kmerseek alphabets,
# each at the k that carries at least the bits of hp_lehninger2 k = 24). Each task indexes
# human GENCODE v49 canonical proteins, fits E-values, searches all 72,617 Botryllus
# proteins and keeps the matches with E <= 10 or E = inf (233_kmerseek_alphabet_arm.py).
#
# kmerseek is v0.4.0, the release notebook 232's run used. Build it first, once:
#
#   sbatch notebooks/233_build_kmerseek.sbatch
#
# then submit this array after it (or alone, if the binary already exists):
#
#   sbatch --dependency=afterok:<build job id> notebooks/233_kmerseek_alphabets.sbatch
#
# Sizing: the only measured arm is hp_lehninger2 k = 24 on the laptop (PR 10, 22 Sept 2026):
# index 10 s, search 1 min 16 s, peak memory not recorded. The other alphabets make more
# chance regions at the same bits (notebook 241's fit: protein20 k = 6 had 3.1 M regions
# against 37 thousand for hp_lehninger2 k = 24), so a search can take up to ~100 times as
# long. 12 h and 64 GB cover that; read the real use afterwards with
#   sacct -j <array job id> -o JobID%20,State,Elapsed,MaxRSS
# and size the next run from it.
#
# Copy back only the small files:
#   rsync -av sherlock-dtn:/scratch/users/olgabot/233-kmerseek-alphabets/{evalue_le_10,ka_survival} \
#       /Users/olga/data/botryllus/kmerseek-0.4-alphabets/
set -euo pipefail

: "${SCRATCH:?}" "${SLURM_SUBMIT_DIR:?}" "${SLURM_ARRAY_TASK_ID:?}"
WORK=$SCRATCH/233-kmerseek-alphabets
export KMERSEEK=${KMERSEEK:-$SCRATCH/kmerseek-v0.4.0/target/release/kmerseek}
export HUMAN=${HUMAN:-$SCRATCH/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa}
export QUERY=${QUERY:-$WORK/Bs_proteins.fa}
export OUT_DIR=${OUT_DIR:-$WORK}
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK

for f in "$KMERSEEK" "$HUMAN" "$QUERY" "$SLURM_SUBMIT_DIR/notebooks/233_arms.tsv"; do
    [ -s "$f" ] || { echo "233: missing input $f" >&2; exit 2; }
done

# gcc for libstdc++ when kmerseek runs, python for the arm script (standard library only).
module load gcc/14.2.0 python/3.12.1
"$KMERSEEK" --version
python3 "$SLURM_SUBMIT_DIR/notebooks/233_kmerseek_alphabet_arm.py" "$SLURM_ARRAY_TASK_ID"
echo "233: task $SLURM_ARRAY_TASK_ID done. Time and memory: sacct -j ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID -o JobID,Elapsed,MaxRSS"
