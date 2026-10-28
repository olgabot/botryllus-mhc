#!/bin/bash
# Notebook 234 on Sherlock: notebook 233's 19-alphabet kmerseek search, repeated with the
# new Botryllus proteome (45,339 FUN proteins, botryllus_new_genome_protein.fasta, 23 Oct
# 2026) as the query. Same target (human GENCODE v49 canonical), same alphabets and k
# (233_arms.tsv), same kmerseek 0.4.0 and settings; only the query and the folder change.
#
# Run on a Sherlock login node from $SCRATCH/234-kmerseek-alphabets-fun, after copying
# there: this notebooks/ folder's 233_* scripts, 233_arms.tsv,
# data/go0002376_immune_system_process_human_quickgo.tsv, and the proteome as
# botryllus_new_genome_protein.fasta. It only submits jobs.
#
# Sizing, from notebook 233's run (array 49626308, 72,617 proteins, 17.0 M residues; the
# FUN proteome has 16.9 M residues, so about the same work):
#   gbmr7 (row 12)              9 h 6 min, peak memory 64 GB = the cap   -> 16 h, 128 GB
#   gbmr4 (row 9)               2 h 28 min, peak 64 GB = the cap          -> 16 h, 128 GB
#   the other 17 alphabets      155 min in all, peak 26 GB per task
#                               -> 2 tasks of ~78 min each                -> 4 h, 48 GB
#   immune filter (49792907)    56 min, 16 GB peak                       -> 6 h, 32 GB
#   E = inf summary (49793305)  7 min, 1.9 GB peak                       -> as written
set -euo pipefail
: "${SCRATCH:?}"
export WORK=$SCRATCH/234-kmerseek-alphabets-fun
export QUERY=$WORK/botryllus_new_genome_protein.fasta
export OUT_DIR=$WORK
cd "$WORK"
mkdir -p logs

# rows of 233_arms.tsv; ';' separates array tasks
heavy=$(ROW_GROUPS="12;9" sbatch --parsable --export=ALL --job-name=234-kmerseek-heavy \
    --array=0-1 --mem=128G --time=16:00:00 \
    --output=logs/234-kmerseek-heavy-%A_%a.out notebooks/233_kmerseek_alphabets.sbatch)
light=$(ROW_GROUPS="0,4,6,10,13,15,16,18;1,2,3,5,7,8,11,14,17" sbatch --parsable --export=ALL \
    --job-name=234-kmerseek-light --array=0-1 --mem=48G --time=04:00:00 \
    --output=logs/234-kmerseek-light-%A_%a.out notebooks/233_kmerseek_alphabets.sbatch)
filter=$(sbatch --parsable --export=ALL --job-name=234-filter-immune --mem=32G \
    --dependency=afterok:$heavy:$light \
    --output=logs/234-filter-immune-%j.out notebooks/233_filter_immune.sbatch)
summary=$(sbatch --parsable --export=ALL --job-name=234-summarise-inf \
    --dependency=afterok:$filter \
    --output=logs/234-summarise-inf-%j.out notebooks/233_summarise_inf.sbatch)
echo "heavy $heavy light $light filter $filter summary $summary"
