#!/bin/bash
# Notebook 235 on Sherlock: notebook 234's kmerseek search (FUN Botryllus proteome against
# human GENCODE v49 canonical, 19 alphabets, kmerseek 0.4.0), with each alphabet at the k
# carrying at least the bits of protein20 k = 7, bits counted from Swiss-Prot residue
# frequencies (235_arms.tsv, from 235_plan_arms.py). Only the k values and the folder
# change; scripts are notebook 233's, reading the arm table from ARMS_TSV.
#
# Run on a Sherlock login node from $SCRATCH/235-kmerseek-alphabets-fun-k7, after copying
# there: this notebooks/ folder's 233_* scripts, 235_arms.tsv,
# data/go0002376_immune_system_process_human_quickgo.tsv, and the cleaned FUN proteome
# (234_prepare_fun_proteome.py output) as botryllus_new_genome_protein.fasta. It only
# submits jobs.
#
# Sizing, from notebook 234's run (same query and target, smaller k). Larger k means fewer
# chance matches, and most of a task's time goes on writing matches (gbmr7 k = 13: search
# 1 h 25 min, task 9 h 34 min for 524 M matches), so these are upper bounds:
#   gbmr7 k 23 (row 12)        was 9 h 34 min, peak 128 GB = the cap    -> 12 h, 128 GB
#   gbmr4 k 22 (row 9)         was 2 h 13 min, peak 111 GB              -> 12 h, 128 GB
#   the other 17, 2 tasks      was 1 h 32 min and 1 h 16 min, peak 44 GB -> 4 h, 64 GB
#   filter + summary, 1 job    was 56 min + 10 min, peak 16.5 GB        -> 4 h, 32 GB
set -euo pipefail
: "${SCRATCH:?}"
export WORK=$SCRATCH/235-kmerseek-alphabets-fun-k7
export QUERY=$WORK/botryllus_new_genome_protein.fasta
export OUT_DIR=$WORK
export ARMS_TSV=$WORK/notebooks/235_arms.tsv
cd "$WORK"
mkdir -p logs

# Rerun, 28 Oct 2026. Jobs 49951877_0, 49951880_0 and 49951880_1 failed for six alphabets
# (gbmr7, hp_kyte_doolittle2, funcgroups8, hp_lehninger_hpc3, hsdm17, sdm12): at these k the
# E-value fit found too few scored regions even with 1,000 calibration sequences, so the
# search refused to run. The other 13 finished. The rerun calibrates on 19,000 of the
# ~19,500 human proteins. Measured on the first try: every task peaked under 9 GB and ran
# under 17 min for up to 9 alphabets, so one task with 16 GB and 2 h covers the six.
# Run as: bash notebooks/235_submit_protein20_k7_swissprot.sh rerun
if [ "${1:-}" = rerun ]; then
    rerun=$(ROW_GROUPS="12,4,13,7,17,14" RETRY_KA_QUERIES=19000 sbatch --parsable --export=ALL \
        --job-name=235-kmerseek-rerun --array=0 --mem=16G --time=02:00:00 \
        --output=logs/235-kmerseek-rerun-%A_%a.out notebooks/233_kmerseek_alphabets.sbatch)
    immune=$(sbatch --parsable --export=ALL --job-name=235-immune --partition=hns \
        --cpus-per-task=8 --mem=32G --time=04:00:00 --dependency=afterany:$rerun \
        --output=logs/235-immune-%j.out \
        --wrap="bash notebooks/233_filter_immune.sbatch && bash notebooks/233_summarise_inf.sbatch")
    echo "rerun $rerun immune $immune"
    exit 0
fi


# rows of 235_arms.tsv (same order as 233_arms.tsv); ';' separates array tasks
heavy=$(ROW_GROUPS="12;9" sbatch --parsable --export=ALL --job-name=235-kmerseek-heavy \
    --array=0-1 --mem=128G --time=12:00:00 \
    --output=logs/235-kmerseek-heavy-%A_%a.out notebooks/233_kmerseek_alphabets.sbatch)
light=$(ROW_GROUPS="0,4,6,10,13,15,16,18;1,2,3,5,7,8,11,14,17" sbatch --parsable --export=ALL \
    --job-name=235-kmerseek-light --array=0-1 --mem=64G --time=04:00:00 \
    --output=logs/235-kmerseek-light-%A_%a.out notebooks/233_kmerseek_alphabets.sbatch)
# the summary alone takes ~10 min, so it runs in the filter's job
immune=$(sbatch --parsable --export=ALL --job-name=235-immune --partition=hns \
    --cpus-per-task=8 --mem=32G --time=04:00:00 --dependency=afterok:$heavy:$light \
    --output=logs/235-immune-%j.out \
    --wrap="bash notebooks/233_filter_immune.sbatch && bash notebooks/233_summarise_inf.sbatch")
echo "heavy $heavy light $light immune $immune"

