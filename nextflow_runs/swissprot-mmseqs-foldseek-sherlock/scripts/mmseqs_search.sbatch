#!/bin/bash
# MMseqs2 search of Botryllus proteins against UniProtKB/Swiss-Prot.
# Time, CPUs and memory come from the sbatch command line in the Makefile.
#   sbatch ... mmseqs_search.sbatch <query.fa> <swissprot_mmseqs_db> <out.tsv>
#SBATCH --job-name=bs-sp-mmseqs
#SBATCH --output=logs/%x.%j.out
set -euo pipefail

QUERY=$1
DB=$2
OUT=$3

ml purge
ml biology mmseqs2/18-8cc5c
mmseqs version
hostname
df -h "$L_SCRATCH_JOB"

# Temporary files go to the node-local disk, as /etc/agents/storage.md asks.
TMP=$L_SCRATCH_JOB/mmseqs_tmp
mkdir -p "$TMP"

# --max-seqs 2000: Botryllus Ig-domain proteins match hundreds of Swiss-Prot
# antibodies and receptors; with the default 300 a weaker MHC match could be
# dropped before alignment. EXTRA_ARGS adds options and EXTRA_COLUMNS output
# columns (see run-mhc-exhaustive).
mmseqs easy-search "$QUERY" "$DB" "$L_SCRATCH_JOB/out.tsv" "$TMP" \
    -s 7.5 -e 1e-3 --max-seqs 2000 ${EXTRA_ARGS:-} \
    --threads "$SLURM_CPUS_PER_TASK" \
    --format-mode 4 \
    --format-output "query,target,theader,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov${EXTRA_COLUMNS:+,${EXTRA_COLUMNS// /,}}"

cp "$L_SCRATCH_JOB/out.tsv" "$OUT"
wc -l "$OUT"
