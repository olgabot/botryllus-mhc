#!/bin/bash
# Foldseek search of the ProstT5 query database against AlphaFold/Swiss-Prot.
# The queries have 3Di letters and amino acids but no coordinates, so the
# default 3Di+AA alignment is used. No TM-score, lDDT or `prob` column: each
# needs the query C-alpha coordinates (convertalis stops on a missing _ca file).
#   sbatch ... foldseek_search.sbatch <query_db> <afdb_swissprot_db> <out.tsv>
#SBATCH --job-name=bs-sp-foldseek
#SBATCH --output=logs/%x.%j.out
set -euo pipefail

QDB=$1
DB=$2
OUT=$3

FOLDSEEK=${FOLDSEEK:-$GROUP_HOME/opt/foldseek-10/foldseek/bin/foldseek}
"$FOLDSEEK" version
hostname
df -h "$L_SCRATCH_JOB"

TMP=$L_SCRATCH_JOB/foldseek_tmp
mkdir -p "$TMP"

"$FOLDSEEK" search "$QDB" "$DB" "$TMP/aln" "$TMP" \
    -e 1e-3 --max-seqs 2000 ${EXTRA_ARGS:-} \
    --threads "$SLURM_CPUS_PER_TASK"
"$FOLDSEEK" convertalis "$QDB" "$DB" "$TMP/aln" "$L_SCRATCH_JOB/out.tsv" \
    --format-mode 4 \
    --format-output "query,target,theader,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov${EXTRA_COLUMNS:+,${EXTRA_COLUMNS// /,}}" \
    --threads "$SLURM_CPUS_PER_TASK"

cp "$L_SCRATCH_JOB/out.tsv" "$OUT"
wc -l "$OUT"
