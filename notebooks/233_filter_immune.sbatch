#!/bin/bash
#SBATCH --job-name=233-filter-immune
#SBATCH --partition=hns
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=233-filter-immune-%j.out
#
# Notebook 233: keep only the matches to human immune genes from each alphabet's
# E <= 10 / E = inf file (233_kmerseek_alphabets.sbatch). Some of those files are too big
# to copy to the laptop (gbmr7 k = 13: 475 M matches, 69 GB gzipped). Immune genes are the
# QuickGO GO:0002376 list notebook 232 uses (notebooks/data/, copied next to this script).
#
# For each alphabet it writes, in $WORK/immune/:
#   <tag>.immune.evalue_le_10.csv.gz   immune matches with a finite E <= 10
#   <tag>.immune.evalue_inf.csv.gz     immune matches kmerseek did not extend (E = inf)
# The human gene symbol is field 7 of target_name, split on "|" (GENCODE header).
#
# Submit from $SCRATCH/233-kmerseek-alphabets:
#   sbatch notebooks/233_filter_immune.sbatch
# Sizing: gzip -dc reads about 100-200 MB/s per core; the biggest file is 363 GB
# uncompressed, so ~1 h for it; the 19 files run 8 at a time.
set -euo pipefail
: "${SCRATCH:?}" "${SLURM_SUBMIT_DIR:?}"
WORK=$SCRATCH/233-kmerseek-alphabets
GENES=$SLURM_SUBMIT_DIR/notebooks/data/go0002376_immune_system_process_human_quickgo.tsv
mkdir -p "$WORK/immune"
[ -s "$GENES" ] || { echo "233: missing $GENES" >&2; exit 2; }

filter_one() {
    local in=$1 genes=$2 out_dir=$3
    local tag
    tag=$(basename "$in" .evalue_le_10_or_inf.csv.gz)
    [ -s "$out_dir/$tag.immune.evalue_inf.csv.gz" ] && { echo "$tag already done"; return 0; }
    gzip -dc "$in" | awk -F, -v genes="$genes" \
        -v fin="$out_dir/$tag.immune.evalue_le_10.csv" -v finf="$out_dir/$tag.immune.evalue_inf.csv" '
        BEGIN { while ((getline line < genes) > 0) { split(line, g, "\t"); immune[g[1]] = 1 } }
        NR == 1 {
            for (i = 1; i <= NF; i++) { if ($i == "target_name") t = i; if ($i == "region_evalue") e = i }
            if (!t || !e) { print "233: header has no target_name or region_evalue" > "/dev/stderr"; exit 3 }
            print > fin; print > finf; next
        }
        { split($t, f, "|"); if (!(f[7] in immune)) next
          if ($e == "inf" || $e == "Infinity") print > finf; else print > fin }'
    gzip -f "$out_dir/$tag.immune.evalue_le_10.csv"
    gzip -f "$out_dir/$tag.immune.evalue_inf.csv"
    echo "$tag done $(date)"
}
export -f filter_one

ls "$WORK"/evalue_le_10/*.evalue_le_10_or_inf.csv.gz | \
    xargs -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'filter_one "$1" "$2" "$3"' _ {} "$GENES" "$WORK/immune"
ls -la "$WORK/immune"
