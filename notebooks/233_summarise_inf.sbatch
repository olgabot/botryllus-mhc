#!/bin/bash
#SBATCH --job-name=233-summarise-inf
#SBATCH --partition=hns
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=233-summarise-inf-%j.out
#
# Notebook 233: reduce each alphabet's immune matches with E = inf (from
# 233_filter_immune.sbatch) to files small enough for the laptop. gbmr7 k = 13 alone has
# 474 M matches kmerseek did not extend, so the full files stay on Sherlock.
#
# For each alphabet, in $WORK/immune/:
#   <tag>.immune.evalue_inf.hla.csv.gz   every E = inf match to an HLA gene, all columns
#   <tag>.immune.evalue_inf.by_gene.tsv  per human gene: matches, Botryllus proteins, and
#                                         the lowest region_tail_probability (kmerseek
#                                         0.4.0 writes no corrected column)
#   <tag>.immune.evalue_inf.totals.tsv   matches, Botryllus proteins, human genes
#
# Submit after the filter job:
#   sbatch --dependency=afterok:<filter job id> notebooks/233_summarise_inf.sbatch
set -euo pipefail
: "${SCRATCH:?}"
IMMUNE=$SCRATCH/233-kmerseek-alphabets/immune

summarise_one() {
    set -euo pipefail
    local in=$1
    local stem=${in%.csv.gz}
    # the gzipped HLA file is written last, so it marks a finished alphabet
    [ -s "$stem.hla.csv.gz" ] && { echo "$(basename "$stem") already done"; return 0; }
    gzip -dc "$in" | awk -F, -v hla="$stem.hla.csv" -v by_gene="$stem.by_gene.tsv" -v totals="$stem.totals.tsv" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "query_name") q = i
                if ($i == "target_name") t = i
                if ($i == "region_tail_probability") p = i
            }
            # awk runs END even after exit, so END checks this flag before writing
            if (!q || !t || !p) { print "233: missing column" > "/dev/stderr"; bad = 1; exit 3 }
            print > hla; next
        }
        { split($t, f, "|"); g = f[7]; n++; nq[g]++
          if (!(($q, g) in seen)) { seen[$q, g] = 1; nprot[g]++ }
          if (!($q in queries)) { queries[$q] = 1; nqueries++ }
          if (!(g in best) || $p + 0 < best[g]) best[g] = $p + 0
          if (g ~ /^HLA-/) print > hla }
        END {
            if (bad) exit 3
            print "human_gene\tn_matches\tn_botryllus_proteins\tbest_region_tail_probability" > by_gene
            for (g in nq) { ngenes++; printf "%s\t%d\t%d\t%.6g\n", g, nq[g], nprot[g], best[g] > by_gene }
            print "n_matches\tn_botryllus_proteins\tn_human_genes" > totals
            printf "%d\t%d\t%d\n", n, nqueries, ngenes > totals
        }'
    gzip -f "$stem.hla.csv"
    echo "$(basename "$stem") done $(date)"
}
export -f summarise_one

ls "$IMMUNE"/*.immune.evalue_inf.csv.gz | xargs -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'summarise_one "$1"' _ {}
