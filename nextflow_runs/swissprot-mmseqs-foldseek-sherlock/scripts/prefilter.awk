# Keep the hits whose Swiss-Prot target is in the reference list (file 1,
# accession in column 1), and add the target accession and the hit's rank
# among that query's hits (1 = best; the search tools write each query's hits
# best first).
#   awk -F'\t' -f prefilter.awk reference.tsv hits.tsv
FNR == NR { if (FNR > 1) ref[$1] = 1; next }
FNR == 1  { print $0 "\ttarget_accession\thit_rank"; next }
{
    if ($1 != prev) { rank = 0; prev = $1 }
    rank++
    acc = $2
    if (acc ~ /^AF-/)     { split(acc, a, "-"); acc = a[2] }
    else if (acc ~ /\|/)  { split(acc, a, "|"); acc = a[2] }
    if (acc in ref) print $0 "\t" acc "\t" rank
}
