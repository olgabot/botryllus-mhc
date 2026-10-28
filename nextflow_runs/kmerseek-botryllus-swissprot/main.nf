#!/usr/bin/env nextflow
// kmerseek search of the whole Botryllus proteome against UniProtKB/Swiss-Prot, with each
// of kmerseek's 19 alphabets at the k-mer sizes that carry the bits of protein20 at k = 6
// and at k = 7 (assets/pairs.tsv, written by bin/plan_pairs.py): 38 alphabet-ksize pairs.
//
// For each pair: one Swiss-Prot index with an E-value fit, then one search per chunk of
// 2,000 Botryllus proteins. Matches with E <= 10 are kept. A pair whose E-value fit fails
// even on 1,000 calibration sequences is searched without extension, and its matches are
// kept only below a Bonferroni-corrected region p-value (see params.exact_max_tail_probability).
//
// Run on Sherlock only: see the Makefile.

nextflow.enable.dsl = 2

// Per-pair first asks for the search, from the 0.4 dark-set run (assets/search_resources.tsv).
def searchResources() {
    def f = file(params.search_resources)
    def rows = f.splitCsv(header: true, sep: '\t')
    rows.collectEntries { r -> ["${r.alphabet}.k${r.ksize}".toString(), [r.memory_gb as int, r.time_h as int]] }
}

process PREPARE_QUERY {
    tag "botryllus"
    publishDir "${params.outdir}/query", mode: 'copy', pattern: '*.tsv'
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    path fasta

    output:
    path "chunks/chunk_*.fasta", emit: chunks
    path "query_renaming.tsv", emit: renaming

    script:
    """
    prepare_query.py ${fasta} ${params.chunk_size} chunks
    mv chunks/query_renaming.tsv .
    """
}

process PREPARE_SWISSPROT {
    tag "swissprot"
    storeDir "${params.index_dir}/swissprot"
    cpus 1
    memory '4 GB'
    time '1h'

    input:
    path fasta_gz

    output:
    tuple path("swissprot.accession_headers.fasta"), path("swissprot.metadata.tsv")

    script:
    """
    prepare_swissprot.py ${fasta_gz} swissprot.accession_headers.fasta swissprot.metadata.tsv
    """
}

process KMERSEEK_INDEX {
    tag "${alphabet}.k${ksize}"
    // Indexes do not depend on the query, so they survive a new work dir or a rerun.
    storeDir "${params.index_dir}/kmerseek"
    cpus 8
    memory { 128.GB * task.attempt }
    time   { 6.h * task.attempt }

    input:
    tuple val(alphabet), val(ksize), val(protein20_k), val(penalty), val(xdrop), path(reference)

    output:
    tuple val(alphabet), val(ksize), val(protein20_k), val(penalty), val(xdrop),
          path("swissprot_2026_03.${alphabet}.k${ksize}.lcremoved.kmerseek.rocksdb")

    script:
    def idx = "swissprot_2026_03.${alphabet}.k${ksize}.lcremoved.kmerseek.rocksdb"
    """
    set -euo pipefail
    export RAYON_NUM_THREADS=${task.cpus}
    kmerseek --version
    kmerseek index --input ${reference} --output ${idx} \\
        --alphabet ${alphabet} --ksize ${ksize} --scaled 1 --remove-low-complexity \\
        --extend-mismatch-penalty ${penalty} --extend-xdrop ${xdrop} \\
        --ka-queries ${params.ka_queries} --ka-seed 1 --ka-survival-out ka_survival.csv \\
        --kmer-stats-out spectrum.csv.gz 2>&1 | tee index.log

    # The `fitted` column of kmerseek's survival CSV, found by name.
    fitted() { awk -F, 'NR==1 { for (i=1;i<=NF;i++) if (\$i=="fitted") c=i; next }
                        NR==2 && c { print tolower(\$c) }' "\$1"; }
    status=\$(fitted ka_survival.csv)
    if [ "\$status" != true ]; then
        # Too few scored regions on ${params.ka_queries} sequences: refit on more, as
        # notebook 235 did. The fit is stored in the index.
        kmerseek calibrate --target ${idx} \\
            --extend-mismatch-penalty ${penalty} --extend-xdrop ${xdrop} \\
            --ka-queries ${params.ka_queries_retry} --ka-seed 1 \\
            --ka-survival-out ka_survival.retry.csv 2>&1 | tee calibrate.log
        status=\$(fitted ka_survival.retry.csv)
    fi
    [ "\$status" = true ] && echo fitted > fit_status.txt || echo unfitted > fit_status.txt
    mv index.log ka_survival*.csv spectrum.csv.gz fit_status.txt ${idx}/
    [ -f calibrate.log ] && mv calibrate.log ${idx}/ || true
    # Searches open the index read-only; read-only files make any write fail loudly.
    # The folder itself stays writable so storeDir can move it into place.
    chmod -R a-w ${idx}
    chmod u+w ${idx}
    """

    stub:
    def idx = "swissprot_2026_03.${alphabet}.k${ksize}.lcremoved.kmerseek.rocksdb"
    """
    mkdir ${idx}
    touch ${idx}/CURRENT
    echo fitted > ${idx}/fit_status.txt
    """
}

process KMERSEEK_SEARCH {
    tag "${alphabet}.k${ksize}.${chunk.baseName}"
    publishDir "${params.outdir}/matches/${alphabet}.k${ksize}", mode: 'link', pattern: '*.csv.zst'
    publishDir "${params.outdir}/best/${alphabet}.k${ksize}", mode: 'copy', pattern: '*.best.tsv'
    publishDir "${params.outdir}/logs/${alphabet}.k${ksize}", mode: 'copy', pattern: '*.log'
    cpus 4
    memory { 1.GB * (mem_gb as int) * task.attempt }
    time   { 1.h * (time_h as int) * task.attempt }

    input:
    tuple val(alphabet), val(ksize), val(protein20_k), val(penalty), val(xdrop), path(index),
          val(mem_gb), val(time_h), path(chunk)

    output:
    tuple val(alphabet), val(ksize), val(protein20_k), path("${slug}.best.tsv"), emit: best
    path "${slug}.matches.csv.zst", emit: matches
    path "${slug}.log", emit: log

    script:
    slug = "${chunk.baseName}.${alphabet}.k${ksize}"
    """
    set -uo pipefail
    export RAYON_NUM_THREADS=${task.cpus}
    # A fitted pair is extended and keeps E <= ${params.max_evalue}. An unfitted pair runs
    # without extension, where every E-value is inf, and keeps region p-values below
    # ${params.exact_max_tail_probability}. awk reads "inf" as 0, so inf and nan are dropped
    # by name before the comparison.
    if grep -qx fitted ${index}/fit_status.txt; then
        statistic=region_evalue; cutoff=${params.max_evalue}
        flags="--extend-mismatch-penalty ${penalty} --extend-xdrop ${xdrop} --ka-queries 0"
    else
        statistic=region_tail_probability; cutoff=${params.exact_max_tail_probability}
        flags="--extend-mismatch-penalty 0"
    fi
    # Nextflow runs this with bash -e; the pipeline's exit codes are read below instead.
    set +e
    kmerseek search --query ${chunk} --target ${index} \\
        --alphabet ${alphabet} --ksize ${ksize} \$flags \\
        --max-query-pvalue ${params.max_query_pvalue} --min-region-score ${params.min_region_score} \\
        2> ${slug}.log \\
      | awk -F, -v col="\$statistic" -v max="\$cutoff" '
            NR==1 { for (i=1;i<=NF;i++) if (\$i==col) c=i
                    if (!c) { print col " column missing" > "/dev/stderr"; exit 3 }
                    print; next }
            tolower(\$c) ~ /inf|nan/ { next }
            (\$c+0) <= (max+0)' \\
      | zstd -q -T2 -o ${slug}.matches.csv.zst
    status=(\${PIPESTATUS[@]})
    set -e
    [ "\${status[0]}" -eq 0 ] && [ "\${status[1]}" -eq 0 ] && [ "\${status[2]}" -eq 0 ] \\
        || { echo "search \${status[0]}, filter \${status[1]}, zstd \${status[2]}" >&2; tail -20 ${slug}.log >&2; exit 1; }

    # One row per Botryllus protein with any kept match: how many Swiss-Prot proteins it
    # matched, how many at ${params.strong_cutoff} or below, and its best match. Columns
    # are found by header name.
    printf 'query_name\\tstatistic\\tn_targets\\tn_targets_strong\\tbest_target\\tbest_value\\n' > ${slug}.best.tsv
    zstd -dc ${slug}.matches.csv.zst \\
    | awk -F, -v col="\$statistic" -v strong=${params.strong_cutoff} '
        NR==1 { for (i=1;i<=NF;i++) { if (\$i=="query_name") q=i; if (\$i=="target_name") t=i; if (\$i==col) c=i }
                next }
        { v = \$c + 0; key = \$q SUBSEP \$t
          if (!(key in seen)) { seen[key] = 1; n[\$q]++ }
          if (v <= strong && !(key in sseen)) { sseen[key] = 1; s[\$q]++ }
          if (!(\$q in best) || v < best[\$q]) { best[\$q] = v; bt[\$q] = \$t } }
        END { for (x in best) printf "%s\\t%s\\t%d\\t%d\\t%s\\t%g\\n", x, col, n[x], s[x] + 0, bt[x], best[x] }' \\
    | sort >> ${slug}.best.tsv
    """

    stub:
    slug = "${chunk.baseName}.${alphabet}.k${ksize}"
    """
    printf 'query_name,target_name,region_evalue\\n' | zstd -q -o ${slug}.matches.csv.zst
    printf 'query_name\\tstatistic\\tn_targets\\tn_targets_strong\\tbest_target\\tbest_value\\nFUN000001_FUN000001\\tregion_evalue\\t1\\t1\\tP12345\\t1e-20\\n' > ${slug}.best.tsv
    touch ${slug}.log
    """
}

process COLLECT {
    tag "summary"
    publishDir "${params.outdir}/summary", mode: 'copy'
    cpus 1
    memory '16 GB'
    time '2h'

    input:
    path pairs_tsv
    path "best/*"
    path index_status
    tuple path(reference), path(metadata)
    path query_chunks

    output:
    path "kmerseek_botryllus_swissprot.*"

    script:
    """
    collect_summary.py --pairs ${pairs_tsv} --best-dir best --fit-status ${index_status} \\
        --metadata ${metadata} --queries ${query_chunks} --strong ${params.strong_cutoff} \\
        --out-prefix kmerseek_botryllus_swissprot
    """
}

workflow {
    def resources = searchResources()
    def pairs = Channel.fromPath(params.pairs)
        .splitCsv(header: true, sep: '\t')
        .map { r -> [r.alphabet, r.ksize, r.protein20_k, r.penalty, r.xdrop] }

    PREPARE_QUERY(file(params.query, checkIfExists: true))
    def reference = PREPARE_SWISSPROT(file(params.swissprot_fasta, checkIfExists: true))

    def indexes = KMERSEEK_INDEX(pairs.combine(reference.map { fa, meta -> fa }))

    def jobs = indexes
        .map { a, k, p20, pen, xd, idx ->
            def r = resources["${a}.k${k}".toString()]
            if (r == null) error "no row for ${a} k${k} in ${params.search_resources}"
            [a, k, p20, pen, xd, idx, r[0], r[1]] }
        .combine(PREPARE_QUERY.out.chunks.flatten())
    KMERSEEK_SEARCH(jobs)

    // One line per pair: alphabet, ksize, fitted or unfitted.
    def fitStatus = indexes
        .map { a, k, p20, pen, xd, idx -> "${a}\t${k}\t${idx.resolve('fit_status.txt').text.trim()}" }
        .collectFile(name: 'fit_status.tsv', newLine: true, sort: true)

    COLLECT(file(params.pairs), KMERSEEK_SEARCH.out.best.map { it[3] }.collect(), fitStatus,
            reference, PREPARE_QUERY.out.chunks.collect())
}
