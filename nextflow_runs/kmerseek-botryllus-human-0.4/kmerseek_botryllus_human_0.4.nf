#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Botryllus schlosseri proteins (query) vs human GENCODE canonical translations
// (target), searched with kmerseek 0.4.0.
//
// What 0.4.0 changes versus the 0.3 pipeline next door:
//   * `--encoding hp-lehninger` became `--alphabet hp_lehninger2`.
//   * Matched regions are grown past their exact k-mer run (ungapped extension,
//     X-drop), so a region can hold positions where query and target disagree.
//     region_n_mismatches counts them.
//   * Every row carries region_ka_bits and region_evalue. The E-value needs a
//     Karlin-Altschul lambda and K, fitted for one (mismatch penalty, X-drop)
//     pair. The fit is done once when the index is built and stored in it; a
//     search that passes a different penalty or X-drop refits on every run, so
//     index and search here pass the same two numbers.
//   * The index refuses to open any pre-0.4 index, so this run cannot reuse the
//     0.3 index directory and writes its own.
params.kmerseek_src      = "${HOME}/code/kmerseek/.claude/worktrees/release-0.4.0"
params.kmerseek_main_src = "${HOME}/code/kmerseek/.claude/worktrees/release-0.4.0"
params.input    = "${HOME}/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa"
params.query    = "${HOME}/data/botryllus/Bs_proteins.fa"
params.outdir   = "${HOME}/data/botryllus/kmerseek-0.4"
params.alphabet = "hp_lehninger2"
params.ksize    = 24

// Mismatch penalty C and X-drop for the extension. C = 1.63 is the optimum for a
// two-letter HP alphabet read off the per-position HP agreement of related pairs
// (2024-kmerseek-analysis notebook 230); X = 4C is the kmerseek default ratio.
params.extend_mismatch_penalty = 1.63
params.extend_xdrop            = 6.52

// Database sequences searched against the index to fit lambda and K. 0 skips the
// fit, and then every region_evalue comes back infinite.
params.ka_queries = 200
params.ka_seed    = 1

params.max_query_pvalue  = 0.05
params.min_region_score  = 1.3
params.visualize_query_name = "BHF"
params.visualize_max_hits   = 20

// Naming stem shared by every published file.
def stem = "kmerseek-botryllus-vs-human-gencode.canonical.v49.${params.alphabet}.k${params.ksize}.v0.4.0"

process KMERSEEK_BUILD {
    tag "build:0.4.0"

    output:
    path "kmerseek", emit: binary

    script:
    // The anaconda `ld` on PATH cannot link rocksdb, and bindgen needs to be
    // pointed at the macOS SDK explicitly. Dropping anaconda from PATH for this
    // one task is what makes the build work on this laptop.
    """
    echo "[BUILD] Start: \$(date)"
    start=\$(date +%s)

    export PATH=\$(echo "\$PATH" | tr ':' '\\n' | grep -v anaconda3 | paste -sd: -)
    export SDKROOT=\$(xcrun --show-sdk-path)
    export BINDGEN_EXTRA_CLANG_ARGS="-isysroot \$SDKROOT"
    export RUSTFLAGS="-C link-arg=-undefined -C link-arg=dynamic_lookup"

    cargo build --release --manifest-path ${params.kmerseek_src}/Cargo.toml

    cp ${params.kmerseek_src}/target/release/kmerseek .
    ./kmerseek --version

    end=\$(date +%s)
    echo "[BUILD] End: \$(date)"
    echo "[BUILD] Elapsed: \$((end - start)) seconds"
    """
}

process KMERSEEK_INDEX {
    tag "index:${params.alphabet}:k${params.ksize}"

    publishDir params.outdir, mode: 'copy', pattern: "*.ka_survival.csv"

    input:
    path binary
    path input_fa

    output:
    path "*.kmerseek.rocksdb",  emit: index
    path "*.ka_survival.csv",   emit: ka_survival

    script:
    """
    echo "[INDEX] Start: \$(date)"
    start=\$(date +%s)

    ./${binary} index \\
        --input ${input_fa} \\
        --alphabet ${params.alphabet} \\
        --ksize ${params.ksize} \\
        --scaled 1 \\
        --extend-mismatch-penalty ${params.extend_mismatch_penalty} \\
        --extend-xdrop ${params.extend_xdrop} \\
        --ka-queries ${params.ka_queries} \\
        --ka-seed ${params.ka_seed} \\
        --ka-survival-out ${stem}.ka_survival.csv

    end=\$(date +%s)
    echo "[INDEX] End: \$(date)"
    echo "[INDEX] Elapsed: \$((end - start)) seconds"
    """
}

process KMERSEEK_SEARCH {
    tag "search:${params.alphabet}:k${params.ksize}"

    input:
    path binary
    path query
    path target_index

    output:
    path "*.raw.csv", emit: raw_csv

    script:
    // No date in this filename: it is a work-dir intermediate, and a date would
    // change this process's cache hash every midnight and re-run the search on
    // any -resume that crosses one.
    def out_csv = "${stem}.raw.csv"
    """
    echo "[SEARCH] Start: \$(date)"
    start=\$(date +%s)

    ./${binary} search \\
        --query ${query} \\
        --target ${target_index} \\
        --alphabet ${params.alphabet} \\
        --ksize ${params.ksize} \\
        --extend-mismatch-penalty ${params.extend_mismatch_penalty} \\
        --extend-xdrop ${params.extend_xdrop} \\
        --ka-queries ${params.ka_queries} \\
        --ka-seed ${params.ka_seed} \\
        --max-query-pvalue ${params.max_query_pvalue} \\
        --min-region-score ${params.min_region_score} \\
        --output ${out_csv}

    end=\$(date +%s)
    echo "[SEARCH] End: \$(date)"
    echo "[SEARCH] Elapsed: \$((end - start)) seconds"
    """
}

process KMERSEEK_FILTER {
    tag "filter:bh:p0.05"

    publishDir params.outdir, mode: 'copy'

    input:
    path raw_csv

    output:
    path "*.bh_corrected.csv", emit: bh_csv
    path "*.filtered.csv",     emit: filtered_csv

    script:
    def out_stem = raw_csv.baseName.replace('.raw', '')
    // 0.4 renamed the region p-value columns and added region_evalue, so this is
    // not the 0.3 filter script. region_tail_probability is the raw Poisson tail
    // the region score is the -log10 of, which is what Benjamini-Hochberg needs.
    """
    /Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python3 << 'PYEOF'
import polars as pl
import numpy as np
from statsmodels.stats.multitest import multipletests

df = pl.scan_csv("${raw_csv}", ignore_errors=True).collect(engine="streaming")
stem = "${out_stem}"

p = np.clip(df["region_tail_probability"].to_numpy().astype(float), 1e-300, 1.0)
_, p_bh, _, _ = multipletests(p, alpha=0.05, method="fdr_bh")
df = df.with_columns(pl.Series("region_tail_probability_bh", p_bh))
df.write_csv(f"{stem}.bh_corrected.csv")

sig = df.filter(pl.col("region_tail_probability_bh") < 0.05)
sig.write_csv(f"{stem}.filtered.csv")

ev = df["region_evalue"].to_numpy().astype(float)
n_finite = int(np.isfinite(ev).sum())
print(f"Total region rows: {len(df):_}")
print(f"BH-significant (region_tail_probability_bh < 0.05): {len(sig):_}")
print(f"Rows with a finite E-value: {n_finite:_} of {len(df):_}")
if n_finite:
    print(f"Rows with E-value < 1: {int((ev < 1).sum()):_}")
    print(f"Rows with E-value < 0.001: {int((ev < 0.001).sum()):_}")
else:
    print("No finite E-value: the Karlin-Altschul fit was refused, check the ka_survival CSV")
PYEOF
    """
}

process KMERSEEK_EXTRACT_QUERY_CSV {
    tag "extract:${params.visualize_query_name}"

    publishDir params.outdir, mode: 'copy'

    input:
    path raw_csv

    output:
    path "*.query_subset.csv", emit: query_csv

    script:
    def out_csv = "${params.visualize_query_name}.${params.alphabet}.k${params.ksize}.v0.4.0.query_subset.csv"
    // visualize_hits.py resolves --query-name by scanning the whole CSV without
    // the streaming engine, which is impractical on the whole-proteome output.
    """
    /Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python3 << 'PYEOF'
import polars as pl

df = (
    pl.scan_csv("${raw_csv}")
    .filter(pl.col("query_name").str.strip_chars() == "${params.visualize_query_name}")
    .collect(engine="streaming")
)
df.write_csv("${out_csv}")
print(f"Extracted {len(df):_} rows for query '${params.visualize_query_name}'")
PYEOF
    """
}

process KMERSEEK_VISUALIZE {
    tag "visualize:${params.visualize_query_name}"

    publishDir params.outdir, mode: 'copy'

    input:
    path query_csv
    path query_fasta

    output:
    path "*.png"
    path "*.svg"

    script:
    """
    /Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python3 \\
        ${params.kmerseek_main_src}/scripts/visualize_hits.py \\
        --csv ${query_csv} \\
        --query-fasta ${query_fasta} \\
        --query-name ${params.visualize_query_name} \\
        --max-hits ${params.visualize_max_hits} \\
        --output-dir .
    """
}

workflow {
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    query_ch = Channel.fromPath(params.query,  checkIfExists: true)

    KMERSEEK_BUILD()
    KMERSEEK_INDEX(KMERSEEK_BUILD.out.binary, input_ch)
    KMERSEEK_SEARCH(KMERSEEK_BUILD.out.binary, query_ch, KMERSEEK_INDEX.out.index)
    KMERSEEK_FILTER(KMERSEEK_SEARCH.out.raw_csv)
    KMERSEEK_EXTRACT_QUERY_CSV(KMERSEEK_SEARCH.out.raw_csv)
    KMERSEEK_VISUALIZE(KMERSEEK_EXTRACT_QUERY_CSV.out.query_csv, query_ch)
}
