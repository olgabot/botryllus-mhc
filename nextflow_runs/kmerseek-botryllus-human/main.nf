#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.kmerseek_src = "${HOME}/code/kmerseek"
params.input    = "${HOME}/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa"
params.query    = "${HOME}/data/botryllus/Bs_proteins.fa"
params.outdir   = "${HOME}/data/botryllus"
params.encoding  = "hp"
params.ksize     = 24
params.timestamp = new Date().format('yyyy-MM-dd')

process KMERSEEK_BUILD {
    tag "build"

    output:
    path "kmerseek-rust", emit: binary

    script:
    """
    echo "[BUILD] Start: \$(date)"
    start=\$(date +%s)

    cargo build --release --manifest-path ${params.kmerseek_src}/Cargo.toml

    cp ${params.kmerseek_src}/target/release/kmerseek-rust .

    end=\$(date +%s)
    echo "[BUILD] End: \$(date)"
    echo "[BUILD] Elapsed: \$((end - start)) seconds"
    """
}

process KMERSEEK_INDEX {
    tag "index:${params.encoding}:k${params.ksize}"

    input:
    path binary
    path input_fa

    output:
    path "*.kmerseek.rocksdb", emit: index

    script:
    """
    echo "[INDEX] Start: \$(date)"
    start=\$(date +%s)

    ./${binary} index \\
        --encoding ${params.encoding} \\
        --ksize ${params.ksize} \\
        --input ${input_fa}

    end=\$(date +%s)
    echo "[INDEX] End: \$(date)"
    echo "[INDEX] Elapsed: \$((end - start)) seconds"
    """
}

process KMERSEEK_SEARCH {
    tag "search:${params.encoding}:k${params.ksize}"

    input:
    path binary
    path query
    path target_index

    output:
    path "*.raw.csv", emit: raw_csv

    script:
    def out_csv = "kmerseek-botryllus-vs-human-gencode.canonical.v49.${params.encoding}.k${params.ksize}.${params.timestamp}.raw.csv"
    """
    echo "[SEARCH] Start: \$(date)"
    start=\$(date +%s)

    ./${binary} search \\
        --query ${query} \\
        --target ${target_index} \\
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
    def stem = raw_csv.baseName.replace('.raw', '')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    from scipy.stats import poisson as _poisson
    from statsmodels.stats.multitest import multipletests

    df = pd.read_csv("${raw_csv}", on_bad_lines="skip")
    stem = "${stem}"

    if "expected_shared_kmers" in df.columns:
        lam = df["expected_shared_kmers"].values.copy().astype(float)
        k   = df["n_intersecting_hashes"].values.astype(int)
        lam = np.where(lam <= 0, 1e-300, lam)
        poisson_p = np.clip(_poisson.sf(k - 1, lam), 1e-300, 1.0)
        _, p_bh, _, _ = multipletests(poisson_p, alpha=0.05, method="fdr_bh")
        df["poisson_p"]    = poisson_p
        df["poisson_p_bh"] = p_bh
    else:
        import warnings
        warnings.warn("expected_shared_kmers not found — p-values set to NaN")
        df["poisson_p"]    = float("nan")
        df["poisson_p_bh"] = float("nan")

    df.to_csv(f"{stem}.bh_corrected.csv", index=False)

    sig = df[df["poisson_p_bh"] < 0.05]
    sig.to_csv(f"{stem}.filtered.csv", index=False)

    print(f"Total hits: {len(df):,}")
    print(f"Significant (poisson_p_bh < 0.05): {len(sig):,}")
    """
}

workflow {
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    query_ch  = Channel.fromPath(params.query,  checkIfExists: true)

    KMERSEEK_BUILD()
    KMERSEEK_INDEX(KMERSEEK_BUILD.out.binary, input_ch)
    KMERSEEK_SEARCH(KMERSEEK_BUILD.out.binary, query_ch, KMERSEEK_INDEX.out.index)
    KMERSEEK_FILTER(KMERSEEK_SEARCH.out.raw_csv)
}
