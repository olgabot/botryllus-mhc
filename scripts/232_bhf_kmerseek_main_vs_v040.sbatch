#!/bin/bash
#SBATCH --job-name=bhf-kmerseek-main-vs-v040
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=bhf_kmerseek_main_vs_v040_%j.log
# BHF vs human GENCODE v49 canonical translations, searched with two kmerseek builds:
# release-0.4.0 (0aed5ca, the build behind the September results) and main (fd549f9,
# 2026-10-24). Both use the September settings: hp_lehninger2 k=24, scaled 1, mismatch
# penalty 1.63, X-drop 6.52, 200 Karlin-Altschul fit queries, seed 1, and both
# low-complexity arms. The v0.4.0 build is the control: it should give back the 66
# September rows for BHF exactly.
#   sbatch 232_bhf_kmerseek_main_vs_v040.sbatch   # from $SCRATCH/bhf-kmerseek-main-vs-v040
# Expects human.fa and bhf.fa in the working directory. Writes out/<build>.<arm>.csv.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK

WORK=$PWD
mkdir -p out
md5sum human.fa bhf.fa > out/inputs.md5

declare -A COMMIT=([v040]=0aed5ca [main]=fd549f9)
for build in v040 main; do
    src=$WORK/src-$build
    [ -d "$src" ] || git clone -q https://github.com/seanome/kmerseek "$src"
    # Sherlock's system git predates `git -C`.
    (cd "$src" && git fetch -q origin --tags && git checkout -q "${COMMIT[$build]}" \
        && git log -1 --format='%h %ad %s' --date=short) > "out/$build.commit.txt"
    cargo build --release --bin kmerseek --manifest-path "$src/Cargo.toml"
done

for build in v040 main; do
    bin=$WORK/src-$build/target/release/kmerseek
    for arm in lckept lcremoved; do
        lc_flag=""
        [ "$arm" = lcremoved ] && lc_flag="--remove-low-complexity"
        stem=out/$build.$arm
        rm -rf "$stem.rocksdb"
        /usr/bin/time -v -o "$stem.index.time" "$bin" index \
            --input human.fa --output "$stem.rocksdb" \
            --alphabet hp_lehninger2 --ksize 24 --scaled 1 $lc_flag \
            --extend-mismatch-penalty 1.63 --extend-xdrop 6.52 \
            --ka-queries 200 --ka-seed 1 --ka-survival-out "$stem.ka_survival.csv" \
            > "$stem.index.log" 2>&1
        /usr/bin/time -v -o "$stem.search.time" "$bin" search \
            --query bhf.fa --target "$stem.rocksdb" \
            --alphabet hp_lehninger2 --ksize 24 \
            --extend-mismatch-penalty 1.63 --extend-xdrop 6.52 \
            --ka-queries 200 --ka-seed 1 \
            --max-query-pvalue 0.05 --min-region-score 1.3 \
            --output "$stem.csv" > "$stem.search.log" 2>&1
        echo "$build $arm rows: $(($(wc -l < "$stem.csv") - 1))"
        grep -h "Maximum resident\|Elapsed (wall" "$stem.index.time" "$stem.search.time"
    done
done
