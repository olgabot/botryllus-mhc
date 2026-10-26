#!/bin/bash
#SBATCH --job-name=233-build-kmerseek
#SBATCH --partition=hns
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=233-build-kmerseek-%j.out
#
# Build kmerseek v0.4.0 on Sherlock for notebook 233, into $SCRATCH/kmerseek-v0.4.0.
# The module set and the CPATH fix are the ones from notebook 241's sbatch in
# seanome/2024-kmerseek-analysis, where they built kmerseek on Sherlock.
set -euo pipefail
: "${SCRATCH:?}"
SRC=$SCRATCH/kmerseek-v0.4.0
module load rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
# The llvm module puts LLVM's libc++ headers on CPATH, and g++ then compiles RocksDB
# against them and fails ("'abort' has not been declared in 'std'"). Keep gcc's headers
# only, and give bindgen libclang directly.
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
[ -d "$SRC/.git" ] || git clone --quiet https://github.com/seanome/kmerseek "$SRC"
(cd "$SRC" && git fetch --quiet --tags && git checkout --quiet v0.4.0 && cargo build --release --bin kmerseek)
"$SRC/target/release/kmerseek" --version
