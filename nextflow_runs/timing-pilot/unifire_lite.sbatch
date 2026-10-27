#!/bin/bash
# UniFIRE 2026.3 lite image (UniRule, ARBA and PIRSR rules, no InterProScan inside) on the
# full InterProScan XML from interproscan_full.sbatch. A whole UniFIRE run is this step plus
# that InterProScan run; the full UniFIRE image only bundles the two.
#
# UniFIRE asks for 24 GB at least (2024-kmerseek-analysis proteome-unifire/Makefile).
#
# Usage: sbatch unifire_lite.sbatch <pilot dir>
#SBATCH --job-name=pilot-unifire
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
set -euo pipefail
PILOT=$1
XML=$PILOT/interproscan_full/pilot500_full.xml
OUT=$PILOT/unifire
mkdir -p "$OUT"
date +%s > "$OUT/start_epoch.txt"
apptainer run --bind "$(dirname "$XML")":/volume --bind "$OUT":/out \
    "$PILOT/containers/unifire-2026.3-lite.sif" \
    -i /volume/"$(basename "$XML")" -t iprscanxml -o /out 2>&1 | tee "$OUT/unifire.log"
date +%s > "$OUT/end_epoch.txt"
