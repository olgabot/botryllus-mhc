#!/bin/bash
# InterProScan 5.59-91.0 with every member database it ships (no --applications), on the
# pilot FASTA with UniFIRE headers (OX=30301). Two uses:
#   1. the time for a full InterProScan run, which proteinannotator's 6-database subset
#      does not give;
#   2. the XML that UniFIRE's lite image reads. UniFIRE takes the taxid from the sequence
#      ids in this XML, so the input must be the OX= FASTA.
# Same container and data as proteinannotator's INTERPROSCAN task, set up the same way as
# that module's script (copy interproscan.properties, point bin.directory at the image).
#
# Usage: sbatch interproscan_full.sbatch <pilot dir> <ips data dir> <image> <fasta>
#SBATCH --job-name=pilot-ips-full
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
set -euo pipefail
PILOT=$1
IPS_DATA=$2
IMAGE=$3
FASTA=$4
OUT=$PILOT/interproscan_full
mkdir -p "$OUT"
cd "$OUT"
ln -sfn "$IPS_DATA" data
date +%s > start_epoch.txt

TMP=${L_SCRATCH:-/tmp}
apptainer exec --bind "$PILOT" --bind "$IPS_DATA" --bind "$(dirname "$FASTA")" --bind "$TMP" \
    --env TMP="$TMP" "$IMAGE" bash -c '
    set -euo pipefail
    INTERPROSCAN_DIR="$(dirname "$(dirname "$(which interproscan.sh)")")"
    cp "$(find "$INTERPROSCAN_DIR/share" -name interproscan.properties)" .
    sed -i "/^bin\.directory=/ s|.*|bin.directory=$INTERPROSCAN_DIR/bin|" interproscan.properties
    export INTERPROSCAN_CONF=interproscan.properties
    interproscan.sh --cpu '"$SLURM_CPUS_PER_TASK"' --disable-precalc \
        --input '"$FASTA"' --formats XML,TSV --output-file-base pilot500_full \
        --tempdir "$TMP/ips-tmp"
' 2>&1 | tee interproscan.log

date +%s > end_epoch.txt
