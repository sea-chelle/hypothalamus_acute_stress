#!/bin/bash -l

#$ -N atac_aggregate
#$ -l h_rt=02:00:00
#$ -l mem=8G
#$ -pe smp 4
#$ -cwd
#$ -j y
#$ -o logs/

set -euo pipefail

# Load project paths
source "$HOME/Scratch/ATACseq/scripts/config.sh"

# Activate conda environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atacseq

POOL_DIR="$PROJECT_DIR/TOBIAS/pooled"
BIN_DIR="$PROJECT_DIR/TOBIAS/BINDetect"
OUTDIR="$PROJECT_DIR/TOBIAS/aggregate"

mkdir -p "$OUTDIR"

MOTIFS=(
    "Atf3_MA1988.2"
    "Irf1_MA0050.4"
    "AhrArnt_MA0006.2"
)

for MOTIF in "${MOTIFS[@]}"
do
    echo "Generating aggregate footprint for $MOTIF"

    # Use the same motif sites across all four conditions
    TFBS="$BIN_DIR/male/$MOTIF/beds/${MOTIF}_all.bed"

    TOBIAS PlotAggregate \
        --TFBS "$TFBS" \
        --signals \
            "$POOL_DIR/male_basal/male_basal_corrected.bw" \
            "$POOL_DIR/male_60/male_60_corrected.bw" \
            "$POOL_DIR/female_basal/female_basal_corrected.bw" \
            "$POOL_DIR/female_60/female_60_corrected.bw" \
        --signal-labels \
            Male_Basal Male_60 Female_Basal Female_60 \
        --title "$MOTIF" \
        --flank 60 \
        --share-y signals \
        --output "$OUTDIR/${MOTIF}_aggregate.pdf" \
        --output-txt "$OUTDIR/${MOTIF}_aggregate.txt"

    echo "Finished $MOTIF"
done

echo "Finished all aggregate footprint plots"
