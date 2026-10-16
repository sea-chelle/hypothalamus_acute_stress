#!/bin/bash -l

#$ -N atac_correct_pool
#$ -l h_rt=08:00:00
#$ -l mem=12G
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

CONDITIONS=(
    "female_basal"
    "female_60"
    "male_basal"
    "male_60"
)

for CONDITION in "${CONDITIONS[@]}"
do
    echo "======================================"
    echo "Running ATACorrect for ${CONDITION}"
    echo "======================================"

    OUTDIR="$POOL_DIR/$CONDITION"
    mkdir -p "$OUTDIR"

    TOBIAS ATACorrect \
        --bam "$POOL_DIR/${CONDITION}.bam" \
        --genome "$MM39_FASTA" \
        --peaks "$PROJECT_DIR/DAR_bed/consensus_background.bed" \
        --blacklist "$MM39_BLACKLIST" \
        --outdir "$OUTDIR" \
        --cores 4

    echo "Finished ATACorrect for ${CONDITION}"
done

echo "Finished pooled TOBIAS ATACorrect analysis"
