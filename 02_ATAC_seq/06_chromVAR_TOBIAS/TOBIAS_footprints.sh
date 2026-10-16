#!/bin/bash -l

#$ -N atac_fp_pool
#$ -l h_rt=03:00:00
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

CONDITIONS=(
    "female_basal"
    "female_60"
    "male_basal"
    "male_60"
)

for CONDITION in "${CONDITIONS[@]}"
do
    echo "======================================"
    echo "Running FootprintScores for ${CONDITION}"
    echo "======================================"

    TOBIAS FootprintScores \
        --signal "$POOL_DIR/$CONDITION/${CONDITION}_corrected.bw" \
        --regions "$PROJECT_DIR/DAR_bed/consensus_background.bed" \
        --output "$POOL_DIR/$CONDITION/${CONDITION}_footprints.bw" \
        --cores 4

    echo "Finished FootprintScores for ${CONDITION}"
done

echo "Finished pooled TOBIAS FootprintScores analysis"
