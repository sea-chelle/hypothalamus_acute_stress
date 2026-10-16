#!/bin/bash -l

#$ -N bindetect_m
#$ -l h_rt=06:00:00
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
OUTDIR="$PROJECT_DIR/TOBIAS/BINDetect/male"

mkdir -p "$OUTDIR"

TOBIAS BINDetect \
    --signals \
        "$POOL_DIR/male_basal/male_basal_footprints.bw" \
        "$POOL_DIR/male_60/male_60_footprints.bw" \
    --cond-names male_basal male_60 \
    --motifs "$PROJECT_DIR/TOBIAS/JASPAR2024_mouse_CORE.jaspar" \
    --genome "$MM39_FASTA" \
    --peaks "$PROJECT_DIR/DAR_bed/consensus_background.bed" \
    --outdir "$OUTDIR" \
    --prefix male_basal_vs_60 \
    --cores 4

echo "Finished male BINDetect analysis"
