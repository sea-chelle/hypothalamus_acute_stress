#!/bin/bash -l

#$ -N atac_homer
#$ -l h_rt=06:00:00
#$ -l mem=8G
#$ -pe smp 4
#$ -cwd
#$ -j y
#$ -o logs/

source "$HOME/Scratch/ATACseq/scripts/config.sh"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atacseq

# force HOMER to use Conda's newer C++ libraries
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# input and output directories
DAR_DIR="$PROJECT_DIR/DAR_bed"
HOMER_DIR="$PROJECT_DIR/HOMER_motifs"

mkdir -p "$HOMER_DIR"

# consensus peak background
BACKGROUND="$DAR_DIR/consensus_background.bed"

# overall UP DARs
findMotifsGenome.pl \
    "$DAR_DIR/overall_UP_DARs.bed" \
    mm39 \
    "$HOMER_DIR/overall_UP" \
    -bg "$BACKGROUND" \
    -size given \
    -p 4

# male UP DARs
findMotifsGenome.pl \
    "$DAR_DIR/male_UP_DARs.bed" \
    mm39 \
    "$HOMER_DIR/male_UP" \
    -bg "$BACKGROUND" \
    -size given \
    -p 4

# female UP DARs
findMotifsGenome.pl \
    "$DAR_DIR/female_UP_DARs.bed" \
    mm39 \
    "$HOMER_DIR/female_UP" \
    -bg "$BACKGROUND" \
    -size given \
    -p 4

# overall DOWN DARs
findMotifsGenome.pl \
    "$DAR_DIR/overall_DOWN_DARs.bed" \
    mm39 \
    "$HOMER_DIR/overall_DOWN" \
    -bg "$BACKGROUND" \
    -size given \
    -p 4

echo "Finished HOMER motif analysis"
