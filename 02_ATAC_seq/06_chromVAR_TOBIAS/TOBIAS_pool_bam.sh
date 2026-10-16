#!/bin/bash -l

#$ -N atac_pool_bams
#$ -l h_rt=02:00:00
#$ -l mem=8G
#$ -pe smp 4
#$ -cwd
#$ -j y
#$ -o logs/

# load project paths
source "$HOME/Scratch/ATACseq/scripts/config.sh"

# activate environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atacseq

OUTDIR="$PROJECT_DIR/TOBIAS/pooled"
mkdir -p "$OUTDIR"

# Female basal: n = 4
samtools merge -@ 4 \
    "$OUTDIR/female_basal.bam" \
    "$BAM_DIR/2128-f-basal_S5_L001.final.bam" \
    "$BAM_DIR/2226-f-basal_S9_L001.final.bam" \
    "$BAM_DIR/2230-f-basal_S10_L001.final.bam" \
    "$BAM_DIR/2471-f-basal_S17_L001.final.bam"

# Female 60 min: n = 2
samtools merge -@ 4 \
    "$OUTDIR/female_60.bam" \
    "$BAM_DIR/2132-f-60_S7_L001.final.bam" \
    "$BAM_DIR/2231-f-60_S11_L001.final.bam"

# Male basal: n = 3
samtools merge -@ 4 \
    "$OUTDIR/male_basal.bam" \
    "$BAM_DIR/2243-m-basal_S12_L001.final.bam" \
    "$BAM_DIR/2249-m-basal_S14_L001.final.bam" \
    "$BAM_DIR/2450-m-basal_S15_L001.final.bam"

# Male 60 min: n = 3
samtools merge -@ 4 \
    "$OUTDIR/male_60.bam" \
    "$BAM_DIR/2224-m-60_S8_L001.final.bam" \
    "$BAM_DIR/2244-m-60_S13_L001.final.bam" \
    "$BAM_DIR/2455-m-60_S16_L001.final.bam"

# index pooled BAMs
for BAM in "$OUTDIR"/*.bam
do
    samtools index -@ 4 "$BAM"
done

echo "Finished pooling and indexing TOBIAS BAMs"
