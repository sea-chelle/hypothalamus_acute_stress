#!/bin/bash -l

#$ -N h3k27ac_consensus
#$ -l h_rt=04:00:00
#$ -l mem=4G
#$ -pe smp 1
#$ -cwd
#$ -j y
#$ -o logs/

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atacseq

H3K27AC_DIR="$HOME/Scratch/ATACseq/260904_PRJ26157_mouse_ARS/ferris_h3k27ac"

ALL_PEAKS="$H3K27AC_DIR/all_h3k27ac_peaks_with_samples.bed"
MERGED_PEAKS="$H3K27AC_DIR/h3k27ac_merged_peaks.bed"
CONSENSUS_PEAKS="$H3K27AC_DIR/h3k27ac_consensus_2of3.bed"

rm -f "$ALL_PEAKS"

for PEAK_FILE in "$H3K27AC_DIR"/*.narrowPeak.gz; do

    SAMPLE=$(basename "$PEAK_FILE" .narrowPeak.gz)

    echo "Adding peaks from: $SAMPLE"

    zcat "$PEAK_FILE" | \
    awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"} {
        print $1, $2, $3, sample
    }' >> "$ALL_PEAKS"

done

sort -k1,1 -k2,2n "$ALL_PEAKS" \
    > "$H3K27AC_DIR/all_h3k27ac_peaks_sorted.bed"

bedtools merge \
    -i "$H3K27AC_DIR/all_h3k27ac_peaks_sorted.bed" \
    -c 4,4 \
    -o distinct,count_distinct \
    > "$MERGED_PEAKS"

awk 'BEGIN{OFS="\t"}
     $5 >= 2 {
         print $1, $2, $3, $4, $5
     }' "$MERGED_PEAKS" \
    > "$CONSENSUS_PEAKS"

echo "H3K27ac consensus complete"

echo "Total merged regions:"
wc -l "$MERGED_PEAKS"

echo "Regions supported by >=2 replicates:"
wc -l "$CONSENSUS_PEAKS"
