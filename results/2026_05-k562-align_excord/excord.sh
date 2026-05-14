#!/usr/bin/env bash

rna_bam=k562-rna-Aligned.sortedByCoord.out.bam
rna_bed=k562-rna-Aligned.sortedByCoord.out.bed.gz

# run_excord () {
#     local input_bam="$1"
#     local output_bed="$2"
#     local discordant_distance="$3"

#     echo "# running excord on $input_bam with discordant distance $discordant_distance..."
#     cat $input_bam | \
#         excord --discordantdistance "$discordant_distance" /dev/stdin | \
#         bgzip -c > "$output_bed"
# }
# # cat PCAWG.7d1f5c61-5f71-406e-a191-7c4d26896662.STAR.v1.bam | excord --discordantdistance 500 /dev/stdin | bgzip -c > test.bed.gz

# # rna
# run_excord $rna_bam $rna_bed 500
excord-lr \
    --bam $rna_bam \
    --out test.bed \
    --verbose &> excord-lr.log

