#!/usr/bin/env bash

# downsample fastqs
conda activate polymerization
./sample_fastqs.sh

# call fusions with star
conda activate star_fusion
./star_fusion.sh

# convert tsv to bedpe
conda activate polymerization
dirs=$(mapfile -t < <(ls -d star_fusion_outdir_*))
for dir in "${dirs[@]}"; do
    fraction="${dir##*_}"
    tsv="${dir}/star-fusion.fusion_predictions.tsv"
    bedpe="k562-${fraction}-star_fusion.bed.gz"
    ./starfusion2bedpe.py --input "${tsv}" --output "${bedpe}" --bgzip --pad
done
