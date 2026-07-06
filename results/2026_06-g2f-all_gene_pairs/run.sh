#!/usr/bin/env bash

# get fusion evidence
./run_g2f.py

# aggregate evidence across genes within modality
./run_agg.py

# stack modalitities
./join_modalities.sh

# split fusion evidence into chunks
./chunk.sh

# score
./score_chunks.py

# recollect
./combine_chunks.sh

# remove self-fusions
cd /data/jake/genefusion/results/2026_06-g2f-all_gene_pairs/chunks_scored
files=(blood_scored.tsv.sorted
bone_scored.tsv.sorted
breast_scored.tsv.sorted
esophagus_scored.tsv.sorted
gallbladder_scored.tsv.sorted
headneck_scored.tsv.sorted
kidney_scored.tsv.sorted
liver_scored.tsv.sorted
ovary_scored.tsv.sorted
pancreas_scored.tsv.sorted
prostate_scored.tsv.sorted
)
for f in "${files[@]}"; do
    out="${f}.no_selfie"
    echo "# processing $f to outfile: $out"
    printf "gene_left\tgene_right\tscore_uniform\n" > $out
    tail -n +2 $f | awk '$1 != $2' >> $out
done

# plot tissue-wise score distributions
