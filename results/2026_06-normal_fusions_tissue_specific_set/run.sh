#!/usr/bin/env bash

bed='../2025_04-gene_bedfile_cln/grch37.genes.bed'

./run.py

./check_bed.py

# substitute some aliases
sed -i 's|FLVCR1-AS1|FLVCR1|g' recurrent_normal_fusion_set.tsv
sed -i 's|LINC00694|C3ORF86|g' recurrent_normal_fusion_set.tsv
sed -i 's|ZNF674-AS1|ZNF674|g' recurrent_normal_fusion_set.tsv
sed -i 's|MARCH _2|MARCHF2|g' recurrent_normal_fusion_set.tsv
sed -i 's|SEPT7P2|SEPTIN7P2|g' recurrent_normal_fusion_set.tsv

# filter out a few fusions with no coords
grep -v -E "IGHV3OR16-6|hsa-mir-6723|RNA28S5|RNA45S5" ./recurrent_normal_fusion_set.tsv > recurrent_normal_fusion_set.filtered.tsv
tail -n +2 recurrent_normal_fusion_set.filtered.tsv | cut -f1,2 > recurrent_normal_fusion_set.filtered.nohead.c1c2.tsv


# get combined bed file
cat $bed missing_genes.bed | sort -k1,1 -k2,2n -k3,3n > grch37.genes.added.bed

# test
./test_merge_bed_fusion_set.py && echo "# successful merge of fusion set and bed"
