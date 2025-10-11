#!/usr/bin/env bash

# assumes wrapper script was already run for 1000G

# onekg g2f dir
cd ../2025_09-onekg_giggle2fusion/g2f

# tsv to bedpe
./fusiontsv2bedpe.py \
    -i pop_normal_fusions_pe_and_sample.tsv \
    -b ..//2025_09-gene_bed/grch37.bed \
    -o pop_normal_fusion_pe_sample.bedpe

# reorder cols
awk -F'\t' \
    '{OFS="\t"; print $5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$3,$4}' \
    pop_normal_fusion_pe_sample.bedpe > z && mv z pop_normal_fusion_pe_sample.bedpe

### cleanup 

# rm trivial self hits
awk '$5 != $10' pop_normal_fusion_pe_sample_sorted.bedpe \
    > pop_normal_fusion_pe_sample_sorted_no_dups1.bedpe

# rm non trivial duplicates a-b, b-a
awk -F'\t' \
    '!seen[($5 <= $10) ? $5 FS $10 : $10 FS $5]++' \
    pop_normal_fusion_pe_sample_sorted_no_dups1.bedpe > pop_normal_fusion_pe_sample_sorted_no_dups2.bedpe
    
# sanity check sort
sort pop_normal_fusion_pe_sample_sorted_no_dups2.bedpe \
        --parallel 8 \
        --buffer-size=62G \
        -k1,1V -k2,2n -k3,3n -k6,6V -k7,7n -k8,8n \
        -T /data/jake/tmp \
        > pop_normal_fusion_pe_sample_final.bedpe
# bgzip 
bgzip pop_normal_fusion_pe_sample_final.bedpe
bgzip -dc pop_normal_fusion_pe_sample_final.bedpe.gz > decompressed.tsv

# tabix index
#tabix -p bed pop_normal_fusion_pe_sample_final.bedpe.gz &> tabix_pop_normal_fusion.log

# error: looks like a newline near end of file
# i use tail -n 10000 | less -SN
# the newline occurs at 5520 in the last 10k lines, so 10k-5519=4481 from the end
# the file has 608830293 lines, so 608830293-4481=608825812
# so take head -n 608825812
head -n 608825812 decompressed.tsv | bgzip > pop_normal_fusion_pe_sample_final2.bedpe.gz

# tabix index
tabix -p bed pop_normal_fusion_pe_sample_final2.bedpe.gz &> tabix_pop_normal_fusion2.log

# rename
mv pop_normal_fusion_pe_sample_final2.bedpe.gz onekg_fusion.bedpe.gz
mv pop_normal_fusion_pe_sample_final2.bedpe.gz.tbi onekg_fusion.bedpe.gz.tbi

echo "upload to zenodo"
