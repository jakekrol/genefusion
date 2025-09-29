#!/usr/bin/env bash

# setup haas for join with 1kg results
./setup_haas.py \
    --haas ../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt \
    --bedfile ../results/2025_09-gene_bed/grch37.bed \
    --alias_map ../2025_09-haas_interval_frac/merged_alias.yaml \
    --output 'haas_left_right_mapped.tsv'

# rm unresolved gene names
grep -v "^-1" haas_left_right_mapped.tsv > haas_left_right_mapped_cln.tsv

# join haas with 1kg results
onekg="../2025_09-onekg_giggle2fusion/g2f/pop_normal_fusions.tsv"

join.py \
    -x haas_left_right_mapped_cln.tsv \
    -y $onekg \
    -t left \
    -o haas_onekg.tsv \
    -k left,right

# duplicates will have the same 1000 genomes evidence
# however, they may have different haas evidence in the junction (J) or spanning (S) reads
# future analysis may want to compare S reads with onekg reads
drop_duplicates.py -i haas_onekg.tsv -o haas_onekg_no_dups.tsv -c 5

# rm previous FusionAnnotator column
cut -f10 --complement haas_onekg_no_dups.tsv > z && mv z haas_onekg_no_dups.tsv

# fill missing 1000 genomes records with 0
script='import pandas as pd;df=pd.read_csv("haas_onekg_no_dups.tsv",sep="\t");\
df["pe_count_normal"].fillna(0,inplace=True);\
df["pe_count_normal"] = df["pe_count_normal"].astype(int);\
df = df.sort_values(by="pe_count_normal",ascending=False);\
cols = df.columns.tolist();\
cols.insert(2, cols.pop(cols.index("pe_count_normal")));\
df = df[cols];\
df.to_csv("haas_onekg_no_dups_filled.tsv",sep="\t",index=False)'
python3 -c "$script"

# re-run FusionAnnotator after resolving name aliases
tail -n +2 haas_onekg_no_dups_filled.tsv | \
    cut -f1,2 | \
    sed 's|\t|--|' > tmp.txt
FusionAnnotator \
    --annotate tmp.txt \
    --no_add_header_column \
    --genome_lib_dir $GENOME_LIB_DIR \
    > tmp2.txt

# reformat for join
sed 's|--|\t|' tmp2.txt > tmp3.txt
sed -i '1ileft\tright\tannotation' tmp3.txt

# join new annotation
join.py \
    -x haas_onekg_no_dups_filled.tsv \
    -y tmp3.txt \
    -t left \
    -o haas_onekg_no_dups_filled_annotated.tsv \
    -k left,right

# mv annotation column 
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$11,$4,$5,$6,$7,$8,$9,$10}' haas_onekg_no_dups_filled_annotated.tsv > z && mv z haas_onekg_no_dups_filled_annotated.tsv

# cleanup
rm tmp.txt tmp2.txt tmp3.txt

# filter out selfies,overlap,homology,neighbor
grep -v -E "SELFIE|OVERLAP|BLAST|NEIGHBOR" haas_onekg_no_dups_filled_annotated.tsv \
    > haas_onekg_no_dups_filled_annotated_filt.tsv

# inspect filtered fusions
grep -E "SELFIE|OVERLAP|BLAST|NEIGHBOR" haas_onekg_no_dups_filled_annotated.tsv > haas_onekg_filtered_out.tsv

# print the count of filtered fusions
echo "Count of filtered fusions:"
grep -c -E "SELFIE|OVERLAP|BLAST|NEIGHBOR" haas_onekg_no_dups_filled_annotated.tsv 

# histogram of 1kg supporting reads for haas fusions
tail -n +2 haas_onekg_no_dups_filled_annotated.tsv \
    | cut -f3 \
    | hist.py -o haas_onekg_readcount_histogram.png \
    --ylog -x "1000 Genomes depth" -y Frequency \
    --title "All benchmark fusions"

# with filtered data
tail -n +2 haas_onekg_no_dups_filled_annotated_filt.tsv \
    | cut -f3 \
    | hist.py -o haas_onekg_readcount_histogram_filt.png \
    --ylog -x "1000 Genomes depth" -y Frequency \
    --title "Post-filter fusions"

# of filtered fusions
tail -n +2 haas_onekg_annotation_filtered_out.tsv \
    | cut -f3 \
    | hist.py -o haas_onekg_readcount_histogram_filtered_out.png \
    --ylog -x "1000 Genomes depth" -y Frequency \
    --title "Annotation filtered out fusions"

# compute num. fusions with >1000 supporting reads in 1kg
# in filtered data
echo "Count of filtered fusions with >100 supporting reads in 1kg:"
x=$(tail -n +2 haas_onekg_no_dups_filled_annotated_filt.tsv \
    | awk '$3 > 100' | wc -l)
n=$(tail -n +2 haas_onekg_no_dups_filled_annotated_filt.tsv | wc -l)
echo "$x out of $n"

# # inspect fusions with >10000 supporting reads in 1kg
# tail -n +2 haas_onekg_no_dups_filled_annotated.tsv \
#     | awk '$3 > 10000' | \
#     # nr numeric reverse
#     sort -k3,3nr > haas_onekg_depth_gt10000.tsv

# # do the same with >1000 supporting reads in 1kg
# tail -n +2 haas_onekg_no_dups_filled_annotated.tsv \
#     | awk '$3 > 1000' | \
#     sort -k3,3nr > haas_onekg_depth_gt1000.tsv

# careful! we are using the no_dups file here
# scatter plot of haas spanning reads vs 1kg depth
# for cases where haas spanning reads is greater than 10
tail -n +2 haas_onekg_no_dups_filled_annotated.tsv | 
    # awk '$7 > 10' | \
    cut -f 7,11 | \
    scatter.py -o haas_onekg_spanning_vs_1kgdepth.png \
    --xlog --ylog -x "Haas spanning reads" -y "1000 Genomes depth"


