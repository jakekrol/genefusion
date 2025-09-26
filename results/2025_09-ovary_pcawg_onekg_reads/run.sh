#!/usr/bin/env bash

# gather ovary gold standard file
# cp $GENEFUSION/logreg/ovary/ovary_gold_standard_pcawg.tsv .

# use alias mapped file to convert gene any gene names
alia=../2025_09-haas_interval_frac/merged_alias.yaml
sourc=ovary_gold_standard_pcawg.tsv
target=../2025_09-gene_bed/grch37.bed
output=ovary_gold_standard_pcawg_mapped.tsv

../../scripts/python/gene2alias.py \
    --alias $alia \
    --source $sourc \
    --target $target \
    --output $output \
    --source_col 1,2 \
    --source_header

# left sort first!
script='import pandas as pd; from genefusion.genefusion import add_left_right_col; df=pd.read_csv("ovary_gold_standard_pcawg_mapped.tsv", sep="\t"); df=add_left_right_col(df,x="left",y="right"); df=df.sort_values(["left","right"]); df.to_csv("ovary_gold_standard_pcawg_mapped_sorted.tsv", sep="\t", index=False)'
python3 -c "$script"

# rm unresolved gene names
grep -v "^-1" ovary_gold_standard_pcawg_mapped_sorted.tsv > z && mv z ovary_gold_standard_pcawg_mapped_sorted.tsv


# join onekg supporting reads
onekg="../2025_09-onekg_giggle2fusion/g2f/pop_normal_fusions.tsv"
join.py \
    -x ovary_gold_standard_pcawg_mapped_sorted.tsv \
    -y $onekg \
    -t left \
    -o ovary_pcawg_onekg.tsv \
    -k left,right

# fill missing with 0
script='import pandas as pd; df=pd.read_csv("ovary_pcawg_onekg.tsv", sep="\t"); df["pe_count_normal"]=df["pe_count_normal"].fillna(0);df = df.sort_values(["pe_count_normal"],ascending=False); df.to_csv("ovary_pcawg_onekg_filled.tsv", sep="\t", index=False)'
python3 -c "$script"

# annotate 
cut -f 1,2 ovary_pcawg_onekg_filled.tsv | tail -n +2 | sed 's|\t|--|' > fusions.tmp
FusionAnnotator --no_add_header_column --genome_lib_dir $GENOME_LIB_DIR --annotate fusions.tmp > fusions.annotated

# combine
sed 's|--|\t|' fusions.annotated > fusions_annotated.tmp
sed '1i left\tright\tannotation' fusions_annotated.tmp > z && mv z fusions_annotated.tmp
join.py -x ovary_pcawg_onekg_filled.tsv -y fusions_annotated.tmp -t left -o ovary_pcawg_onekg_annotated.tsv -k left,right

# cleanup
rm fusions.tmp fusions.annotated fusions_annotated.tmp

# make filtered version
grep -v -E "SELFIE|OVERLAP|NEIGHBORS|BLAST" ovary_pcawg_onekg_annotated.tsv \
    > ovary_pcawg_onekg_annotated_filtered.tsv

# histogram of onekg supporting reads
tail -n +2 ovary_pcawg_onekg_annotated.tsv | \
    cut -f 3 | \
    hist.py -o ovary_pcawg_onekg_all.hist.png --ylog --bins 30 \
    -x "1KG supporting reads" -y "Frequency" --title "Ovary PCAWG fusions"

# histogram filtered
tail -n +2 ovary_pcawg_onekg_annotated_filtered.tsv | \
    cut -f 3 | \
    hist.py -o ovary_pcawg_onekg_filtered.hist.png --ylog --bins 30 \
    -x "1KG supporting reads" -y "Frequency" --title "Ovary PCAWG fusions"

# count num filtered
n=$(wc -l ovary_pcawg_onekg_annotated.tsv)
x=$(wc -l ovary_pcawg_onekg_annotated_filtered.tsv)
echo "Total fusions: $n"
echo "Fusions after filtering: $x"

# count num filtered fusiosn with >10 supporting reads
tail -n +2 ovary_pcawg_onekg_annotated_filtered.tsv | awk '$3>10' | wc -l




