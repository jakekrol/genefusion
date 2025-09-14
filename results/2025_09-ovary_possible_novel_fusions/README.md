# Ovary possible novel fusions

- Top scoring fusions scored by logreg model

```
# start w/ top 1000
head -n 1000 ovary_rand_inferred_sort.tsv > rand_inferred_top1000.tsv
cut -f 1,2 rand_inferred_top1000.tsv > z && mv z rand_inferred_top1000.tsv
sed -i 's|\t|--|' rand_inferred_top1000.tsv
# filter using FusionAnnotator
FusionAnnotator --annotate rand_inferred_top1000.tsv --genome_lib_dir $GENOME_LIB_DIR --max_neighbor_dist 100000 > rand_inferred_top1000_annotated.txt
grep -v "NEIGHBOR" rand_inferred_top1000_annotated.txt | grep -v "BLASTPAIR" | grep -v "PARALOG" > rand_inferred_top1000_annotated_filt.txt
# clean up duplicates (caused by neighbors)
script='import pandas as pd;f="./rand_inferred_top1000_annotated_filt.txt";df=pd.read_csv(f,sep="\t");df["x"]=df["left--right"].apply(lambda x: set(x.split("--")));df.drop_duplicates(subset=["x"],inplace=True);df.drop(columns=["x"],inplace=True);df.to_csv("rand_inferred_top1000_annotated_filt_no_dups.txt",sep="\t",index=False)'
python -c "$script"
```