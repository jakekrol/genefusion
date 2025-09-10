## goal

Compute fraction of haas fusions that I have intervals for on GRCh37

### purpose

Haas fusions will be quantified for their evidence in 1kg

### approach

- get set of haas fusions
- count how many, after alias mapping, I have interval data for

### expected outcomes

- >80% good to run 1kg experiment
- >50% & <80% still good to run, but will need to add more intervals later
- <50% don't run the experiment yet


### code

```
echo "building hgnc alias map"
./build_hgnc_alias_map.py \
    --hgnc ../../data/2025_09-gene_bed/genenames.tsv \
    --gtf ../../data/2025_09-gene_bed/gencode.v19.annotation.gtf.gz \
    --output hgnc_alias.yaml
# haas
echo "building haas alias map"
./build_haas_alias_map.py \
    --input ../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt \
    --output haas_alias.yaml

# merge
echo "merging hgnc and haas alias maps"
./merge_hgnc_haas_aliases.py \
    --hgnc hgnc_alias.yaml \
    --haas haas_alias.yaml \
    --output merged_alias.yaml

# get set of genes with interval data
echo "getting genes with interval data"
cut -f 4 ../2025_09-gene_bed/grch37.bed | sort > genes_bed.txt

# get haas fusions
echo "getting haas fusions"
tail -n +2 ../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt | \
    cut -f 3 > haas_fusions.txt

# count fraction of haas fusions with interval data
echo "counting fraction of haas fusions with interval data"
./count_frac_haas_fusion_w_interval_data.py --fusions haas_fusions.txt \
    --alias merged_alias.yaml \
    --genes genes_bed.txt \
    --output haas_fusion_interval_fraction.txt
```

### Result
$> 0.96$ fusions from haas I have interval data for both genes.

$>$ because the final counting script cannot handle commas in the haas fusions.

However, this amount is already sufficient to continue with the planned experiment. 
No need to refine.