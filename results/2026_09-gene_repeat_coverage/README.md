# goal

quantify the fraction of base pairs overlapping a repeat region for each gene in grch37

# approach

```
bedtools coverage -a gene.bed -b repeat.bed
```