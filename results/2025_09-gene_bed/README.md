# Goal

A robust system to map GRCh37 gene names to coordinates, and vice versa.

1. Download HGNC gene names
- https://www.genenames.org/download/custom/
- columns: HGNC ID, Alias symbols, Approved symbol, Approved name, Previous symbols, Ensembl gene ID
- chr: 1-23, X, Y

2. Get coordinates from GENCODE
`wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz`

3. Build a BED with coords for all gene names and aliases

```
../../scripts/python/build_gene_bed.py \
    --hgnc ../../data/2025_09-gene_bed/genenames.tsv \
    --gtf ../../data/2025_09-gene_bed/gencode.v19.annotation.gtf.gz \
    --output grch37.chr.bed
```
