# 1000G fusion evidence

Use 1000G fusion evidence for normal fusion filtering.

## Download

https://doi.org/10.5281/zenodo.17317215

```
# download index
wget 'https://zenodo.org/records/17317215/files/onekg_fusion.bedpe.gz.tbi?download=1'
# download table (<4GB)
# this will take some time 
wget 'https://zenodo.org/records/17317215/files/onekg_fusion.bedpe.gz?download=1'
```


## Usage

Gene locations can be found in a provided [bed file](results/2025_09-gene_bed/grch37.bed).

For a fusion of interest, query for the gene with the lesser chromosome number and/or base pair coordinates.

Example: ERG-TMPRSS2 fusion evidence in 1000G.

```
erg_location="21:39751949-40033704"
tabix onekg_fusion.bedpe.gz $erg_location \
    grep ERG | grep TMPRSS2
# 21      39751948        40033704        neg     ERG     21      42836477        42903043        neg     TMPRSS2    8       8
```

Here, 8 samples in 1000G each have 1 supporting read pair for the ERG-TMPRSS2 fusion.

Example 2: LRRC4C-KIAA0645 fusion evidence in 1000G.

```
time (
    tabix onekg_fusion.bedpe.gz 11:40135753-41481323 \
        | grep LRRC4C | grep KIAA0645
)
# 11      40135752        41481323        neg     LRRC4C  22      32149943        32303012        pos     KIAA0645   79      73
#
# real    0m0.032s
# user    0m0.035s
#sys     0m0.025s
```
Here, 73 unique samples have at least 1 supporting read pair for LRRC4C-KIAA0645, and the total amount of supporting read pairs is 79.

## Details

Column descriptions
```
1: left gene chromosome
2: left gene start coordinate
3: left gene end coordinate
4: left gene strand
5: left gene name
6: right gene chromosome
7: right gene start coordinate
8: right gene end coordinate
9: right gene strand
10: right gene name
11: read count
12: sample count
```


