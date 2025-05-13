# GIGGLE2fusion pipeline

- gargs is just a parallelization tool, like xargs or gnu parallel
    - stdin redirection is usually a tab delimited file
    - {#} uses the #-th column value in the input file

1. GIGGLE search

- i: index
- g: gene name
- c: chromosome
- l: start base pair
- r: stop base pair
- o: outfile

```
gargs -p 60 --log=g.log -o "scripts/shell/genefusion_giggle.sh -i {0} -g {4} -c {1} -l {2} -r {3} -s {5} -o {6}" < i.txt
```

2. clean GIGGLE/excord file

- remove non-standard and "-1" chromosomes
- input {0}: GIGGLE/excord file
- output {1}: GIGGLE/excord file with some rows filtered

```
gargs -p 60 --log=g.log -o "./cln_excord.sh {0} {1}" < i.txt
```

3. temporarily swap the left and right hand interval

- fact: we only use the right-hand interval to find fusions; since, each GIGGLE file corresponds to a source gene
- purpose: 

swapping allows bedtools intersect of the right-hand interval with the left-hand preserved as metadata.
i want to keep both left/right intervals for an experiment comparing intervals of TP fusions and random PE evidence of genes

- input {0}: GIGGLE clean file
- output {1}: GIGGLE clean file with left and right hand intervals swapped

```
gargs -p 60 --log=g.log -o "./swap_intervals.sh {0} {1}" < i.txt
```

4. intersect

- purpose: 
    - fact: giggle search result returns all discordant paired end reads
    - goal: want to know if the other read (not the one in query gene) lands on another gene (a fusion!)

- input {0}: GIGGLE swapped file
- output {1}: BED file of intersections of right-hand interval and a gene bed file

```
gargs -p 60 --log=g.log -o "./intersect_swapped.sh {0} {1}" < i.txt
```

5. unswap intervals

- purpose: remove any confusion downstream by unswapping the left and right hand columns

- details:

- c1-5 is from genefile (the fusion target)
- c6-9 is right-hand (the fusion target)
- c10-13 is left-hand (the fusion source)
- c14 is ? (metadata)
- c15 is sample (metadata)

```
gargs -p 60 --log=g.log -o "./unswap_intervals.sh {0} {1}" < i.txt
```

6. clean sample names

- purpose: remove dir prefix and extension suffix from sample IDs

```
dinter_unswap='dir with unswapped intersect files'
tissue=ovary
cd $dinter_unswap
# test without -i flag first if not confident
# prefix
ls | gargs -p 60 --log=../g.log -o "sed -i 's|ovary_sort/||' {0}"
# suffix
ls | gargs -p 60 --log=../g.log -o "sed -i 's/\.excord\.bed\.gz//' {0}"
```

7. compute file id to specimen type mapping

- purpose: map PCAWG file ids to tumor/normal specimen
- input: PCAWG metadata file
- output: 2 column table mapping file ids to specimen type
- details:
    - only need to do once. creates mapping for all BAM files in all PCAWG projects
 
```
./build_fid2type.py -m <metadata_file> -o <output_mapping table>
```

8. add specimen type column

- input: 1) GIGGLE file (`-i`), 2) sample/file-id column index (`-s`), and 3) path to file id to specimen type lookup table (`-l`) (step 7)
- output: GIGGLE file with sample column appended (`-o`)

```
./add_sampletype.py -i <giggle_file> -s <sample_col_idx> -l <path_to_lookup_tbl> -o <outfile>
```

9. split files by specimen

- input: GIGGLE intersected file with specimen column added
- output: split files by the specimen column (expect 2 * number of input files, bc split by tumor/normal)
- details:
    - {0} input GIGGLE file
    - {1} out dir for split files (usually constant)
```
gargs -p 60 --log=g.log -o "./specimensplit_intersect.sh {0} {1}" < w.txt
```

10. remove the tumor/normal prefix

- input: filenames
- output: renamed files
- details
    - do this for `giggleinter_(tumor/normal)_final` and `pop_(tumor/normal)_fusion_counts`
```
# example
cd <dir>
ls > ../z.txt
cd ../
# if the prefix is 'tumour.'
sed "s|tumour\.||" z.txt > y.txt
paste z.txt y.txt > w.txt
cd <dir>
gargs -p 60 -o "mv {0} {1}" ../w.txt
```

11. count fusions

- purpose: for each gene count its fusion partners from paired-end (PE) read evidence
- input {0}: path to gene-wise intersected GIGGLE file (population or sample is fine)
- output {1}: path to count file of partner fusions
```
gargs -p 60 --log=g.log -o "./count_fusions.sh {0} {1}" < i.txt
```

12. build a fusion PE count table

- purpose: combine fusion counts (population-wise) into a single 3 column table
- input: path to directory of gene-wise fusion count files
- output: path to 3 column gene fusion table (left, right, PE-read count)
- details:
    - make sure you only use files for a distinct specimen type (tumor/normal), no mixing
    - long runtime
        - could speed up w/ map-reduce paradigm s.t. each file is processed in parallel, but this requires a bit of setup
```
./build_fusionpe_tbl2.py -i pop_normal_fusion_counts -o pop_normal_fusions.tsv -l build_fusion_pe_tbl2-normal.log
```

13. count number of samples with any evidence for a fusion

- purpose: for each fusion get the number of samples with any PE evidence for it
- input: path to directory of gene-wise, speciment-split, unswapped intersect files
- output: path to 3 column table (left, right, # of samples w/ >=1 PE-evidence)
- details:
    - make sure you only use files for a distinct specimen type (tumor/normal), no mixing
    - long runtime
        - could speed up w/ map-reduce paradigm s.t. each file is processed in parallel, but this requires a bit of setup
```
./count_samples_w1.py -i giggleinter_final_normal -o sample_counts_normal.tsv -l count_samples_w1-normal.log
```

14. compute gene total burden

- purpose: for each gene, sum its PE-evidence in all non-self genes (total burden)
- input: path to fusion PE count table (from 12)
- output: path to 2 column table (gene, burden)

```
./burden_total.py -i pop_tumor_fusions.tsv -o burden_total_tumor.tsv
```
