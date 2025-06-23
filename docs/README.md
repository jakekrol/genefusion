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

- input: giggle intersect unswapped file
- output: sample column with index prefix and .excord.bed.gz extension removed. leaving only the PCAWG fileid
- purpose: remove dir prefix and extension suffix from sample IDs

```
./cln_sample_name.py -i {0} -o {1}
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
./add_specimentype.py -i <giggle_file> -s <sample_col_idx> -l <path_to_lookup_tbl> -o <outfile>
```

9. split files by specimen

- input: GIGGLE intersected file with specimen column added
- output: split files by the specimen column (expect 2 * number of input files, bc split by tumor/normal)
- details:
    - {0} input GIGGLE file
    - {1} out dir for split files (usually constant)
    - specimencolidx is usually 16 (1-indexed)
```
gargs -p 60 --log=g.log -o "./specimensplit_intersect.sh -i {0} -o {1} -s <specimencolidx>" < w.txt
```

10. move by speciment and remove the tumor/normal prefix

- input: dir of files with specimen prefixes
- output: moved files to dirs by specimen type
```
./migrate_specimen.py -i <input_dir> -t <tumor_dest> -n <normal_dest> -p <n_procs>
```

11. count fusions

- purpose: for each gene count its fusion partners from paired-end (PE) read evidence
- input {0}: path to gene-wise intersected GIGGLE file (population or sample is fine)
- output {1}: path to count file of partner fusions
- details:
    - `-z`: hack to extract left gene from file name
    - `-r`: right gene column index (1-indexed,usually 4)
```
gargs -p 60 --log=g.log -o "./count_fusions.py -i {0} -o {1} -z -r 4" < i.txt
```

12. build a fusion PE count table

- purpose: combine fusion counts from count fusions step (population-wise) into a single table
- input: path to directory of gene-wise fusion count files
- output: path to aggregate table

```
./agg_pe_counts.py -i pop_normal_fusion_counts -o pop_normal_fusions.tsv
```

13. count number of samples with any evidence for a fusion

- purpose: for each fusion get the number of samples with any PE evidence for it
- input: giggle intersect unswapped file (ideally with cleaned sample name)
- output: path to 3 column table (left, right, # of samples w/ >=1 PE-evidence)
- details:
    - make sure you only use files for a distinct specimen type (tumor/normal), no mixing
    - `-z` will extract left gene from input file name
    - operates file-wise, requires an aggregation step
```
./distinct_sample_counts.py -i <left_giggle> -o <out> -r 4 -s 15 -z 
```

can aggregate with tail -n +2 -q pop_tumor_fusion_sample_counts/* >> pop_tumor_fusion_sample_counts.tsv

14. compute gene total burden

- purpose: for each gene, sum its PE-evidence in all non-self genes (total burden)
- input: path to fusion PE count table (from 12)
- output: path to 2 column table (gene, burden)

```
./burden_total.py -i pop_tumor_fusions.tsv -o burden_total_tumor.tsv --header
```

15. join sample count onto pe count

```
./join.py -x pe_counts_tumor.tsv -y sample_counts_tumor_new.tsv -t left -o tumor_pe_and_sample_new.tsv -k left,right
```

16. add total burden for left and right gene columns

```
# left
./add_burden_col.py -f tumor_pe_and_sample_new.tsv -b burden_total_tumor_new.tsv -o tumor_pe_sample_and_burden_new.tsv -k1 0 -k2 0 -n burden_total_left -hf -hb
# right
./add_burden_col.py -f tumor_pe_sample_and_burden_new.tsv -b burden_total_tumor_new.tsv tumor_pe_sample_and_burden_new.tsv -k1 1 -k2 0 -n burden_total_right -hf -hb
```

17. add sample density column

```
./add_sample_density.py -i tumor_pe_sample_burden_both_new.tsv -o tumor_pe_sample_burden_both_density.tsv
```

18. add burden product column

```
./add_burden_product.py -i tumor_pe_sample_burden_both_density.tsv -o tumor_pe_sample_burden_both_density_product.tsv
```

19. breakpoint rand test (Clark Evans R)

```
gargs -p 60 --log=g.log -o "./clark_evans_R -i {0} -o {1} -z -n 10" < clark_evans_R.input
```

