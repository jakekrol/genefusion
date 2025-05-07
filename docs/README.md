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
ls | gargs -p 60 --log=../g.log -o "sed -i '|ovary_sort/||' {0}"
# suffix
ls | gargs -p 60 --log=../g.log -o "sed -i 's/\.excord\.bed\.gz//' {0}"
```

7. sample splitting

- purpose: split unswapped, gene-intersected, GIGGLE BED file by the sample column. Each sample gets its own file
- input {0}: in file
- output {1}: directory to house sample-wise, gene-wise intersected GIGGLE files
- details:
    - usually {1} is constant

```
gargs -p 60 --log=g.log -o "./samplesplit_intersect.sh {0} {1}" < i.txt
```

8. migrate samples

- purpose: move sample-wise, gene-wise, gene-intersected GIGGLE BED files by looking up whether the sample is tumor/normal
- parameters:
    - `-i` input dir of sample GIGGLE files
    - `-n` output normal dir
    - `-t` output tumor dir
    - `-s` path to `fileid2sample_type.py` script, for lookup
    - `-l` logfile path
 
- details:
    - slow startup as the id2sampletype cache is built, then fast
    - serial
        - consider using a pre-computed id2sampletype cache to parallelize this

 
```
./migrate_sample2.py -i $din -n $dnormal -t $dtumor -s /data/jake/genefusion/scripts/python/fileid2sample_type.py -l /data/jake/genefusion/logs/migrate_sample2.log
```

9. aggregate by sample type

- purpose: combine sample-wise, gene-wise intersected GIGGLE files of a sample type (tumor/normal) s.t. each gene gets one file again
    - useful for computing pop. scale stats
- parameters:
    - `-i` input dir of norm/tumor sample-wise, gene-wise intersected GIGGLE files of a specific sample type (from step 8 output)
    - `-o` output dir of the aggregated files
    - `-l` path to logfile
- details:
    - serial, not parallelized
    - uses file name to parse sample and gene
```
./agg_sample_fusions.py -i $din -o $dout -l /data/jake/genefusion/logs/agg_sample_fusions.log
```

10. count fusions

- purpose: for each gene count its fusion partners from paired-end (PE) read evidence
- input {0}: path to gene-wise intersected GIGGLE file (population or sample is fine)
- output {1}: path to count file of partner fusions
```
gargs -p 60 --log=g.log -o "./count_fusions.sh {0} {1}" < i.txt
```
