# Genefusion

## Install
```
https://github.com/jakekrol/genefusion.git && cd genefusion && pip install -e .
```

## GIGGLE interface

`gf-giggle_sharded` CLI script performs sharded giggle queries
```
# serial
gf-giggle_sharded /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/scratch/2024-12-15-giggle_sharded

# parallel
gf-giggle_sharded /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/scratch/2024-12-15-giggle_sharded True shard True 4
```

- Dependency: `genefusion_giggle.sh` script

- Output file name pattern: `shard_#.chrm.strand.gene.start.end.giggle`

## DNA gene2gene evidence

### Population

- Input: 1) giggle output file and 2) a BED file of GRCh37 genes
- Output: target regions

```
# pipe or provide giggle files as stdin to gargs
ls | gargs --log ../gargs.log -p 64 -o "$bedtools intersect -a $genefile -b <(cut -f 5-7 {0}) > $outdir/{0}.g2g"
```

### Sample-wise gene2gene evidence

Same as population, except split GIGGLE output files by the sample column, then do intersect.