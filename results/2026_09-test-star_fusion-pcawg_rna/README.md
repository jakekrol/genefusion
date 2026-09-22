example code to download, convert bam2fastq, and run star-fusion on pcawg bams from EGA

EGAC00001000010: ICGC PCAWG Dataset for RNA-Seq BAM aligned using Star. Project: CLLE-ES.

example bam: EGAF00001719756

```
# login
egafetch-linux-amd64 auth login
# download
egafetch-linux-amd64 download EGAF00001719756
# sort by query name
samtools sort -@ 8 -n -o PCAWG.d0e5bdc0-f887-11e4-a9de-4b74fb958afb.STAR.v1.namesort.bam PCAWG.d0e5bdc0-f887-11e4-a9de-4b74fb958afb.STAR.v1.bam
# bamtofastq
bedtools bamtofastq \
    -i PCAWG.d0e5bdc0-f887-11e4-a9de-4b74fb958afb.STAR.v1.namesort.bam \
    -fq r1.fastq \
    -fq2 r2.fastq
conda activate star_fusion
STAR-Fusion --CPU 8 --genome_lib_dir $GENOME_LIB_DIR \
                        --left_fq r1.fastq \
                        --right_fq r2.fastq \
                        --output_dir star_fusion_outdir
```
