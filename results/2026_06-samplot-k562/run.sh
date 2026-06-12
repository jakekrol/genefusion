#!/usr/bin/env bash

bam=../2026_05-k562-zenodo/star_fusion_outdir2/Aligned.out.bam
annotations=../2026_06-grch37-transcripts/gencode.v19.annotation.gtf.exons.bed
fusion=../2026_05-k562-zenodo/star_fusion_outdir2/star-fusion.fusion_predictions.tsv

chrom=10
start=51606988
end=51732772

cp $bam $annotations $fusion .
bam=$(basename $bam)
annotations=$(basename $annotations)
fusion=$(basename $fusion)

# sort bam
samtools sort -o sorted_$bam -m 4G --threads 16 $bam
bam=sorted_$bam
# index bam
# samtools index $bam

# prepare transcript bed - sort before indexing
sort -k1,1 -k2,2n $annotations > ${annotations}.sorted
bgzip -c ${annotations}.sorted > $annotations.gz
tabix -p bed $annotations.gz
annotations=$annotations.gz

# subset bed
zcat $annotations | grep -E "TIMM|PARG" | bgzip -c > TIMM_PARG.bed.gz
tabix -p bed TIMM_PARG.bed.gz

# plot
samplot plot \
	-b $bam \
	-n k562 \
	-c $chrom \
	-s $start \
	-e $end \
	-t DEL \
	-o p.png \
	-A TIMM_PARG.bed.gz \
	--max_coverage 100
