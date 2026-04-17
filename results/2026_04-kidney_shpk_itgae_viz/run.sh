#!/usr/bin/env bash
# bam=../2026_02-samplot_validation/plot_data/bams/liver_POLR3GL--LIX1L_DEL_PCAWG.4dd18422-e203-11e4-93ca-dd97a5ed4642.STAR.v1.bam
bam=../2026_02-samplot_validation/plot_data/bams/kidney_SHPK--ITGAE_DEL_PCAWG.8631d5c8-dd36-11e4-b9d1-4999c254ba06.STAR.v1.bam
cp $bam ${bam}.bai .
bam=$(basename $bam)
# transcript_file=..//2026_02-samplot_validation/plot_data/annotations/POLR3GL--LIX1L.bed.gz
 transcript_file=../2026_02-samplot_validation/plot_data/annotations/SHPK--ITGAE.bed.gz
# cp $transcript_file ${transcript_file}.tbi .
transcript_file=$(basename $transcript_file)

# chr=1
# start=145456236
# end=145501669
chr=17
start=3468738
end=3705537
length=$((end - start))
window=$((length / 100))
name=SHPK--ITGAE

samplot plot \
	-b $bam -c $chr -s $start -e $end -o plot.png -n $name \
	-A $transcript_file \
	--max_coverage 100  \
	--window $window \
	-W 7 -H 4