#!/usr/bin/env bash
# conda activate samplot_env
# args
bucket=${1}
# const
CPUS=30
BAMS='bams'
OUTDIR='samplots'
OUTDIR=$(realpath $OUTDIR)
FONTSIZE=20
CHR=17
START=3468738
END=3705537
LENGTH=$((END - START))
WINDOW=$((LENGTH / 100))
PAD=100000
PAD_START=$((START - PAD))
PAD_END=$((END + PAD))
REGION="${CHR}:${PAD_START}-${PAD_END}"
ANNOTATION_FILE="SHPK--ITGAE.bed.gz"
ANNOTATION_FILE=$(realpath $ANNOTATION_FILE)
MAX_COVERAGE=100
W=12
H=6
ANNOTATION_FILENAMES="."
echo "bucket: $bucket"
echo "bams: $BAMS"
echo "outdir: $OUTDIR"
echo "chr: $CHR"
echo "start: $START"
echo "end: $END"
echo "length: $LENGTH"
echo "window: $WINDOW"
echo "W: $W"
echo "H: $H"
echo "PAD: $PAD"
mkdir -p $OUTDIR
mkdir -p $BAMS
# list kidney rna-seq files
# aws s3 ls s3://$bucket/rna_cancerdata/kidney/ | grep PCAWG | awk '{print $4}' > kidney_rna_seq_files.txt
# get indices
# mapfile -t bais < <(grep bai$ kidney_rna_seq_files.txt )
# echo "#bais: ${bais[@]}"
# for bai in "${bais[@]}"; do
#   echo "bai: $bai"
#   aws s3 cp s3://$bucket/rna_cancerdata/kidney/$bai $BAMS/
# done

# # get bam regions
# grep bam$ kidney_rna_seq_files.txt > bams.txt
# # mapfile -t bams < <(grep bam$ kidney_rna_seq_files.txt )
# cd $BAMS
# pwd
# cat ../bams.txt | gargs -p $CPUS --log=../gargs.log "samtools view -b s3://$bucket/rna_cancerdata/kidney/{0} $REGION > {0}"




cd $BAMS
pwd

for bam in $(ls *.bam); do
    echo "bam: $bam"
    name=$bam
    samplot plot \
        -b $bam -c $CHR -s $START -e $END -o $OUTDIR/${name}.png -n $name \
        -A $ANNOTATION_FILE \
        --max_coverage $MAX_COVERAGE  \
        --window $WINDOW \
        -W $W -H $H \
        --annotation_fontsize $FONTSIZE \
        --annotation_filenames $ANNOTATION_FILENAMES \
        --xaxis_label_fontsize $((FONTSIZE-10)) \
        --yaxis_label_fontsize $((FONTSIZE-10))
done

