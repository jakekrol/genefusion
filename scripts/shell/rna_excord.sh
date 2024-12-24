#!/usr/bin/env bash
bedtools='/data/jake/bedtools.static.binary'
set -u
# try using excord on RNA seq

# getting smallest bam file
# aws s3 ls s3://layerlabcu/rna_cancerdata/ovary/ | grep 'bam$'  | sort -k3n | head
# 3.8 GB

# copy down a single bam and its index
# s3='s3://layerlabcu/rna_cancerdata/ovary/'
# bam='PCAWG.b112f618-c81a-4639-b24a-6521ee527d31.STAR.v1.bam'
# idx='PCAWG.b112f618-c81a-4639-b24a-6521ee527d31.STAR.v1.bam.bai'
# aws s3 cp "${s3}${bam}" .
# aws s3 cp "${s3}${idx}" .

# excord
# cat *.bam | excord --discordantdistance 500 /dev/stdin > test.out
# cat *.bam | excord --discordantdistance 600 /dev/stdin | bgzip -c > test.bed

dl_bam () {
    bam=$1
    out=$2
    s3_path=${3:-'s3://layerlabcu/rna_cancerdata/ovary/'}
    aws s3 cp "${s3_path}${bam}" $out
    aws s3 cp "${s3_path}${bam}.bai" $out.bai
}
export -f dl_bam
# bam='PCAWG.e32387b8-8245-4e8d-ad8f-23b1bc7060c9.TopHat2.v1.bam'
# out='test.bam'

### downloading 10 randomly sampled bams
# input='ten_bams.txt'

# gargs --log gargs.log -p 10 -o "dl_bam {0} {0}" < $input

do_excord () {
    bam=$1
    out="${bam%.*}.bed"
    cat $bam | excord --discordantdistance 500 /dev/stdin | bgzip -c > $out
}
export -f do_excord

### do excord
# for bam in $(ls *.bam); do
#     echo $bam
#     do_excord $bam
# done

### clean excord files
cln_excord ()
{
    in=$1
    out=$2
    zcat $in | awk '($1 !~ /^GL/ && $1 != "MT" && $1 != "*" && $2 != "-1" && $3 != "-1")' > "$out"
    # Explanation:
    # $1 != "*"  -> Keep lines where the first column does not begin with "GL"
    # $1 != "*"  -> Keep lines where the first column is not "*"
    # $2 != "-1" -> Keep lines where the second column is not "-1"
    # $3 != "-1" -> Keep lines where the third column is not "-1"
}
export -f cln_excord

# ls *.bed.gz | gargs --log gargs.log -p 10 -o "cln_excord {0} {0}.cln"

sort_and_zip ()
{
    in=$1
    out=$2
    $bedtools sort -i $in | bgzip -c > $out
}
export -f sort_and_zip

# ls *.bed.gz.cln | gargs --log sort_and_zip.log -p 10 -o "sort_and_zip {0} {0}.sorted.gz"

cln_chrm ()
{
    in=$1
    out=$2
    zcat $in | awk '($1 != "MT" && $1 !~ /^GL/ && $1 !~ /^hs/)' | bgzip -c > $out
}
export -f cln_chrm

ls *.sorted.gz | gargs --log cln_chrm.log -p 10 -o "cln_chrm {0} {0}.cln.gz"

# rename
# awk -F'.' '{print $1 "." $2 "." $3 "." $4 "." $8 "." $5 ".gz"}'
