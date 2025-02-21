#!/usr/bin/env bash

# gene-wise giggle search intersect excord hits with all other genes

# Usage: ./genefusion_giggle.sh -i index -f Homo_sapiens.GRCh37.82.chr.genes.bed -g ERG -c 21 -s neg -l 39751949 -r 40033704
echo "$0 $@"

while getopts "i:f:g:c:s:l:r:o:b" opt; do
  case $opt in
    i)
      index=$OPTARG
      ;;
    f)
      genefile=$OPTARG
      ;;
    g)
      gene=$OPTARG
      ;;
    c)
      chrm=$OPTARG
      ;;
    s)
      strand=$OPTARG
      ;;
    l)
      left=$OPTARG
      ;;
    r)
      right=$OPTARG
      ;;
    o)
      outfile=$OPTARG
      ;;
    # b)
    #   outintersect=$OPTARG
    #   ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
date
if [ -z "${outfile}" ]; then
  outfile="$chrm.$strand.$gene.$left.$right.giggle"
fi
# if [ -z "${outintersect}" ]; then
#   outintersect="$chrm.$strand.$gene.$left.$right.intersect"
# fi

giggle search -v -i $index -r $chrm:$left-$right > $outfile # | cut -f 5-7 > $outfile
# rm excord hits with all -1 interval
# adjust to get bedtools intersect to work
#sed -i '/^-1/d' $outfile
# intersect with gene file
#bedtools intersect -a $genefile -b $outfile > $outintersect
date

#date
#region="21:39751949-40033704"
#index="/data/jake/genefusion/data/prostate/shards/shard_0/index"
#gene_file="/data/jake/genefusion/data/prostate/shards/shard_0/chr21.neg.test"
#hits=/data/jake/genefusion/data/prostate/shards/shard_0/ERG.hits
#giggle search -v -i $index -r $region | cut -f 5-7 > $hits
#sed -i '/^-1/d' $hits
#bedtools intersect -a $gene_file -b $hits > bedtools_intersect.out
## cut -f4 bedtools_intersect.out | sort | uniq -c | sort -n
#date
