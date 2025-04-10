#!/usr/bin/env bash

# gene-wise giggle search intersect excord hits with all other genes

# Usage: ./genefusion_giggle.sh -i index -f Homo_sapiens.GRCh37.82.chr.genes.bed -g ERG -c 21 -s neg -l 39751949 -r 40033704
echo "$0 $@"

while getopts "i:f:g:c:s:l:r:o:b" opt; do
  case $opt in
    i)
      index=$OPTARG
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

echo "searching  $gene ($chrm:$left-$right) in index $index"
(
  cd $index/..
  giggle search -v -i $index -r $chrm:$left-$right > $outfile # | cut -f 5-7 > $outfile
)

