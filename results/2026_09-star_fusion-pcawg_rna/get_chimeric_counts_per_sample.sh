#!/usr/bin/env bash

DIR_STAR_FUSION=completed
FILE_CHIMERIC=Chimeric.out.junction
OUTFILE="chimeric_counts.tsv"



printf "sample\tchimeric_count\n" > $OUTFILE
for f in "${DIR_STAR_FUSION}"/*/star-fusion-out/"${FILE_CHIMERIC}"; do
    echo $f 
    sample=$(echo $f | cut -f 2 -d'/')
    chimeric_count=$(tail -n +2 $f | wc -l)
    printf "${sample}\t${chimeric_count}\n" >> ${OUTFILE}
done