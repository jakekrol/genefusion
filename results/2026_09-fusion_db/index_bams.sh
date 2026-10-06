#!/usr/bin/env bash
DIRBAM='bams'
THREADS=4
mapfile -t bams < <(ls $DIRBAM/*/*.bam)
echo "${bams[@]}"

for f in "${bams[@]}"; do
    rm "${f}.bai" || echo "# no index to remove for $f"
    samtools index -@ ${THREADS} $f || echo "# failed to index $f"
done