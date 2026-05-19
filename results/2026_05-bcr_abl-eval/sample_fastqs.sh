#!/usr/bin/env bash

k562_fastq1="../../data/2026_05-k562-ccle/SRR521460_1.fastq.gz"
k562_fastq2="../../data/2026_05-k562-ccle/SRR521460_2.fastq.gz"
cp "${k562_fastq1}" ./k562-1.0-1.fastq.gz
cp "${k562_fastq2}" ./k562-1.0-2.fastq.gz
fractions=(0.2 0.4 0.6 0.8)

for fraction in "${fractions[@]}"; do
    echo "Sampling ${fraction} of the reads..."
    seqtk sample -s100 "${k562_fastq1}" "${fraction}" > "k562-${fraction}-1.fastq"
    seqtk sample -s100 "${k562_fastq2}" "${fraction}" > "k562-${fraction}-2.fastq"
done