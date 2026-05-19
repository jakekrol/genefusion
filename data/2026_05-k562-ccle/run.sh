#!/usr/bin/env bash

url1='https://zenodo.org/records/13363154/files/SRR521460_1.fastq.20M.fq.gz?download=1'
url2='https://zenodo.org/records/13363154/files/SRR521460_2.fastq.20M.fq.gz?download=1'

wget -O SRR521460_1.fastq.gz $url1
wget -O SRR521460_2.fastq.gz $url2