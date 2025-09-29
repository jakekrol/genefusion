#!/usr/bin/env bash

url="https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
wget $url
unzip rna_tissue_consensus.tsv.zip
echo "Downloaded and unzipped rna_tissue_consensus.tsv.zip"
