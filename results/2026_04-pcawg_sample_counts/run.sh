#!/usr/bin/env bash
set -euo pipefail
./get_sample_counts.py \
    -d /data/jake/stix-pcawg-dna \
    -r /data/jake/stix-pcawg-rna \
    -o .
./sample_counts_plot.py

