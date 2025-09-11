#!/usr/bin/env bash

# make giggle2fusion template
./bed2input.py --bed ../2025_09-onekg_bed/grch37.bed \
    --output giggle2fusion.template
# setup pipeline
mkdir g2f
cd g2f
ln -s ../../../scripts/python/wrapper_giggle2fusion.py
cp ../giggle2fusion.template .
# run pipeline
./wrapper_giggle2fusion.py \
    --template giggle2fusion.template \
    -d $(pwd) \
    -g /data/stix/1kg/low_coverage/alt_sort_b \
    -s 0,1,2,3,4,5,6 \
    -p 60 \
    -t normal
