#!/usr/bin/env bash

mapfile -t files < <(cat /data/jake/gene-fusion/data/meta/tumour_samples.txt)
(
    cd /data/jake/gene-fusion/data/cancer_data/prostate/prostate_sort
    for file in ${files[@]}; do
        cp $file ../prostate_tumour_sort/$file
    done
)
