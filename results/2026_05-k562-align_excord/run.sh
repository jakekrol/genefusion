#!/usr/bin/env bash

# rna processing
./star_index.sh
./star_align.sh

# index bams for excord
./index_bams.sh

# excord
./excord.sh

