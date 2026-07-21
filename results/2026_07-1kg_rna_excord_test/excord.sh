#!/usr/bin/env bash

input=SRR19762165.star_fusion.sort.bam
output=${input%.bam}.excord.bed.gz
cat $input | excord --discordantdistance 500 /dev/stdin | bgzip -c > $output