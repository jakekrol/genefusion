#!/usr/bin/env bash

echo "# begin test_wrapper_giggle2fusion.sh"
t=$(mktemp -d "./test_wrapper_giggle2fusion.XXXXXX")
echo "# created temp outdir $t"
script=../scripts/python/wrapper_giggle2fusion.py

bed='../results/2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed'
base_dir=$t
index='/data/jake/stix-pcawg-dna/prostate_tumor.index.0'
steps="0"
dir_executables=$(realpath ../executables)
cpus=4
type="tumor"
template='../results/2025_10-fusion_search_template/giggle_search.template'
modality="dna"

echo "# testing parameters:"
echo "bed: $bed"
echo "base_dir: $base_dir"
echo "giggle_index: $index"
echo "steps: $steps"
echo "dir_executables: $dir_executables"
echo "processes: $cpus"
echo "type: $type"
echo "template: $template"
echo "modality: $modality"

cmd="$script \
    --bed $bed \
    --base_dir $base_dir \
    --giggle_index $index \
    --steps $steps \
    --dir_executables $dir_executables \
    --processes $cpus \
    --type $type \
    --template $template \
    --modality $modality"

echo "# running command:"
echo "$cmd"

$cmd

echo "# end test_wrapper_giggle2fusion.sh"
