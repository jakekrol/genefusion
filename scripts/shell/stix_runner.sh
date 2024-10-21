#!/usr/bin/env bash
# for query gene, search for gene fusions; same chr and strand

# ex: stix_runner.sh > stix_runner.log 2>&1

# args
set -u
gene_i=${1:-ERG}
chr=${2:-21}
l_s=${3:-39751949}
l_e=${4:-40033704}
gene_file=${5:-chr21.neg}
outdir=${6:-/data/jake/gene-fusion/data/fusions}
[ -d "$outdir" ] || { echo "Directory '$outdir' does not exist. Exiting."; exit 1; }
echo "BEGIN"
date
echo -e "#gene: $gene_i\tchr: $chr\tl_s: $l_s\tl_e: $l_e\tgene_file: $gene_file\toutdir: $outdir" 
# query stix with gene_i=left and gene_j:= rows in gene bed file
while read -r _ r_s r_e gene_j _; do
  echo "#gene_j: $gene_j"
  if [ "$r_s" -gt "$l_e" ]; then
    stix -i alt_sort_b -d 1kg.ped.db -s 500 -t DEL -l $chr:$l_s-$l_e -r $chr:$r_s-$r_e > $outdir/$gene_i.$gene_j.stix
  else
    echo "#skipping $gene_j (preceeds query)"
  fi
done < $gene_file
date
echo "END"

# next:
# perform query for each unique gene,chromosome,strand 
