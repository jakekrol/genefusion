#!/usr/bin/env bash


outdir=$(pwd)
r_erg=21:39751949-40033704
r_tmprss2=21:42836478-42903043
index=/data/jake/stix-pcawg-dna
# giggle search
outsearch=$outdir/giggle.erg.out
( 
    cd $index
    mapfile -t shards < <(ls $index | grep prostate | grep tumor | grep index)
    for shard in "${shards[@]}"; do
        echo "Processing $shard"
        giggle search -i $shard -r $r_erg -v >> $outsearch
    done
)
outsearchfmt=$outdir/giggle.erg.fmt
gf-fmt_giggle_fusion -g ERG -r $r_erg -i $outsearch -o $outsearchfmt -m "index=$index"

# clean
script=../../scripts/shell/cln_excord.sh
outclean=$outdir/giggle.erg.fmt.cln
$script -i $outsearchfmt -o $outclean

# swap intervals
script=../../scripts/shell/swap_intervals.sh
outswap=$outdir/giggle.erg.fmt.cln.swap
$script -i $outclean -o $outswap

# intersect
script=../../scripts/shell/intersect_swapped.sh
genebed=../2025_04-gene_bedfile_cln/grch37.genes.bed
outintersect=$outdir/giggle.erg.fmt.cln.swap.intersect
$script -i $outswap -g $genebed -o $outintersect

# unswap intervals
script=../../scripts/shell/unswap_intervals.sh
outunswap=$outdir/giggle.erg.fmt.cln.swap.intersect.unswap
$script -i $outintersect -o $outunswap

# filter for tmprss2 only
outtmprss2=$outdir/giggle.erg.fmt.cln.swap.intersect.unswap.tmprss2
awk -F'\t' -v OFS='\t' '
    /^#/ { print; next }
    !header_seen { print; header_seen=1; next }
    $4 == "TMPRSS2" { print }
' "$outunswap" > "$outtmprss2"



