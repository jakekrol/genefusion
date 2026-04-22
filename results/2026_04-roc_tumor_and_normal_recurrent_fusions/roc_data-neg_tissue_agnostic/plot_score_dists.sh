#!/usr/bin/env bash

MEAN_SCRIPT=/data/jake/rl-tools/analysis/mean
# HIST_SCRIPT=/data/jake/rl-tools/plot/hist.py
DENSITY_SCRIPT=/data/jake/rl-tools/plot/density.py
XMIN=-1
XMAX=1

mapfile -t files < <(ls roc_input*.tsv)

for f in "${files[@]}"; do
    # pos
    tmp_pos=$(mktemp --suffix=.pos.tmp)
    trap "rm -f $tmp_pos" EXIT
    outfile="${f%.tsv}.density.png"
    tail -n +2 $f | awk '$4 == 1 {print $3}' > $tmp_pos
    mean="$(cat $tmp_pos | $MEAN_SCRIPT | cut -f 1 -d ',')"
    # cat $tmp_pos | \
    #     $DENSITY_SCRIPT -o $outfile --x_min $XMIN --x_max $XMAX \
    #         --axvline $mean

    # neg
    tmp_neg=$(mktemp --suffix=.neg.tmp)
    trap "rm -f $tmp_neg" EXIT
    tail -n +2 $f | awk '$4 == 0 {print $3}' > $tmp_neg
    mean="$(cat $tmp_neg | $MEAN_SCRIPT | cut -f 1 -d ',')"
    $DENSITY_SCRIPT -i "$tmp_pos,$tmp_neg" -o $outfile --show_median --names "Positive,Negative" \
        --ylabel "Score" --title $f
done
