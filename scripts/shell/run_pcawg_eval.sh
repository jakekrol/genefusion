#!/usr/bin/env bash

# Comment out strict error handling temporarily for debugging
# set -euo pipefail

# Default values
n=10000
seed=0

# Usage function
usage() {
    echo "Usage: $0 -d DATAFILE -t TISSUE -o OUTDIR [-n SAMPLE_SIZE] [-s SEED]"
    echo "  -d, --datafile    Input data file (required)"
    echo "  -t, --tissue      Tissue type (required)" 
    echo "  -o, --outdir      Output directory (required)"
    echo "  -n, --sample-size Sample size for random sampling (default: 10000)"
    echo "  -s, --seed        Random seed (default: 0)"
    echo "  -h, --help        Show this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--datafile)
            datafile="$2"
            shift 2
            ;;
        -t|--tissue)
            tissue="$2"
            shift 2
            ;;
        -o|--outdir)
            outdir="$2"
            shift 2
            ;;
        -n|--sample-size)
            n="$2"
            shift 2
            ;;
        -s|--seed)
            seed="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "${datafile:-}" ]]; then
    echo "Error: Datafile is required (-d/--datafile)"
    exit 1
fi

if [[ -z "${tissue:-}" ]]; then
    echo "Error: Tissue is required (-t/--tissue)"
    exit 1
fi

if [[ -z "${outdir:-}" ]]; then
    echo "Error: Output directory is required (-o/--outdir)"
    exit 1
fi

# Validate numeric arguments
if ! [[ "$n" =~ ^[0-9]+$ ]]; then
    echo "Error: Sample size must be a positive integer"
    exit 1
fi

if ! [[ "$seed" =~ ^[0-9]+$ ]]; then
    echo "Error: Seed must be a non-negative integer"
    exit 1
fi

# Create output directory if it doesn't exist
if [[ ! -d "$outdir" ]]; then
    echo "Creating output directory: $outdir"
    mkdir -p "$outdir"
fi

# Change to output directory first
cd "$outdir" || {
    echo "Error: Could not change to output directory $outdir"
    exit 1
}

# files
pcawgfile=/data/jake/genefusion/results/2025_06-fusion_gold_standard/fusion_gold_standard.tsv

# scripts
script_join="join_ddb.py"
script_scatter_read_sample="/data/jake/genefusion/scripts/python/scatter_pcawg_eval.py"
script_rand="/data/jake/rl-tools/sim/rand.sh"
script_hist_pcawg="/data/jake/genefusion/scripts/python/hist_pcawg_fusions.py"

# File existence checks (datafile assumed to be in outdir)
if [[ ! -f "$datafile" ]]; then
    echo "Error: Datafile '$datafile' does not exist in $outdir"
    exit 1
fi

if [[ ! -f "$pcawgfile" ]]; then
    echo "Error: PCAWG file '$pcawgfile' does not exist"
    exit 1
fi

# Script existence checks
for script in "$script_scatter_read_sample" "$script_rand" "$script_hist_pcawg"; do
    if [[ ! -x "$script" ]]; then
        echo "Error: Script '$script' not found or not executable"
        exit 1
    fi
done

echo "Starting PCAWG evaluation for $tissue tissue..."
echo "Parameters: datafile=$datafile, tissue=$tissue, outdir=$outdir, n=$n, seed=$seed"


echo "Filtering PCAWG fusions for $tissue tissue..."
# filter pcawg fusions for tissue
if ! grep "pcawg" "$pcawgfile" | grep "$tissue" > "pcawg-$tissue-fusions.tsv"; then
    echo "Error: No PCAWG fusions found for tissue '$tissue'"
    exit 1
fi

pcawg_count=$(wc -l < "pcawg-$tissue-fusions.tsv")
echo "Found $pcawg_count PCAWG fusions for $tissue"

echo "Joining data with PCAWG fusions..."
# subset data for pcawg fusions
# we have to join with unfiltered file to avoid missing fusion data
join_ddb.py \
    --left <(printf "left\tright\n" && cut -f1,2 "pcawg-$tissue-fusions.tsv") \
    --right scored_dna1.0_t0.5_r0.5_u50_expr_anno_sort_no_dups_anno_filt.tsv \
    --type left \
    --keys left,right \
    --output "pcawg-$tissue-fusions-scored.tsv"

/data/jake/rl-tools/wrangle/fillna.sh \
    -i "pcawg-$tissue-fusions-scored.tsv" \
    -o "pcawg-$tissue-fusions-scored.tsv.tmp" \
    -f 0
mv "pcawg-$tissue-fusions-scored.tsv.tmp" "pcawg-$tissue-fusions-scored.tsv"

scored_count=$(tail -n +2 "pcawg-$tissue-fusions-scored.tsv" | wc -l)
echo "Scored $scored_count PCAWG fusions"

echo "Creating scatter plot data..."
# read v sample tumor
tail -n +2 "$datafile" | head -n "$n" | cut -f 3,4 > "top${n}_reads_samples-scatter.txt"
sed -i 's/$/\ttop'"$n"'-score/' "top${n}_reads_samples-scatter.txt"
tail -n +2 "pcawg-$tissue-fusions-scored.tsv" | cut -f 3,4 > pcawg-scatter.txt
sed -i 's/$/\tpcawg/' pcawg-scatter.txt
cat "top${n}_reads_samples-scatter.txt" pcawg-scatter.txt > tumor-read-v-sample-scatter.input

echo "Plotting scatter plot..."
echo "Running: $script_scatter_read_sample --input tumor-read-v-sample-scatter.input --output $tissue-tumor-read-v-sample-scatter.png"
if ! "$script_scatter_read_sample" \
    --input tumor-read-v-sample-scatter.input \
    --output "$tissue-tumor-read-v-sample-scatter.png" \
    --xlabel "Num. supporting reads" \
    --ylabel "Num. supporting samples" \
    --title "$tissue fusion tumor evidence"; then
    echo "Error: Scatter plot failed"
    exit 1
fi

echo "Sampling random scores..."
# sample random scores
if ! tail -n +2 "$datafile" | cut -f 9 | "$script_rand" "$n" "$seed" > "${n}-rand-scores.txt"; then
    echo "Error: Random sampling failed"
    exit 1
fi

echo "Creating score histogram..."
if ! "$script_hist_pcawg" \
    --pcawg <(tail -n +2 "pcawg-$tissue-fusions-scored.tsv" | cut -f 9) \
    --null "${n}-rand-scores.txt" \
    --output rand_v_pcawg_hist.png \
    --title "$tissue scores" \
    --xlabel "Fusion Score" \
    --ylabel "Density" \
    --density; then
    echo "Error: Histogram creation failed"
    exit 1
fi

echo "PCAWG evaluation completed successfully!"
echo "Output files created in: $outdir"
echo "  - $tissue-tumor-read-v-sample-scatter.png"
echo "  - rand_v_pcawg_hist.png"


