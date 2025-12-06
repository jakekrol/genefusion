#!/usr/bin/env python3

# visualize scoring experiment
# heatmap of median calibration rank under diff params

import argparse
import matplotlib.pyplot as plt
import os,sys
import shutil
import subprocess
import pandas as pd
import glob
import re
import numpy as np
from matplotlib.patches import FancyArrowPatch
import matplotlib.colors as mcolors
import duckdb as ddb

PARAMS = {
    'w_dna': [0.1, 0.5, 0.9],
    'w_tumor': [0.1, 0.5, 0.9],
    'w_read': [0.1, 0.5, 0.9],
    'upper': [10,50,100]
}

parser = argparse.ArgumentParser(description='Visualize fusion scoring')
parser.add_argument('-i', '--input', type=str, required=True, help='Input directory of scoring experiment')
parser.add_argument('--cache', type=str, default=None, help='Use cached summary file if exists')
parser.add_argument('--fontsize', type=int, default=12, help='Font size for plots')
args = parser.parse_args()

def rank_finder(score, scorefile, col):
    # maybe use score instead of x,y to speed up
    cmd = f"tail -n +2 {scorefile} | " \
        f"cut -f {col} | " \
        f"awk -v t={score} '$1 >= t {{rank=NR}} END {{print rank}}'"
    print(f"Running command: {cmd}")
    # launch cmd using subprocess
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    rank_str = result.stdout.strip()
    if rank_str:
        return int(rank_str) - 1  # minus 1 for header
    else:
        return -1

def fname2params(scorefile):
    # extract params from filename
    match = re.search(r'dna([^_]*)_t([^_]*)_r([^_]*)_u([^_]*)', scorefile)
    if match:
        w_dna = float(match.group(1))
        w_tumor = float(match.group(2))
        w_read = float(match.group(3))
        upper = int(match.group(4))
    return w_dna, w_tumor, w_read, upper

        
def cache_check(args):
    # check for cached results
    if args.cache:
        if os.path.exists(args.cache):
            df_results = pd.read_csv(args.cache, sep='\t')
            print(f"Loaded cached results from {args.cache}")
            return df_results
        else:
            raise ValueError(f"Cache file {args.cache} not found")
    else:
        return None


def data_gather(args):
    # gather all score files
    scorefiles = []
    for file in os.listdir(args.input):
        if file.endswith('_sort.tsv'):
            scorefile = os.path.join(args.input, file)
            scorefiles.append(scorefile)
    print(f"Found {len(scorefiles)} score files")

    # for each score file
    # for each calibration fusion
    # get rank
    data = {
        'w_dna': [],
        'w_tumor': [],
        'w_read': [],
        'upper': [],
        'median_rank': []
    }
    for scorefile in scorefiles:
        print(f"Processing score file: {scorefile}")
        w_dna, w_tumor, w_read, upper = fname2params(scorefile)
        print("extracted params:", "w_dna:", w_dna, "w_tumor:", w_tumor, "w_read:", w_read, "upper:", upper)
        outfile = os.path.basename(scorefile).replace('_sort.tsv', '_cal_rank.tsv')
        # sanity check
        assert outfile != scorefile, "Output file should be different from scorefile"
        print(f"Loading {scorefile} into duckdb")
        conn = ddb.connect(database=':memory:')
        df_score = conn.execute(f'SELECT fusion_score FROM read_csv_auto("{scorefile}", delim="\t")').df()
        cal_file =  f"scored_dna{w_dna}_t{w_tumor}_r{w_read}_u{upper}_sort_cal_all.tsv"
        df_cal = pd.read_csv(os.path.join(args.input, cal_file), sep='\t')
        # only consider low onekg freq fusions for ranking
        df_cal_low = df_cal[df_cal['stratum'] == 'low_freq']
        df_cal_low = df_cal_low.sort_values(by='fusion_score', ascending=False).reset_index(drop=True)
        ranks = []
        for idx, c in df_cal_low['fusion_score'].items():
            rank_c = (df_score['fusion_score'] >= c).sum()
            ranks.append(rank_c)
        median_rank = int(pd.Series(ranks).median())
        data['median_rank'].append(median_rank)
        data['w_dna'].append(w_dna)
        data['w_tumor'].append(w_tumor)
        data['w_read'].append(w_read)
        data['upper'].append(upper)

    # make dataframe
    df_results = pd.DataFrame(data)
    # write results to file
    df_results.to_csv(os.path.join(args.input, 'scoring_summary.tsv'), sep='\t')
    return df_results

def plot(df_results, args):
    rows = df_results['w_dna'].nunique()
    columns = df_results['upper'].nunique()
    # squeeze forces 2d ax obj even if rows or columns = 1
    fig,ax = plt.subplots(rows,columns,figsize=(11,7), squeeze=False)
    # get global min/max for color scaling
    rmin = np.log10(df_results['median_rank'].min() + 1)
    rmax = np.log10(df_results['median_rank'].max() + 1)

    # map w_dna and upper to subplot indices
    w_dna_values = sorted(df_results['w_dna'].unique())
    upper_values = sorted(df_results['upper'].unique())
    w_dna_to_i = {v: i for i, v in enumerate(w_dna_values)}
    upper_to_j = {v: j for j, v in enumerate(upper_values)}
    for ((w_dna,upper), df_group) in df_results.groupby(['w_dna','upper']):
        i = w_dna_to_i[w_dna]
        j = upper_to_j[upper]

        pivot = df_group.pivot(columns='w_tumor', index='w_read', values='median_rank')
        # log10 scale the values
        pivot_scaled = np.log10(pivot + 1)  # shift by 1 to avoid log(0)
        cmap = plt.cm.Reds_r
        cmap = mcolors.LinearSegmentedColormap.from_list(
            'Reds_r_light',
            cmap(np.linspace(0.3, 1.0, 256))  # Skip darkest 30%
        )
        im = ax[i,j].imshow(
            pivot_scaled.values,
            cmap=cmap, aspect='auto',
            vmin=rmin, vmax=rmax,
            origin='lower') # reverse colormap so low ranks are highlighted
        # set axes labels only for edge subplots
        # left col
        if j == 0:
            ax[i,j].set_ylabel('$w_{read}$', fontsize=args.fontsize)
            ax[i,j].set_yticks(range(len(pivot.index)), labels=pivot.index)
        else:
            ax[i,j].set_yticks(range(len(pivot.index)), labels=[])
        # bottom row
        if i == rows - 1:
            ax[i,j].set_xlabel('$w_{tumor}$', fontsize=args.fontsize)
            ax[i,j].set_xticks(range(len(pivot.columns)), labels=pivot.columns)
        else:
            ax[i,j].set_xticks(range(len(pivot.columns)), labels=[])
        ax[i,j].set_title(f'$\\mathbf{{w_{{dna}}={w_dna}, upper={upper}}}$', fontweight='bold', fontsize=args.fontsize)
        for (k,l),val in np.ndenumerate(pivot.values):
            text_val = int(val) if not np.isnan(val) else 'N/A'
            ax[i,j].text(l,k,text_val,ha='center',va='center',color='black',fontweight='bold', fontsize=args.fontsize-2)
    # add major axis labels with arrows
    if columns > 1:
        fig.text(0.5, 0.02, 'Upper Factor', ha='center', fontsize=args.fontsize, fontweight='bold')
        # Horizontal arrow at bottom (for columns/upper)
        arrow_x = FancyArrowPatch((0.15, 0.08), (0.9, 0.08),
                            transform=fig.transFigure,
                            arrowstyle='<->', lw=2, color='black',
                            mutation_scale=20)
        fig.add_artist(arrow_x)
    if rows > 1:
        fig.text(0.02, 0.5, '$\\mathbf{w_{dna}}$', va='center', rotation=90, fontsize=args.fontsize, fontweight='bold')
        # Vertical arrow on left (for rows/w_dna)
        arrow_y = FancyArrowPatch((0.08, 0.15), (0.08, 0.9),
                                transform=fig.transFigure,
                                arrowstyle='<->', lw=2, color='black',
                                mutation_scale=20)
        fig.add_artist(arrow_y)
    # colorbar
    # cbar=fig.colorbar(im, ax=ax, orientation='vertical', fraction=.1)
    # cbar.set_ticks([])
    # cbar.set_label('Median Rank')
    # cbar.ax.invert_yaxis()  # Flip colorbar direction
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.savefig(os.path.join(args.input, 'scoring_heatmap.png'))
    return os.path.join(args.input, 'scoring_heatmap.png')

def main():
    df_results = cache_check(args)
    if df_results is not None:
        print("Using cached results")
    else:
        df_results = data_gather(args)
    print("Plotting results")
    outfile = plot(df_results, args)
    print(f"Done. Output saved to {outfile}")

if __name__ == "__main__":
    main()
