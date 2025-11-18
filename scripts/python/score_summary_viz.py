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
print("Exiting: need to make this script faster")
sys.exit(1)

PARAMS = {
    'w_dna': [0.1, 0.5, 0.9],
    'w_tumor': [0.1, 0.5, 0.9],
    'w_read': [0.1, 0.5, 0.9],
    'upper': [10,50,100]
}

parser = argparse.ArgumentParser(description='Visualize fusion scoring')
parser.add_argument('-i', '--input', type=str, required=True, help='Input directory of scoring experiment')
parser.add_argument('-c', '--calibration', type=str, required=True, help='Calibration fusion list')
parser.add_argument('-p', '--processes', type=int, default=5, help='Number of processes to use')
parser.add_argument('--cache', type=str, default=None, help='Use cached summary file if exists')
parser.add_argument('--fontsize', type=int, default=12, help='Font size for plots')
args = parser.parse_args()

def rank_finder(x,y, scorefile):
    # maybe use score instead of x,y to speed up
    cmd = f"awk '$1 == \"{x}\" && $2 == \"{y}\" {{ print NR; exit }}' {scorefile}"
    print(f"Running command: {cmd}")
    # launch cmd using subprocess
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    rank_str = result.stdout.strip()
    if rank_str:
        return int(rank_str) - 1  # minus 1 for header
    else:
        return -1

# check for cached results
if args.cache:
    if os.path.exists(args.cache):
        df_results = pd.read_csv(args.cache, sep='\t')
        print(f"Loaded cached results from {args.cache}")
    else:
        raise ValueError(f"Cache file {args.cache} not found")


if not args.cache:
    # read calibration fusions
    df_cal = pd.read_csv(args.calibration, sep='\t')
    # gather all score files
    scorefiles = []
    for file in os.listdir(args.input):
        if file.endswith('_sort.tsv'):
            scorefile = os.path.join(args.input, file)
            scorefiles.append(scorefile)

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
        # extract params from filename
        match = re.search(r'dna([^_]*)_t([^_]*)_r([^_]*)_u([^_]*)', scorefile)
        if match:
            w_dna = float(match.group(1))
            w_tumor = float(match.group(2))
            w_read = float(match.group(3))
            upper = int(match.group(4))
        print("extracted params:", "w_dna:", w_dna, "w_tumor:", w_tumor, "w_read:", w_read, "upper:", upper)
        outfile = os.path.basename(scorefile).replace('_sort.tsv', '_cal_rank.tsv')
        assert outfile != scorefile, "Output file should be different from scorefile"
        ranks = []
        # make a deep copy of df_cal
        # we will assign ranks for each scorefile and write
        # the calibration rankings
        df_cal_cp = df_cal.copy(deep=True)
        for idx, row in df_cal_cp.iterrows():
            left = row['left']
            right = row['right']
            rank = rank_finder(left, right, scorefile)
            if rank == -1:
                print(f"Warning: Fusion ({left}, {right}) not found in scorefile {scorefile}")
                continue
            df_cal_cp.at[idx, 'rank'] = rank
            ranks.append(rank)
        # write calibration ranks to file
        df_cal_cp.to_csv(os.path.join(args.input, outfile), sep='\t', index=False)
        # get median rank for this scorefile
        ranks = df_cal_cp['rank'].tolist()
        median_rank = pd.Series(ranks).median()
        data['w_dna'].append(w_dna)
        data['w_tumor'].append(w_tumor)
        data['w_read'].append(w_read)
        data['upper'].append(upper)
        data['median_rank'].append(median_rank)

    # make dataframe
    df_results = pd.DataFrame(data)
    # write results to file
    df_results.to_csv(os.path.join(args.input, 'scoring_summary.tsv'), sep='\t')
else:
    df_results = pd.read_csv(args.cache, sep='\t')

# do heatmap visualization
rows = df_results['w_dna'].nunique()
columns = df_results['upper'].nunique()
# squeeze forces 2d array even if rows or columns = 1
fig,ax = plt.subplots(rows,columns,figsize=(11,7), squeeze=False)

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
        ax[i,j].set_ylabel('w_read', fontsize=args.fontsize)
    # bottom row
    if i == rows - 1:
        ax[i,j].set_xlabel('w_tumor', fontsize=args.fontsize)
    ax[i,j].set_title(f'w_dna={w_dna}, upper={upper}')
    ax[i,j].set_xticks(range(len(pivot.columns)), labels=pivot.columns)
    ax[i,j].set_yticks(range(len(pivot.index)), labels=pivot.index)
    for (k,l),val in np.ndenumerate(pivot.values):
        text_val = int(val) if not np.isnan(val) else 'N/A'
        ax[i,j].text(l,k,text_val,ha='center',va='center',color='black',fontweight='bold', fontsize=args.fontsize-2)
# add colorbar
if columns > 1:
    fig.text(0.5, 0.02, 'Upper Factor', ha='center', fontsize=12, fontweight='bold')
    # Horizontal arrow at bottom (for columns/upper)
    arrow_x = FancyArrowPatch((0.15, 0.08), (0.9, 0.08),
                          transform=fig.transFigure,
                          arrowstyle='<->', lw=2, color='black',
                          mutation_scale=20)
    fig.add_artist(arrow_x)
if rows > 1:
    fig.text(0.02, 0.5, '$w_dna$', va='center', rotation=90, fontsize=12, fontweight='bold')
    # Vertical arrow on left (for rows/w_dna)
    arrow_y = FancyArrowPatch((0.08, 0.15), (0.08, 0.9),
                            transform=fig.transFigure,
                            arrowstyle='<->', lw=2, color='black',
                            mutation_scale=20)
    fig.add_artist(arrow_y)
# cbar=fig.colorbar(im, ax=ax, orientation='vertical', fraction=.1)
# cbar.set_ticks([])
# cbar.set_label('Median Rank')
# cbar.ax.invert_yaxis()  # Flip colorbar direction
plt.subplots_adjust(left=0.15, bottom=0.15)
plt.savefig(os.path.join(args.input, 'scoring_heatmap.png'))


#    w_dna w_tumor w_read upper  median_rank
# 0    1.0     0.1    0.5   200         -1.0
# 1    1.0     0.1    0.5    50         -1.0
# 2    1.0     0.1    0.1   200         -1.0
# 3    1.0     0.1    0.1   100         -1.0
# 4    1.0     0.5    0.9   200         -1.0
# 5    1.0     0.5    0.5    50         -1.0

