#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import seaborn as sns
parser=argparse.ArgumentParser()
parser.add_argument("--input",default="cosmic_fusions.evidence.tsv")
parser.add_argument("--output",default="cosmic_fusions.scatter.png")
parser.add_argument("--output_data",default="cosmic_fusions.scatter.data.tsv")
parser.add_argument("--bins",default=15)
parser.add_argument("--title", default="COSMIC")
args = parser.parse_args()
args.bins = int(args.bins)

THOUSG_DNA='reads_thousg_high_coverage_dna'
THOUSG_RNA='reads_thousg_mage_short_read_rna'
ANNOTATIONS_SCATTER=[
    ('GOPC', 'ROS1', 0.85, 0.15),
    ('NAB2', 'STAT6', 0.6, 0.9),
    ('CD74', 'NRG1', 0.2, 0.4)
]
ANNOTATION_BAR=('GOPC', 'ROS1', 0.83, 0.15)
PATTERNS_PCAWG_NORMAL_DNA = ['reads', 'pcawg', 'normal', 'dna']

def main():
    df=pd.read_csv(args.input,sep="\t").dropna()
    x=df[THOUSG_DNA]
    y=df[THOUSG_RNA]
    df[['gene_left', 'gene_right', THOUSG_DNA, THOUSG_RNA]].to_csv(args.output_data,sep="\t",index=False)
    pcawg_normal_cols = []
    for col in df.columns:
        matches=[]
        for pattern in PATTERNS_PCAWG_NORMAL_DNA:
            matches.append(pattern in col)
        if all(matches):
            pcawg_normal_cols.append(col)
    df['pcawg_normal_dna_total'] = df[pcawg_normal_cols].sum(axis=1)

    fig,ax=plt.subplots(nrows=1,ncols=2,figsize=(9,4))

    ### scatter thousg
    divider=make_axes_locatable(ax[0])
    ax_histx=divider.append_axes("top",size="20%",pad=0.1,sharex=ax[0])
    ax_histy=divider.append_axes("right",size="20%",pad=0.1,sharey=ax[0])

    ax[0].scatter(x,y,color='black',alpha=0.3,s=1)
    # ax_histx.hist(x,bins=30,color='black')
    bins_x=np.logspace(np.log10(x[x>0].min()),np.log10(x.max()),args.bins)
    ax_histx.hist(x,bins=bins_x,color='black')
    ax_histy.hist(y,bins=args.bins,orientation='horizontal',color='black')

    for a in ANNOTATIONS_SCATTER:
        ax[0].text(
            a[2], a[3], '--'.join(a[:2]), ha='left', va='top', fontsize=6,
            transform=ax[0].transAxes
        )

    ax[0].set_xlabel("1KG DNA reads")
    ax[0].set_ylabel("1KG RNA reads")

    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax_histx.spines['top'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    ax_histy.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)

    ax[0].set_xscale('log')
    # ax_histx.set_yscale('log')
    ax_histy.set_xscale('log')

    ax_histx.tick_params(axis='x',labelbottom=False)
    ax_histy.tick_params(axis='y',labelleft=False)

    ### pcawg normal hist
    bins=[0,1,10,100,1000,10000,100000,np.inf]
    data=df['pcawg_normal_dna_total']

    counts,_=np.histogram(data,bins=bins)

    ax[1].bar(range(len(counts)),counts,color='black')
    ax[1].set_xticks(range(len(counts)))
    ax[1].set_xticklabels(
        [
            '0','[1,10)','[10,100)',
            '[100,1k)','[1k,10k)','[10k,100k)'
            ,'>100k'],rotation=45,ha='right')

    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].set_xlabel('PCAWG DNA normal reads')
    ax[1].set_ylabel('Fusion count')
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)

    ax[1].text(
        ANNOTATION_BAR[2],
        ANNOTATION_BAR[3],
        '--'.join(ANNOTATION_BAR[:2]),
        ha='left',
        va='top',
        fontsize=6,
        transform=ax[1].transAxes
    )

    fig.suptitle(args.title,x=0.12,ha='center')
    plt.tight_layout()
    plt.savefig(args.output)

if __name__=="__main__":
    main()