#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
parser=argparse.ArgumentParser()
parser.add_argument("--input",default="cosmic_fusions.evidence.tsv")
parser.add_argument("--output",default="cosmic_fusions.scatter.png")
parser.add_argument("--bins",default=10)
args = parser.parse_args()
args.bins = int(args.bins)

THOUSG_DNA='reads_thousg_high_coverage_dna'
THOUSG_RNA='reads_thousg_mage_short_read_rna'

def main():
    df=pd.read_csv(args.input,sep="\t").dropna()
    x=df[THOUSG_DNA]
    y=df[THOUSG_RNA]

    fig,ax=plt.subplots(figsize=(6,4))
    divider=make_axes_locatable(ax)
    ax_histx=divider.append_axes("top",size="20%",pad=0.1,sharex=ax)
    ax_histy=divider.append_axes("right",size="20%",pad=0.1,sharey=ax)

    ax.scatter(x,y,color='black',alpha=0.3,s=1)
    ax_histx.hist(x,bins=30,color='black')
    ax_histy.hist(y,bins=30,orientation='horizontal',color='black')

    ax.set_xlabel("1KG DNA reads")
    ax.set_ylabel("1KG RNA reads")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax_histx.spines['top'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    ax_histy.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)

    ax.set_xscale('log')
    ax_histx.set_yscale('log')
    ax_histy.set_xscale('log')

    ax_histx.tick_params(axis='x',labelbottom=False)
    ax_histy.tick_params(axis='y',labelleft=False)
    plt.tight_layout()
    plt.savefig(args.output)

if __name__=="__main__":
    main()