#!/usr/bin/env python3
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from polymerization.analysis import intersect_file2data
parser=argparse.ArgumentParser()
parser.add_argument("--dir_g2f", default="../2026_06-g2f-all_gene_pairs/g2f_out")
parser.add_argument("--output", default="PARGP1_TIMM23-breakpoint-strand.png")
parser.add_argument("--breakpoint_bins",default=10)
parser.add_argument("--jitter_variance", default=500, help="base pair units")
parser.add_argument("--alpha",default=0.3)
parser.add_argument("--pad",default=1000)
# parser.add_argument("--out_breakpoint", default="PARGP1_TIMM23-breakpoint.png")
# parser.add_argument("--out_strand", default="PARGP1_TIMM23-strand.png")
args = parser.parse_args()
args.jitter_variance=float(args.jitter_variance)
args.pad = int(args.pad)
# args.breakpoint_bins=int(args.breakpoint_bins)

GENE_LEFT='PARGP1'
GENE_LEFT_START=46350283 -args.pad
GENE_LEFT_END=46411621
# GENE_LEFT_BINS = np.linspace(GENE_LEFT_START, GENE_LEFT_END,num=args.breakpoint_bins, endpoint=True)
GENE_RIGHT='TIMM23'
GENE_RIGHT_START=51592080
GENE_RIGHT_END=51623365 + args.pad
# GENE_RIGHT_BINS = np.linspace(GENE_RIGHT_START, GENE_RIGHT_END,num=args.breakpoint_bins, endpoint=True)
INTERSECT_SUFFIX='.giggle.clean.swap.intersect.bed.gz'
NORMAL_PATTERNS=['normal', '1000g']
STRAND_COLORS={'+-': 'blue', '++': 'purple', '--': 'green', '-+': 'orange'}
STRAND_NUMS2SYM={'11':'++', '1-1': '+-', '-11': '-+', '-1-1': '--'}
MODALITY2TITLE={
    'tumor_dna': 'Tumor DNA', 'tumor_rna': 'Tumor_RNA',
    'normal_dna': 'Normal DNA', 'normal_rna': 'Normal_RNA'
}
# STRAND2SYM={'strand_pp': "++", "strand_pn": "+-", "strand_np": "-+", "strand_nn": "--"}


def modality_dict():
    dct = {
        "++": 0,
        "+-": 0,
        "-+": 0,
        "--": 0,
        "reads": []
    }
    return dct

def main():
    data = {
        "tumor_dna": modality_dict(),
        "tumor_rna": modality_dict(),
        "normal_dna": modality_dict(),
        "normal_rna": modality_dict()
    }
    for dir_sub in os.listdir(args.dir_g2f):
        # modality
        modality = ''
        if ('tumor' in dir_sub) and ('dna' in dir_sub):
            modality='tumor_dna'
        if ('tumor' in dir_sub) and ('rna' in dir_sub):
            modality='tumor_rna'
        if (any([pattern in dir_sub for pattern in NORMAL_PATTERNS])) and ('dna' in dir_sub):
            modality='normal_dna'
        if (any([pattern in dir_sub for pattern in NORMAL_PATTERNS])) and ('rna' in dir_sub):
            modality='normal_rna'
        # read
        intersect_file = os.path.join(args.dir_g2f,dir_sub,f"{GENE_LEFT}{INTERSECT_SUFFIX}")
        df = intersect_file2data(
            path=intersect_file,
            gene_right_query=GENE_RIGHT,
            group_by_sample=False,
            bgzip=True,
        )
        # strand
        for group, df_group in df.groupby(['strand_left','strand_right']):
            if group == (1,1):
                data[modality]['++']+=df_group.shape[0]
            if group == (1,-1):
                data[modality]['+-']+=df_group.shape[0]
            if group == (-1,1):
                data[modality]['-+']+=df_group.shape[0]
            if group == (-1,-1):
                data[modality]['--']+=df_group.shape[0]
        # reads
        df['left_breakpoint'] = (df['start_left'] + df['end_left']) // 2
        df['right_breakpoint'] = (df['start_right'] + df['end_right']) // 2
        for i,row in df.iterrows():
            data[modality]['reads'].append(
                (row['left_breakpoint'], row['right_breakpoint'], f'{row['strand_left']}{row['strand_right']}')
            )
    fig,ax = plt.subplots(nrows=2,ncols=4,figsize=(12,4))
    for i,(modality,sub_dict) in enumerate(data.items()):
        bar_labels=[]
        bar_values=[]
        for k,v in sub_dict.items():
            if k=='reads':
                continue
            else:
                bar_labels.append(k)
                bar_values.append(v)
        df_bar = pd.DataFrame(list(zip(bar_labels,bar_values)),columns=['strand', 'count'])
        df_bar = df_bar.sort_values('strand')
        df_bar['color'] = df_bar['strand'].map(STRAND_COLORS)
        ax[0,i].bar(df_bar['strand'],df_bar['count'],color=df_bar['color'])
        ax[0,i].spines['top'].set_visible(False); ax[0,i].spines['right'].set_visible(False)
        ax[0,i].set_title(MODALITY2TITLE[modality], loc="left")
        ax[0,i].set_ylabel("Count")
        starts=[]
        ends=[]
        strand_syms=[]
        for start,end,strands in sub_dict['reads']:
            starts.append(start)
            ends.append(end)
            strand_syms.append(STRAND_NUMS2SYM[strands])
        # jitter
        starts = starts + np.random.normal(0,args.jitter_variance,size=len(starts))
        ends = ends + np.random.normal(0,args.jitter_variance,size=len(starts))
        colors = [STRAND_COLORS[strands] for strands in strand_syms]
        ax[1,i].scatter(starts,ends,color=colors,s=3,alpha=args.alpha)
        ax[1,i].spines['top'].set_visible(False); ax[1,i].spines['right'].set_visible(False)
        ax[1,i].set_xlim(GENE_LEFT_START, GENE_LEFT_END); ax[1,i].set_ylim(GENE_RIGHT_START,GENE_RIGHT_END)
        ax[1,i].set_xlabel(GENE_LEFT); ax[1,i].set_ylabel(GENE_RIGHT)
        ax[1,i].set_title(MODALITY2TITLE[modality], loc="left")
    plt.tight_layout(); plt.savefig(args.output)

if __name__ == "__main__":
    main()
