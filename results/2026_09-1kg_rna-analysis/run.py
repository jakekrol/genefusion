#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
from polymerization.datasets import get_thousg_rna_star_fusion_calls
from polymerization.io import read_bed
import sqlite3

parser = argparse.ArgumentParser()
# parser.add_argument(
#     "--star_fusion_calls",
#     default="../2026_07-1kg-rna/thousg_rna-star_fusion_calls.sort.tsv"
#     help="path to 1000G rna star fusion results")
parser.add_argument(
    "--chimeric_counts",
    default="../2026_07-1kg-rna/chimeric_counts.tsv"
)
parser.add_argument(
    "--bed",
    default="../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed"
)
parser.add_argument(
    "--out_bar_top_star_sample_counts",
    default="topk-star_fusion-sample_count-bar.png"
)
parser.add_argument(
    "--star_bars",
    default=10
)
parser.add_argument(
    "--out_hist_star_samples_per_fusion",
    default="star_fusion-samples_per_fusion-hist.png"
)
parser.add_argument(
    "--out_hist_star_fusions_per_sample",
    default="star_fusion-fusions_per_sample-hist.png"
)
parser.add_argument(
    "--out_scatter_star_chimeric_v_fusions",
    default="star_fusion-chimeric_v_fusions-scatter.png"
)
parser.add_argument(
    "--db",
    default="../2026_09-fusion_db/fusion.db",
    help="path to sqlite db with giggle evidence"
)
parser.add_argument(
    "--mage_sqlite_table",
    default="thousg_mage_short_read_rna",
    help="name of 1kg rna mage table in sqlite db"
)
parser.add_argument(
    "--num_mage_samples_indexed",
    default=767,
    help="number of mage 1kg samples we indexed for our analysis"
)
parser.add_argument(
    "--out_sample_frac_compare_scatter",
    default="compare_sample_frac-scatter.png"
)
parser.add_argument(
    "--out_data_fusions_per_sample",
    default="fusions_per_sample.tsv"
)
parser.add_argument(
    "--out_data_samples_per_fusion",
    default="samples_per_fusion.tsv"
)
args = parser.parse_args()
args.num_mage_samples_indexed = int(args.num_mage_samples_indexed)
assert os.path.exists(args.db)

def main():
    # inputs
    df_bed = read_bed(args.bed, gene_col_idx=3)
    df_chimeric_counts = pd.read_csv(args.chimeric_counts,sep="\t")
    df_thousg_star = get_thousg_rna_star_fusion_calls()
    df_thousg_star['fusion'] = df_thousg_star['gene_left'] + '--' + df_thousg_star['gene_right']

    # fusion set
    thousg_rna_fusion_set = set(df_thousg_star['fusion'].unique())
    thousg_rna_fusion_set_sub = set()
    gene_set = set(df_bed['gene_name'].unique())
    n_fusions=len(thousg_rna_fusion_set)
    for fusion in thousg_rna_fusion_set:
        gene_a, gene_b =tuple(fusion.split('--'))
        if (gene_a in gene_set) and (gene_b in gene_set):
            thousg_rna_fusion_set_sub.add(fusion)
    k_fusions=len(thousg_rna_fusion_set_sub)
    # subset
    mask = df_thousg_star['fusion'].isin(thousg_rna_fusion_set_sub)
    df_thousg_star = df_thousg_star[mask].reset_index(drop=True)

    ### Samples per fusion
    # count samples per fusion
    n_samples = len(df_thousg_star['sample'].unique())
    samples_per_fusion = []
    for fusion, df_group in df_thousg_star.groupby('fusion'):
        samples_per_fusion.append((fusion, len(df_group['sample'].unique()) ))
    df_samples_per_fusion = pd.DataFrame(samples_per_fusion, columns = ['fusion', 'count'])
    df_samples_per_fusion['sample_fraction'] = df_samples_per_fusion['count'] / n_samples
    df_samples_per_fusion = df_samples_per_fusion.sort_values('count', ascending=False).reset_index(drop=True)
    df_samples_per_fusion['gene_left'] = df_samples_per_fusion['fusion'].apply(lambda x: x.split('--')[0])
    df_samples_per_fusion['gene_right'] = df_samples_per_fusion['fusion'].apply(lambda x: x.split('--')[1])

    # get our sample count per fusion
    conn = sqlite3.connect(args.db)
    df_samples_per_fusion.to_sql("temp", conn, if_exists="replace", index=False)
    df_samples_per_fusion = pd.read_sql_query(f"""
        SELECT temp.*, COALESCE({args.mage_sqlite_table}.samples_thousg_mage_short_read_rna,0) AS samples_thousg_mage_short_read_rna
        FROM temp
        LEFT JOIN {args.mage_sqlite_table}
            ON temp.gene_left = {args.mage_sqlite_table}.gene_left
            AND temp.gene_right = {args.mage_sqlite_table}.gene_right
    """, conn)
    df_samples_per_fusion['frac_samples_thousg_mage_short_read_rna'] = \
        df_samples_per_fusion['samples_thousg_mage_short_read_rna'] / args.num_mage_samples_indexed

    # bar most commmon star fusion calls
    bar_data = df_samples_per_fusion.iloc[0:args.star_bars, :]
    bar_data = bar_data.sort_values('count').reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.barh(bar_data['fusion'], bar_data['count'], color='black')
    for i, (count, fraction) in enumerate(zip(bar_data['count'], bar_data['sample_fraction'])):
        ax.text(count, i, f"{fraction*100:.2f}%", va="center")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('1KG RNA MAGE STAR-Fusion', loc='left')
    ax.set_xlabel('Samples per fusion')
    ax.tick_params(axis="y", labelsize=5.5, rotation=60)
    fig.savefig(args.out_bar_top_star_sample_counts)
    plt.close()

    # hist star fusion sample count distribution
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(df_samples_per_fusion['count'], color='black',bins=20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('1KG RNA MAGE STAR-Fusion', loc='left')
    ax.set_xlabel('Samples per fusion')
    ax.set_ylabel('Fusion count')
    ax.set_yscale('log')
    fig.savefig(args.out_hist_star_samples_per_fusion)
    plt.close()

    ### Fusions and chimeric reads per sample
    fusions_per_sample = []
    for sample, df_group in df_thousg_star.groupby('sample'):
        fusions_per_sample.append((sample, len(df_group['fusion'].unique()) ))
    
    df_fusions_per_sample = pd.DataFrame(fusions_per_sample, columns=['sample', 'fusion_count'])

    # hist fusions called per sample
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(df_fusions_per_sample['fusion_count'], color='black',bins=20)
    mean = df_fusions_per_sample['fusion_count'].mean()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(mean,color='red',linestyle='--')
    ax.set_title('1KG RNA MAGE STAR-Fusion', loc='left')
    ax.set_xlabel('Fusions per sample')
    ax.set_ylabel('Sample count')
    # ax.set_yscale('log')
    fig.savefig(args.out_hist_star_fusions_per_sample)
    plt.close()

    # scatter chimeric reads vs fusions called per sample
    fig, ax = plt.subplots(figsize=(6,4))
    samples_intersection = set(df_fusions_per_sample['sample'].unique()).intersection(
        set(df_chimeric_counts['sample'].unique())
    )
    mask = df_fusions_per_sample['sample'].isin(samples_intersection)
    df_fusions_per_sample_sub = df_fusions_per_sample[mask]
    mask = df_chimeric_counts['sample'].isin(samples_intersection)
    df_chimeric_counts_sub = df_chimeric_counts[mask]
    x = df_fusions_per_sample_sub['fusion_count'].values
    y = df_chimeric_counts_sub['chimeric_count'].values
    ax.scatter(x, y, color='black',s=3, alpha = 0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('1KG RNA MAGE STAR-Fusion', loc='left')
    ax.set_xlabel('Fusions per sample')
    ax.set_ylabel('Chimeric reads per sample')
    ax.set_yscale('log')
    ax.set_ylim(1,max(y))
    fig.savefig(args.out_scatter_star_chimeric_v_fusions)
    plt.close()

    ## scatter compare sample fractions
    fig, ax = plt.subplots(figsize=(6,4))
    x = df_samples_per_fusion['sample_fraction']
    y = df_samples_per_fusion['frac_samples_thousg_mage_short_read_rna']
    ax.scatter(x, y, color='black',s=3, alpha = 0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('1KG RNA MAGE STAR-Fusion', loc='left')
    ax.set_xlabel('STAR-Fusion sample fraction')
    ax.set_ylabel('STIX sample fraction')
    fig.savefig(args.out_sample_frac_compare_scatter)
    plt.close()

    ## write data
    df_samples_per_fusion.to_csv(args.out_data_samples_per_fusion,sep="\t",index=False)
    (pd.merge(df_fusions_per_sample, df_chimeric_counts_sub,how='outer', on='sample')).to_csv(
        args.out_data_fusions_per_sample,sep="\t",index=False)
    

if __name__ == "__main__":
    main()