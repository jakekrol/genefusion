#!/usr/bin/env python3
import argparse
from collections import defaultdict
import os
import pandas as pd
from polymerization.stix2fusion import left_sort_fusion_set
from polymerization.io import read_bed
parser=argparse.ArgumentParser()
parser.add_argument("--phenotypes", default="../../data/2025_03-cosmic/Cosmic_Classification_v105_GRCh37.tsv.gz")
parser.add_argument("--fusions", default="../../data/2025_03-cosmic/Cosmic_Fusion_v105_GRCh37.tsv.gz")
parser.add_argument("--output",default="cosmic_fusions.tsv")
parser.add_argument("--output_subset",default="cosmic_fusions.subset.tsv")
parser.add_argument("--bed", default="../2026_09-living_bed/grch37.genes.sort.bed")
args = parser.parse_args()

TISSUES={'blood', 'kidney', 'liver', 'ovary'}
# classification columns
COL_PRIMARY_SITE='PRIMARY_SITE'
COL_SITE_SUBTYPE1='SITE_SUBTYPE_1'
COL_PHENOTYPE_ID='COSMIC_PHENOTYPE_ID' # shared
# fusion columns
COL_FUSION_ID='COSMIC_FUSION_ID'
COL_5_PRIME_GENE='FIVE_PRIME_GENE_SYMBOL'
COL_3_PRIME_GENE='THREE_PRIME_GENE_SYMBOL'
SITE_BLOOD_LYMPH='haematopoietic_and_lymphoid_tissue'
SITE_SUB1_BLOOD='blood_vessel'
GENES_EXCLUDE={
    'JPT1',
    'LHFPL6',
    'WDCP',
    'DUX4L1',
    'SHTN1',
    'RELCH',
    'AFDN',
    'CIP2A',
    'KNL1',
    'CASP8AP2', # lookupable, should add to bed
    'ERO1A'
}

def set_blood_tissue(primary_site: str, site_subtype: str):
    if (primary_site == SITE_BLOOD_LYMPH) and (site_subtype == SITE_SUB1_BLOOD):
        return 'blood'
    else:
        return primary_site

def main():
    df_bed=read_bed(args.bed,gene_col_idx=3)
    df_fusion = pd.read_csv(args.fusions,sep="\t")
    df_pheno = pd.read_csv(args.phenotypes,sep="\t")
    df_pheno[COL_PRIMARY_SITE] = df_pheno.apply(
        lambda row: set_blood_tissue(row[COL_PRIMARY_SITE], row[COL_SITE_SUBTYPE1]), axis=1
    )
    mask = ~df_fusion[COL_5_PRIME_GENE].isin(GENES_EXCLUDE)
    df_fusion = df_fusion[mask]
    mask = ~df_fusion[COL_3_PRIME_GENE].isin(GENES_EXCLUDE)
    df_fusion = df_fusion[mask]
    ### all
    df_fusion_all = df_fusion.drop_duplicates(subset=[COL_5_PRIME_GENE, COL_3_PRIME_GENE])
    df_fusion_all.dropna(subset=[COL_5_PRIME_GENE, COL_3_PRIME_GENE], inplace=True)
    df_fusion_all.reset_index(drop=True,inplace=True)
    df_fusion_all=df_fusion_all[[COL_5_PRIME_GENE, COL_3_PRIME_GENE]]
    df_fusion_all = left_sort_fusion_set(df_fusion_all,df_bed)
    df_fusion_all.drop(columns=['gene_x','gene_y'],inplace=True)
    df_fusion_all['fusion'] = df_fusion_all.apply(
        lambda row: ''.join(sorted(row['gene_left'] + '--' + row['gene_right'])), axis=1
    )
    df_fusion_all.drop_duplicates(subset=['fusion'],inplace=True)
    df_fusion_all.drop(columns=['fusion'],inplace=True)
    df_fusion_all.reset_index(drop=True)
    mask = df_fusion_all['gene_left'] != df_fusion_all['gene_right']
    df_fusion_all = df_fusion_all[mask]
    df_fusion_all.to_csv(args.output,sep='\t',index=False)

    ### subset
    # fusion to phenotypes
    fusion2phenos=defaultdict(set)
    for i,row in df_fusion.iterrows():
        fusion=row[COL_5_PRIME_GENE] + '--' + row[COL_3_PRIME_GENE]
        pheno=row[COL_PHENOTYPE_ID]
        fusion2phenos[fusion].add(pheno)
    # phenotype to tissue
    pheno2tissue=dict()
    for i,row in df_pheno.iterrows():
        pheno=row[COL_PHENOTYPE_ID]
        tissue=row[COL_PRIMARY_SITE]
        pheno2tissue[pheno]=tissue
    # fusion to tissue
    fusion2tissues=defaultdict(set)
    for fusion,phenos in fusion2phenos.items():
        for pheno in phenos:
            tissue = pheno2tissue[pheno]
            fusion2tissues[fusion].add(tissue)
    n=0
    fusion_subset=set()
    for fusion, tissues in fusion2tissues.items():
        for tissue in tissues:
            if tissue in TISSUES:
                n+=1
                fusion_subset.add(fusion)
    df_fusion['fusion'] = df_fusion[COL_5_PRIME_GENE] + '--' + df_fusion[COL_3_PRIME_GENE]
    data = []
    for fusion in fusion_subset:
        gene_left = fusion.split('--')[0]
        gene_right = fusion.split('--')[1]
        tissues=','.join(list(fusion2tissues[fusion].intersection(TISSUES)))
        data.append((gene_left,gene_right,tissues))
    df_out=pd.DataFrame(data,columns=['gene_x','gene_y', 'tissue'])
    df_out_sort = left_sort_fusion_set(df_out,df_bed)
    df_out_sort.drop_duplicates(inplace=True)
    df_out_sort.to_csv(args.output_subset,sep='\t',index=False)

if __name__=="__main__":
    main()
