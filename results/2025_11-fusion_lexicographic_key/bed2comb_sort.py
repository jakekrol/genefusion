#!/usr/bin/env python3

# create a 2 column file with gene pairs sorted by genomic position
# lexicographically by chromosome, numerically by start position
# ascending order

import pandas as pd
import sys
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Sort gene pairs by genomic position')
parser.add_argument('-i', '--input', help='Input BED file with gene positions',
                    default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed')
parser.add_argument('-o', '--output', help='Output file for sorted gene pairs',
                    default='gene_pairs_sorted_by_genomic_position.tsv')
args = parser.parse_args()
import pandas as pd


# Read bedfile
df = pd.read_csv(args.input, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene', 'strand'])

# Force types
df['chrom'] = df['chrom'].astype(str)  # String for lexicographic comparison
df['start'] = pd.to_numeric(df['start'])  # Numeric for numeric comparison

# Sort by chromosome (lexicographic) then by start (numeric)
df = df.sort_values(by=['chrom', 'start']).reset_index(drop=True)

genes = df['gene'].tolist()
n = len(genes)
total_pairs = n * (n - 1) // 2

# Generate pairs: i < j in sorted order means gene[i] is "smaller"
with open(args.output, 'w') as f:
    f.write("left\tright\n")
    with tqdm(total=total_pairs, desc="Generating pairs", unit="pairs") as pbar:
        for i in range(n):
            for j in range(i + 1, n):
                f.write(f"{genes[i]}\t{genes[j]}\n")
                pbar.update(1)

print(f"Generated {total_pairs} gene pairs")
print(f"Output written to {args.output}")
