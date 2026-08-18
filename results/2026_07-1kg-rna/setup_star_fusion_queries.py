#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('--dir_fastq', default='output', help='Directory containing FASTQ files')
parser.add_argument('--output', default='star_fusion_queries.tsv', help='Output file path for the STAR-Fusion queries')
args = parser.parse_args()

def main():
    files = os.listdir(args.dir_fastq)
    outdata = []
    df = pd.DataFrame({'fastq_file': files})
    df['sample'] = df['fastq_file'].apply(lambda x: x.split('_')[0])
    for sample, df_group in df.groupby('sample'):
        read_1 = os.path.join(args.dir_fastq, f'{sample}_1.fastq.gz')
        read_2 = os.path.join(args.dir_fastq, f'{sample}_2.fastq.gz')
        # make sure both read files exist
        read_1_exists = os.path.exists(read_1)
        read_2_exists = os.path.exists(read_2)
        if read_1_exists and read_2_exists:
            outdata.append({'sample': sample, 'read_1': read_1, 'read_2': read_2})
    df_out = pd.DataFrame(outdata)
    df_out.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()