#!/usr/bin/env python3

import yaml
import matplotlib.pyplot as plt
import argparse
import os,sys
import glob
import pandas as pd

parser = argparse.ArgumentParser(description='pcawg sample counts')
parser.add_argument('-d', '--dna_index_dir', help='dna index directory', required=True)
parser.add_argument('-r', '--rna_index_dir', help='rna index directory', required=True)
parser.add_argument('-o', '--output_dir', help='output directory', required=True)

def yamlfname2data(x):
    tissue = x.split('_')[0]
    specimen = x.split('_')[1]
    return tissue, specimen

def main():
    args = parser.parse_args()
    plot_data = {'tissue': [], 'specimen': [], 'samplecount': [], 'modality': []}

    yamls_dna = glob.glob(os.path.join(args.dna_index_dir, '*.yaml'))
    yamls_dna = [os.path.basename(x) for x in yamls_dna]
    for y in yamls_dna:
        tissue, specimen = yamlfname2data(y)
        plot_data['tissue'].append(tissue)
        plot_data['specimen'].append(specimen)
        plot_data['modality'].append('dna')
        samples=0
        with open(os.path.join(args.dna_index_dir, y)) as f:
            data = yaml.safe_load(f)
            for k, v in data.items():
                samples += len(v)
        plot_data['samplecount'].append(samples)
    yamls_rna = glob.glob(os.path.join(args.rna_index_dir, '*.yaml'))
    yamls_rna = [os.path.basename(x) for x in yamls_rna]
    for y in yamls_rna:
        tissue, specimen = yamlfname2data(y)
        plot_data['tissue'].append(tissue)
        plot_data['specimen'].append(specimen)
        plot_data['modality'].append('rna')
        samples=0
        with open(os.path.join(args.rna_index_dir, y)) as f:
            data = yaml.safe_load(f)
            for k, v in data.items():
                samples += len(v)
        plot_data['samplecount'].append(samples)
    df = pd.DataFrame(plot_data)
    df = df[['tissue', 'specimen', 'modality', 'samplecount']]
    df = df.sort_values(by=['tissue', 'specimen', 'modality', 'samplecount'])
    df.to_csv(os.path.join(args.output_dir, 'sample_counts.tsv'), sep='\t', index=False)
    total_samples = df['samplecount'].sum()
    with open(os.path.join(args.output_dir, 'total_sample_count.txt'), 'w') as f:
        f.write(f'Total sample count: {total_samples}\n')


    


if __name__ == "__main__":
    main()
