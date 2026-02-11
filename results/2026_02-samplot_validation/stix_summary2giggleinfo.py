#!/usr/bin/env python

import argparse
import pandas as pd
import re
import yaml
import os

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-d", "--dirstix", 
    default='/data/jake/stix-pcawg-dna')
args = parser.parse_args()

def filename2tissue(x):
    tissues = {
        'blood','esophagus', 'kidney', 'liver', 'ovary'
    }
    for tissue in tissues:
        m = re.search(tissue, x)
        if m:
            return tissue
    raise ValueError('Could not determine tissue from filename: {}'.format(x))
# def filename2genepair(x):
#     y = x.split('.')[3:]
#     y = '.'.join(y)
#     y.replace('.out', '')
#     z = y.split('--')
#     assert len(z) == 2, 'Expected exactly two genes in filename: {}'.format(x)
#     left,right = z
#     return left, right

def tissue2shardmeta(tissue):
    return f'{args.dirstix}/{tissue}_tumor_shard_metadata.yaml'

def sample2shard(sample_in, shardmeta):
    # search the shardmeta dictionary for sample and return cooresponding shard
    for shard, sample_list in shardmeta.items():
        sample_list = [os.path.basename(s) for s in sample_list]
        sample_list = [x.split('.')[0] for x in sample_list]
        for sample in sample_list:
            if re.match(sample_in, sample):
                result = re.sub('data', 'index', shard)
                return result
    raise ValueError('Could not find shard for sample: {}'.format(sample_in))


def main():
    tissue = filename2tissue(args.input)
    shardmeta_file = tissue2shardmeta(tissue)
    # read yaml
    with open(shardmeta_file) as f:
        shardmeta = yaml.safe_load(f)
    df = pd.read_csv(args.input, sep='\t')
    df['filename'] = args.input
    df['tissue'] = tissue
    # df['left_gene'], df['right_gene'] = zip(*df['filename'].apply(filename2genepair))
    # assign shard based on sample name
    df['shard'] = df['sample'].apply(sample2shard, shardmeta=shardmeta)
    # remove filename column
    df = df.drop(columns=['filename'])
    df.to_csv(args.output, sep='\t', index=False)


if __name__ == "__main__":
    main()