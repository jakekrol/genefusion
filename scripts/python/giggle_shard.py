#!/usr/bin/env python3
import os,sys
import argparse
import re
import math
import shutil
import pandas as pd
import subprocess
import yaml

parser = argparse.ArgumentParser(description='Giggle shard script')
parser.add_argument('-i', '--input', type=str, required=True, help='Input dir with sorted, zipped excord files')
parser.add_argument('-o', '--output', type=str, required=True, help='Output dir parent of shards')
parser.add_argument('-n', '--n_shards', type=int, required=True, help='Number of shards to create')
parser.add_argument('-s', '--shard_meta', type=str, required=True, help='Output file for shard metadata')
parser.add_argument('--pattern', type=str, default=r'.*\.bed\.gz$', help='Pattern to match input files ')
parser.add_argument('--index', type=str, default='sort_b', help='Index name for giggle shards (default: sort_b)')
parser.add_argument('-r', '--run_giggle', action='store_true', help='Run the giggle index command')
args = parser.parse_args()

assert os.path.isdir(args.input), f'Input directory {args.input} does not exist'

# check if giggle is installed by checking return code is 0
if args.run_giggle:
    retcode = subprocess.call(['giggle'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert retcode == 0, 'Giggle is not installed or not in PATH. Please install giggle before running this script.'


# make outputs
os.makedirs(args.output, exist_ok=True)
# make shards
for i in range(args.n_shards):
    os.makedirs(os.path.join(args.output, f'shard_{i}'), exist_ok=True)
    os.makedirs(os.path.join(args.output, f'shard_{i}', 'beds'), exist_ok=True)

pattern = re.compile(args.pattern)
inputs = [f for f in os.listdir(args.input) if pattern.match(f)]
inputs.sort()
inputs = [os.path.join(args.input, f) for f in inputs]
n = len(inputs)
assert n > 0, f'No input files found matching pattern {args.pattern} in {args.input}'

# assign files to shards
shards = [os.path.join(args.output, f'shard_{i}') for i in range(args.n_shards)]
meta = {f'shard_{i}': [] for i in range(args.n_shards)}
while len(inputs) > 0: # until all inputs are assigned
    for i in range(args.n_shards): # loop over shards
        if len(inputs) == 0: # if no more inputs, break
            break
        f = inputs.pop(0) # get the next input file and remove it from the list
        shard_path = os.path.join(shards[i], 'beds', os.path.basename(f))
        shutil.copy(f, shard_path) # assign it to the shard
        print(f'Assigned {f} to {shard_path}')
        meta[f'shard_{i}'] += [shard_path]
# write metadata to yaml
with open(args.shard_meta, 'w') as f:
    yaml.dump(meta, f, default_flow_style=False)
print(f'Wrote shard metadata to {args.shard_meta}')

if args.run_giggle:
    # run giggle index command
    raise NotImplementedError("Giggle indexing is not implemented in this script.")
    for i in range(args.n_shards):
        shard_path = os.path.join(shards[i], 'beds')
        cmd = [
            'giggle', 'index',
            '-i', f"{shard_path}/*",
            '-o', os.path.join(args.output, f'shard_{i}', args.index),
        ]
        print(f'Running command: {" ".join(cmd)}')
        subprocess.run(cmd, cwd=shards[i],check=True)
    print('Giggle indexing completed.')
        








