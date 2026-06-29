#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description='Plot a depth versus num. fusions')
parser.add_argument('--input', default='./fusion_counts.tsv', help='Input TSV file')
parser.add_argument('--output', default='./depth_v_fusion_calls.png', help='Output PNG file')
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=None)
df.columns = ['Depth', 'Fusions']
fig, ax = plt.subplots()
ax.plot(df['Depth'], df['Fusions'], marker='o', color = 'black')
ax.set_xlabel('Depth', fontsize=14)
ax.set_ylabel('Number of Fusion Calls', fontsize=14)
ax.set_xticks(df['Depth'])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(args.output)

