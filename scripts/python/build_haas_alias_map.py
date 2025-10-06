#!/usr/bin/env python3
import pandas as pd
import os,sys
import yaml
import argparse

parser = argparse.ArgumentParser(description="Build alias table for Haas 2019 fusions")
parser.add_argument("--input", "-i", type=str,
                    default="../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt", help="Input Haas fusion table")
parser.add_argument("--output", "-o", type=str,
                    default="haas_aliases.yaml", help="Output YAML file name")
args = parser.parse_args()

df=pd.read_csv(args.input,sep='\t', usecols=[2,5,6], encoding='latin1')

def f(x):
    return x.upper() if type(x) == str else str(x).upper()

# gather gene keys
data = {}
for i,row in df.iterrows():
    # fusion column
    fusion = row['fusion']
    a,b = fusion.split('--')
    for j in [a,b]:
        if "," in j:
            c = j.split(',')
            for k in c:
                data[f(k)] = set()
        else:
            data[f(j)] = set()
    # gencode columns
    if pd.notna(row.get("mapped_gencode_A_gene_list")):
        for d in row["mapped_gencode_A_gene_list"].split(","):
            data[f(d)] = set()
    if pd.notna(row.get("mapped_gencode_B_gene_list")):
        for e in row["mapped_gencode_B_gene_list"].split(","):
            data[f(e)] = set()

# build alias map
for i,row in df.iterrows():
    fusion = row['fusion']
    a,b = fusion.split('--')
    # left,right
    for j in [a,b]:
        # aliases
        if "," in j:
            c = j.split(',')
            # for each alias add its synonyms to list
            for k in c:
                z = c.copy()
                z.remove(k)
                for l in z:
                    data[f(k)].add(f(l))
        else:
            # do nothing
            pass
    # gencode columns
    gencode_a = []
    gencode_b = []
    if pd.notna(row.get("mapped_gencode_A_gene_list")):
        gencode_a = [f(d) for d in row["mapped_gencode_A_gene_list"].split(",")]
    if pd.notna(row.get("mapped_gencode_B_gene_list")):
        gencode_b = [f(e) for e in row["mapped_gencode_B_gene_list"].split(",")]
    # add gencode_a genes as aliases to all aliases of a
    if "," in a:
        c = a.split(',')
        for k in c:
            for d in gencode_a:
                data[f(k)].add(d)
    else:
        for d in gencode_a:
            data[f(a)].add(d)
    # add gencode_b genes as aliases to all aliases of b
    if "," in b:
        c = b.split(',')
        for k in c:
            for e in gencode_b:
                data[f(k)].add(e)
    else:
        for e in gencode_b:
            data[f(b)].add(e)

            
for k,v in data.items():
    data[k] = list(data[k])

# rm self from alias list
for k,v in data.items():
    if k in v:
        v.remove(k)
    v.sort()
    data[k] = v



with open(args.output, 'w') as f:
    yaml.dump(data,f)


# for i,val in df['fusion'].items():
#     a,b = val.split('--')
#     for j in [a,b]:
#         if "," in j:
#             c = j.split(',')
#             for k in c:
#                 X.add(f(k))
#         else:
#             X.add(f(j))

# d = {i: set() for i in X}
# for i,val in df['fusion'].items():
#     a,b = val.split('--')
#     # left,right
#     for j in [a,b]:
#         # aliases
#         if "," in j:
#             c = j.split(',')
#             # for each alias add its synonyms to list
#             for k in c:
#                 z = c.copy()
#                 z.remove(k)
#                 for l in z:
#                     d[f(k)].add(f(l))
# for k,v in d.items():
#     d[k] = list(d[k])

# with open(args.output, 'w') as f:
#     yaml.dump(d,f)
    
