# Genefusion

## Install
```
git clone https://github.com/jakekrol/genefusion.git && cd genefusion && pip install -e .
```

## GIGGLE interface

`gf-giggle_sharded` CLI script performs sharded giggle queries
```
# serial
gf-giggle_sharded /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/scratch/2024-12-15-giggle_sharded

# parallel
gf-giggle_sharded /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/scratch/2024-12-15-giggle_sharded True shard True 4
```

- Dependency: `genefusion_giggle.sh` script

- Output file name pattern: `shard_#.chrm.strand.gene.start.end.giggle`

```
# unsharded example
# uses sh script, not py pkg
# construct a tsv input file
# such that each {#} matches the expected input of genefusion_giggle.sh
/data/jake/genefusion/scripts/shell/genefusion_giggle.sh -i alt_sort_b -f /data/jake/genefusion/data/gene_file.txt -c {0} -s {1} -g {2} -l {3} -r {4} -o {5}/{0}.{1}.{2}.{3}.{4}.giggle
```

## STIX interface

```
# parallel over 12 cores
genefusion-stix_sharded \
    /data/jake/genefusion/data/prostate/shards \
    21 39751949 40033704 42836478 42903043 \
    /data/jake/genefusion/scratch/2024-12-15-stix_sharded \
    erg.tmprss2.stix_sharded.stix DEL 12 
```

## DNA gene2gene evidence

### Population

- Input: 1) giggle output file and 2) a BED file of GRCh37 genes
- Output: target regions

```
# pipe or provide giggle files as stdin to gargs
ls | gargs --log ../gargs.log -p 64 -o "$bedtools intersect -a $genefile -b <(cut -f 5-7 {0}) > $outdir/{0}.g2g"
```

### Sample-wise

Same as population, except split GIGGLE output files by the sample column, then do intersect.

## Data

### Gene2gene matrices

loading
```
# load the npz parts as dataframes and stack
def reassemble_dataframe_from_npz(parts_prefix, num_parts):
    dfs = []
    
    for part_num in range(1, num_parts + 1):
        # Load the .npz file
        part_file = f"{parts_prefix}_part{part_num}.npz"
        data = np.load(part_file,allow_pickle=True)
        
        # Convert the dictionary of arrays back to a DataFrame, including the index
        index = data['index']
        df_part = pd.DataFrame({key: data[key] for key in data if key != 'index'}, index=index)
        
        # Append the DataFrame part to the list
        dfs.append(df_part)
    
    # Concatenate all DataFrame parts
    df = pd.concat(dfs)
    
    return df
C = reassemble_dataframe_from_npz("C",10)
```

saving
```
# npz format is much faster to load than a gzipped tsv
def split_dataframe_to_npz(df, output_prefix, num_parts):
    # Calculate the number of rows per part
    rows_per_part = len(df) // num_parts
    
    for part_num in range(num_parts):
        start_row = part_num * rows_per_part
        end_row = (part_num + 1) * rows_per_part if part_num < num_parts - 1 else len(df)
        
        # Slice the DataFrame
        df_part = df.iloc[start_row:end_row]
        
        # Convert the DataFrame to a dictionary of arrays, including the index
        data_dict = {col: df_part[col].values for col in df_part.columns}
        data_dict['index'] = df_part.index.values
        
        # Save the part as an .npz file
        output_file = f"{output_prefix}_part{part_num + 1}.npz"
        np.savez_compressed(output_file, **data_dict)
split_dataframe_to_npz(C, "C", 10)
```



- When comparing two population matrices, shape can differ because some genes lacked any evidence as either source (row) or target (col).
- Resolve with row/col 0 pad

```
# example H padding
g_C = set(C.index)
g_H = set(H.index)
g = g_C.union(g_H)
def pad_mat(df, genes):
    d = df.copy()
    index = set(df.index)
    new = genes.difference(index)
    pad = pd.DataFrame(columns = df.columns)
    for g in new:
        pad.loc[g] = np.zeros(d.shape[1], dtype=np.int64)
    d = pd.concat([d,pad])
    return d
H_padded = pad_mat(H,g)
```
