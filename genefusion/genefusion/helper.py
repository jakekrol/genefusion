import os,sys
import pandas as pd
import numpy as np

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

def pad_mat(df, genes):
    d = df.copy()
    index = set(df.index)
    new = genes.difference(index)
    pad = pd.DataFrame(columns = df.columns)
    for g in new:
        pad.loc[g] = np.zeros(d.shape[1], dtype=np.int64)
    d = pd.concat([d,pad])
    return d
