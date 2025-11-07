#!/usr/bin/env python3

import os
import argparse
import yaml
import pandas as pd

# score fusions
# with reasonable param search
# visualize score dist
# use calibration fusions as reference

SCRIPT_JOINDDB="join_ddb.py"

### args
parser = argparse.ArgumentParser(description="Score fusions")
parser.add_argument('-i', '--input', required=True, help='Input fusion feature file')
parser.add_argument('-c', '--calibration', required=False, help='Calibration fusions for scoring')
# parser.add_argument('-y', '--yaml',
#                     default='/data/jake/genefusion/data/2025_10-score_yaml_template/score_dna1_t0.5_r0.5_u100.yaml',
#                     help='YAML config file for scoring parameters')
parser.add_argument('-o', '--output_dir', required=True, help='Output directory for experiment')
parser.add_argument('--score_script', default='score_fusions.py', help='Path to scoring script')
parser.add_argument('--size_dna_normal', type=int, default=0, help='Population size for DNA normal')
parser.add_argument('--size_dna_tumor', type=int, default=0, help='Population size for DNA tumor')
parser.add_argument('--size_rna_normal', type=int, default=0, help='Population size for RNA normal')
parser.add_argument('--size_rna_tumor', type=int, default=0, help='Population size for RNA tumor')
parser.add_argument('--size_dna_onekg', type=int, default=2536, help='Population size for DNA 1K Genomes')
args = parser.parse_args()

### setup
def setup():
    # check input columns
    with open(args.input, 'r') as f:
        in_header = f.readline().strip().split('\t')
        print(f"Input columns: {in_header}")
    # check calibration file
    if args.calibration:
        with open(args.calibration, 'r') as f:
            cal_header = f.readline().strip().split('\t')
            print(f"Calibration columns: {cal_header}")
    # script
    assert os.path.exists(args.score_script), f"Scoring script not found: {args.score_script}"
    # yaml
    # assert os.path.exists(args.yaml), f"YAML config file not found: {args.yaml}"
    os.makedirs(args.output_dir, exist_ok=True)
    return in_header

def generate_param_set(
    in_header,
    size_dna_normal,
    size_dna_tumor,
    size_rna_normal,
    size_rna_tumor,
    size_dna_onekg,
    match_read = 'pe'
):
    """Generate a set of parameter combinations for scoring experiments"""
    # list of param dictionaries
    param_list = []
    # read
    rna_norm_read_col = None
    rna_tum_read_col = None
    dna_norm_read_col = None
    dna_tum_read_col = None
    onekg_read_col = None
    # sample
    rna_norm_samp_col = None
    rna_tum_samp_col = None
    dna_norm_samp_col = None
    dna_tum_samp_col = None
    onekg_samp_col = None
    # parse columns
    for col in in_header:
        if ('rna' in col) and ('normal' in col) and (match_read in col):
            rna_norm_read_col = col
        if ('rna' in col) and ('tumor' in col) and (match_read in col):
            rna_tum_read_col = col
        if ('dna' in col) and ('normal' in col) and (match_read in col):
            dna_norm_read_col = col
        if ('dna' in col) and ('tumor' in col) and (match_read in col):
            dna_tum_read_col = col
        if 'onekg' in col and (match_read in col):
            onekg_read_col = col
        if ('rna' in col) and ('normal' in col) and ('sample' in col):
            rna_norm_samp_col = col
        if ('rna' in col) and ('tumor' in col) and ('sample' in col):
            rna_tum_samp_col = col
        if ('dna' in col) and ('normal' in col) and ('sample' in col):
            dna_norm_samp_col = col
        if ('dna' in col) and ('tumor' in col) and ('sample' in col):
            dna_tum_samp_col = col
        if 'onekg' in col and ('sample' in col):
            onekg_samp_col = col
    # hyperparameters to try
    rna_cols = [rna_norm_read_col, rna_tum_read_col, rna_norm_samp_col, rna_tum_samp_col]
    dna_cols = [dna_norm_read_col, dna_tum_read_col, dna_norm_samp_col, dna_tum_samp_col]
    onekg_cols = [onekg_read_col, onekg_samp_col]
    rna = any(col is not None for col in rna_cols)
    dna = any(col is not None for col in dna_cols)
    normal = any(col is not None for col in [rna_norm_read_col, dna_norm_read_col, rna_norm_samp_col, dna_norm_samp_col])
    tumor = any(col is not None for col in [rna_tum_read_col, dna_tum_read_col, rna_tum_samp_col, dna_tum_samp_col])
    onekg = any(col is not None for col in onekg_cols)
    read = any(col is not None for col in [rna_norm_read_col, rna_tum_read_col, dna_norm_read_col, dna_tum_read_col, onekg_read_col])
    sample = any(col is not None for col in [rna_norm_samp_col, rna_tum_samp_col, dna_norm_samp_col, dna_tum_samp_col, onekg_samp_col])
    # weight dna/rna
    if rna and dna:
        w_dna = [0.1,0.5,0.9]
    elif dna and (not rna):
        w_dna = [1.0]
    elif rna and (not dna):
        w_dna = [0.0]
    else:
        raise ValueError("No RNA nor DNA columns found in input")
    # weight tumor/normal
    if normal and tumor:
        w_tumor = [0.1,0.5,0.9]
    elif tumor and (not normal):
        w_tumor = [1.0]
    elif normal and (not tumor):
        w_tumor = [0.0]
    else:
        raise ValueError("No tumor nor normal columns found in input")
    # weight read/sample
    if read and sample:
        w_read = [0.1,0.5,0.9]
    elif read and (not sample):
        w_read = [1.0]
    elif sample and (not read):
        w_read = [0.0]
    else:
        raise ValueError("No read nor sample columns found in input")
    # upper factor
    upper_factors = [50, 100, 200]
    # generate param combinations
    # and print total amount of combos
    n = len(w_dna) * len(w_tumor) * len(w_read) * len(upper_factors)
    print(f"Generating {n} parameter combinations for scoring experiments")
    for w_d in w_dna:
        for w_t in w_tumor:
            for w_r in w_read:
                for u_f in upper_factors:
                    param_dict = {
                        'read_rna_normal': rna_norm_read_col if rna_norm_read_col is not None else 0,
                        'read_rna_tumor': rna_tum_read_col if rna_tum_read_col is not None else 0,
                        'read_dna_normal': dna_norm_read_col if dna_norm_read_col is not None else 0,
                        'read_dna_tumor': dna_tum_read_col if dna_tum_read_col is not None else 0,
                        'read_dna_onekg': onekg_read_col if onekg_read_col is not None else 0,
                        'sample_rna_normal': rna_norm_samp_col if rna_norm_samp_col is not None else 0,
                        'sample_rna_tumor': rna_tum_samp_col if rna_tum_samp_col is not None else 0,
                        'sample_dna_normal': dna_norm_samp_col if dna_norm_samp_col is not None else 0,
                        'sample_dna_tumor': dna_tum_samp_col if dna_tum_samp_col is not None else 0,
                        'sample_dna_onekg': onekg_samp_col if onekg_samp_col is not None else 0,
                        'pop_size_rna_normal': size_rna_normal,
                        'pop_size_rna_tumor': size_rna_tumor,
                        'pop_size_dna_normal': size_dna_normal,
                        'pop_size_dna_tumor': size_dna_tumor,
                        'pop_size_dna_onekg': size_dna_onekg,
                        'weight_dna': w_d,
                        'weight_tumor': w_t,
                        'weight_read': w_r,
                        'upper_factor': u_f
                    }
                    param_list.append(param_dict)
    return param_list, n

def score(X, yml_path, outfile):
    """Run scoring script with given input and yaml, output to specified directory"""
    input_file = X
    cmd = f"python {args.score_script} -i {input_file} -o {outfile} --score_yaml {yml_path}"
    print(f"Running command: {cmd}")
    os.system(cmd)
    return outfile
def sort_tbl(infile, outfile, sort_col='fusion_score', descending=True):
    """Sort a TSV file by specified column"""
    sort_flag = 'DESC' if descending else 'ASC'
    cmd = f"sort_tbl.sh -i {infile} -o {outfile} -c {sort_col} -s {sort_flag}"
    print(f"Running command: {cmd}")
    os.system(cmd)
    return outfile

def score_calibration(cal, scored):
    # do a join ddb
    cmd = f"join_ddb.py -l {cal} -r {scored} --type left -k left,right -o {scored}_cal.tsv"
    print(f"Running command: {cmd}")
    os.system(cmd)
    return f"{scored}_cal.tsv"

def viz(cal, scored):
    ### do a histogram of scores with vlines at calibration points
    # make a csv of calibration scores
    df_cal = pd.read_csv(cal, sep='\t')
    df_cal['x'] = df_cal['fusion_score'].apply(lambda x: f"{x},r,0.25")
    df_cal['x'].to_csv(scored + '_cal_scores.csv', index=False, header=False)
    ### find fusion score column in scored efficiently
    with open(scored, 'r') as f:
        header = f.readline().strip().split('\t')
    if 'fusion_score' in header:
        score_col_index = header.index('fusion_score') + 1  # +1 for 1-based cut
    else:
        raise ValueError("fusion_score column not found in scored file")
    ### extract null distribution
    cmd = f"tail -n +2 {scored} | cut -f {score_col_index} | shuf -n 1000000 > {scored}.null"
    print(f"Running command: {cmd}")
    os.system(cmd)
    ### plot histogram
    cmd = f"cat {scored}.null | hist.py -o {scored}_score_dist.png --bins 30 " \
        f"--axvline {scored}_cal_scores.csv -y 'Frequency' -x 'Fusion score' " \
        f"--ylog --alpha 0.5 " \
        f"--title 'Fusion Score Distribution'"
    print(f"Running command: {cmd}")
    os.system(cmd)
    return f"{scored}_score_dist.png"

def main():
    in_header = setup()
    param_list, n= generate_param_set(
        in_header,
        args.size_dna_normal,
        args.size_dna_tumor,
        args.size_rna_normal,
        args.size_rna_tumor,
        args.size_dna_onekg,
        match_read='pe'
    )
    # write params to yamls
    # use params in filename with pattern score_dna{w_d}_t{w_t}_r{w_r}_u{u_f}.yaml
    # use same pattern for output scored file
    for i, param_dict in enumerate(param_list):
        param_yaml_path = \
            os.path.join(
                args.output_dir,
                f'score_dna{param_dict["weight_dna"]}_t{param_dict["weight_tumor"]}_r{param_dict["weight_read"]}_u{param_dict["upper_factor"]}.yaml'
            )
        with open(param_yaml_path, 'w') as f:
            yaml.dump(param_dict, f)
        print(f"Running scoring with parameter set {i+1}/{n}: {param_yaml_path}")
        scored_output = score(args.input, param_yaml_path,
                              os.path.join(
                                  args.output_dir,
                                  f'scored_dna{param_dict["weight_dna"]}_t{param_dict["weight_tumor"]}_r{param_dict["weight_read"]}_u{param_dict["upper_factor"]}.tsv'
                              ))
        print(f"Scored output written to: {scored_output}")
        # sort
        scored_output_sorted = sort_tbl(scored_output,
                                        scored_output.replace('.tsv', '_sort.tsv'),
                                        sort_col='fusion_score',
                                        descending=True)
        print(f"Sorted scored output written to: {scored_output_sorted}")
        if args.calibration:
            scored_cal = score_calibration(args.calibration, scored_output_sorted)
            print(f"Scored calibration output written to: {scored_cal}")
            outhist = viz(args.calibration, scored_output_sorted)
            print(f"Score distribution plot written to: {outhist}")
    print("All scoring experiments completed.")
            
if __name__ == "__main__":
    main()
