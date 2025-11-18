#!/usr/bin/env python3
import sys

import os
import argparse
import yaml
import pandas as pd
import re
import math

# score fusions
# with reasonable param search
# visualize score dist
# use calibration fusions as reference

# to-do
# add random seed to shuf for reproducibility
# parallelize param runs

SCRIPT_JOINDDB="join_ddb.py"
PARAMS = {
    'w_dna': [0.1, 0.5, 0.9],
    'w_tumor': [0.1, 0.5, 0.9],
    'w_read': [0.1, 0.5, 0.9],
    'upper': [25, 50, 100]
}

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
parser.add_argument('--size_dna_onekg', type=int, default=2535, help='Population size for DNA 1K Genomes')
parser.add_argument('--onekg_sample_col', type=str, default='sample_count_onekg', help='Column name for onekg sample count')
args = parser.parse_args()

onekg_sample_count_thresh = math.ceil(args.size_dna_onekg * 0.01)  # 1% frequency threshold

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
        w_tumor = PARAMS['w_tumor']
    elif tumor and (not normal):
        w_tumor = [1.0]
    elif normal and (not tumor):
        w_tumor = [0.0]
    else:
        raise ValueError("No tumor nor normal columns found in input")
    # weight read/sample
    if read and sample:
        w_read = PARAMS['w_read']
    elif read and (not sample):
        w_read = [1.0]
    elif sample and (not read):
        w_read = [0.0]
    else:
        raise ValueError("No read nor sample columns found in input")
    # upper factor
    upper_factors = PARAMS['upper']
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

def fname2params(fname):
    match = re.search(r'dna([^_]*)_t([^_]*)_r([^_]*)_u([^_]*)', fname)
    if match:
        w_dna = float(match.group(1))
        w_tumor = float(match.group(2))
        w_read = float(match.group(3))
        upper = int(match.group(4))
    assert all(v is not None for v in [w_dna, w_tumor, w_read, upper]), "Failed to extract all parameters from filename"
    return w_dna, w_tumor, w_read, upper

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
    cmd = f"join_ddb.py -l {cal} -r {scored} --type left -k left,right -o {os.path.join(args.output_dir, f'{scored}_cal.tsv')}"
    print(f"Running command: {cmd}")
    os.system(cmd)
    return f"{scored}_cal.tsv"

def viz(cal, scored):
    ### do a boxplot of scores stratified by 
    # i) null
    # ii) calibration >0.01 onekg freq
    # iii) calibration <=0.01 onekg freq

    # find fusion score column in scored efficiently
    with open(scored, 'r') as f:
        header = f.readline().strip().split('\t')
    if 'fusion_score' in header:
        score_col_index = header.index('fusion_score') + 1  # +1 for 1-based cut
    else:
        raise ValueError("fusion_score column not found in scored file")

    # get parameters for title
    w_dna, w_tumor, w_read, upper = fname2params(scored)

    # null
    # --random-source for reproducibility
    cmd = f"tail -n +2 {scored} | cut -f {score_col_index} | shuf -n 1000000 --random-source=<(yes 13) > {scored}.null"
    os.system(cmd)
    
    # calibration
    df_cal = pd.read_csv(cal, sep='\t')
    # stratify
    df_cal['stratum'] = df_cal[args.onekg_sample_col].apply(
        lambda x: 'high_freq' if x >= onekg_sample_count_thresh else 'low_freq'
    )
    assert df_cal['fusion_score'].notna().all(), "NA values found in fusion_score column of calibration data"
    # group by stratum and write to files
    high_freq_file = f"{scored}".replace('.tsv', '_cal_high_freq.tsv')
    low_freq_file = f"{scored}".replace('.tsv', '_cal_low_freq.tsv')
    df_cal[df_cal['stratum'] == 'high_freq']['fusion_score'].to_csv(high_freq_file, index=False, sep='\t', header=False)
    df_cal[df_cal['stratum'] == 'low_freq']['fusion_score'].to_csv(low_freq_file, index=False, sep='\t', header=False)
    # make boxplot
    cmd = f"boxplot.py --files {scored}.null,{high_freq_file},{low_freq_file} " \
        f"-o {scored}_score_boxplot.png " \
        f"--xticklabels 'Rand. gene pairs, Cal. >1% pop. freq 1KGP, Cal. <=1% pop. freq 1KGP' " \
        f"-y 'Fusion score' --title 'w_dna={w_dna};w_tum={w_tumor};w_read={w_read};u={upper}'"
    os.system(cmd)

    # finally write full calibration tsv with fusion scores
    scored_calibration_file = os.path.join(args.output_dir, f"{os.path.basename(scored).replace('.tsv', '_call_all.tsv')}")
    df_cal.to_csv(scored_calibration_file, index=False, sep='\t')

    # ### do a histogram of scores with vlines at calibration points
    # # make a csv of calibration scores
    # df_cal = pd.read_csv(cal, sep='\t')
    # # color by onekg_pop_freq
    # df_cal['color'] = df_cal[args.onekg_sample_col].apply(
    #     lambda x: 'blue' if x >= onekg_sample_count_thresh else 'red'
    # )
    # # drop rows with NA fusion score for now
    # b = df_cal.shape[0]
    # df_cal = df_cal.dropna(subset=['fusion_score'])
    # a = df_cal.shape[0]
    # # write num dropped nas to logfile
    # with open(os.path.join(args.output_dir, 'log.txt'), 'a') as logf:
    #     logf.write(f"Dropped {b - a} calibration fusions with NA fusion_score from {b} total\n")
    # df_cal['alpha'] = 0.5
    # df_cal[['fusion_score','color','alpha']].to_csv(scored + '_cal_scores.csv', index=False, header=False, sep='\t')
    # ### find fusion score column in scored efficiently
    # with open(scored, 'r') as f:
    #     header = f.readline().strip().split('\t')
    # if 'fusion_score' in header:
    #     score_col_index = header.index('fusion_score') + 1  # +1 for 1-based cut
    # else:
    #     raise ValueError("fusion_score column not found in scored file")
    # ### extract null distribution
    # # --ranomdom-source for reproducibility
    # cmd = f"tail -n +2 {scored} | cut -f {score_col_index} | shuf -n 1000000 --random-source=<(yes 13) > {scored}.null"
    # os.system(cmd)
    # # ensure calibration scores are appended
    # df_cal['fusion_score'].to_csv(f"{scored}.null", index=False, header=False, mode='a')
        
    # ### plot histogram
    # # get params for title
    # w_dna, w_tumor, w_read, upper = fname2params(scored)
    # cmd = f"cat {scored}.null | hist.py -o {scored}_score_dist.png --bins 30 " \
    #     f"--axvline {scored}_cal_scores.csv -y 'Frequency' -x 'Fusion score' " \
    #     f"--ylog --alpha 0.5 " \
    #     f"--title f'w_dna={w_dna},w_tumor={w_tumor},w_read={w_read},u={upper}'"
    # print(f"Running command: {cmd}")
    # os.system(cmd)
    # return f"{scored}_score_dist.png"

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
        score_outfile = os.path.join(
            args.output_dir,
            f'scored_dna{param_dict["weight_dna"]}_t{param_dict["weight_tumor"]}_r{param_dict["weight_read"]}_u{param_dict["upper_factor"]}.tsv'
        )
        # cache check
        if os.path.exists(score_outfile):
            print(f"Scored output already exists, skipping: {score_outfile}")
        else:
            # score
            score(args.input, param_yaml_path, score_outfile)
        print(f"Scored output written to: {score_outfile}")
        # sort
        scored_output_sorted = sort_tbl(score_outfile,
                                        score_outfile.replace('.tsv', '_sort.tsv'),
                                        sort_col='fusion_score',
                                        descending=True)
        print(f"Sorted scored output written to: {scored_output_sorted}")
        # viz 
        if args.calibration:
            scored_cal = score_calibration(args.calibration, scored_output_sorted)
            print(f"Scored calibration output written to: {scored_cal}")
            outhist = viz(scored_cal, scored_output_sorted)
            print(f"Score distribution plot written to: {outhist}")
    print("All scoring experiments completed.")
            
if __name__ == "__main__":
    main()
