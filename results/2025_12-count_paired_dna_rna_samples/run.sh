#!/usr/bin/env bash
set -euo pipefail
t_0=$(date +%s)

join_script='/data/jake/rl-tools/wrangle/join_ddb.py'
mapping_file="../2025_12-pcawg_donor_fid_map/icgc25k-wgs_rna-file_id-file_name_donor_map.tsv"
data_dir=../../data/dna_pad
outfile="paired_dna_rna_samples.tsv"
intermediate="tissue_dna_fids.tsv"

printf "tissue\tfile_id\n" > "${intermediate}"
tissues=(blood esophagus kidney liver ovary)
for tissue in "${tissues[@]}"; do
    echo "# processing tissue: ${tissue}"
    sample_yaml="${data_dir}/${tissue}/shards.yaml"
    mapfile -t fids < <(grep "^- " ${sample_yaml} | sed 's|^- ||')
    for fid in "${fids[@]}"; do
        fid=$(basename $fid)
        fid=$(echo "$fid" | cut -d'.' -f1)
        printf "%s\t%s\n" "${tissue}" "${fid}" >> "${intermediate}"
    done
done

echo "# joining to get donor ids"
$join_script \
    --left "${intermediate}" \
    --right "${mapping_file}" \
    --type left \
    --key file_id \
    --output "${outfile}"
echo "# output written to ${outfile}"

# # join with donor ids to find paired samples
$join_script \
    --left <(cut -f 1,3 $outfile) \
    --right ${mapping_file} \
    --type left \
    --key donor_id \
    --output "${outfile}"


echo "# split by tissue"
tail -n +2 "${outfile}" | \
    awk -F'\t' '{print > $1"_samples.tsv"}'
script="
import pandas as pd
files = ['blood_samples.tsv', 'esophagus_samples.tsv', 'kidney_samples.tsv', 'liver_samples.tsv', 'ovary_samples.tsv']
for f in files:
    df = pd.read_csv(f, sep='\t', header=None)
    df.columns = ['tissue', 'donor_id', 'file_id', 'file_name', 'specimen_type','project_code', 'data_type']
    df = df[['donor_id', 'data_type', 'specimen_type']]
    df['specimen_type'] = df['specimen_type'].apply(lambda x: 'normal' if 'normal' in x.lower() else x)
    df['specimen_type'] = df['specimen_type'].apply(lambda x: 'tumor' if 'tumour' in x.lower() else x)
    df_counts = df.groupby(['donor_id', 'data_type'])['specimen_type'].value_counts().unstack(fill_value=0)
    df_counts.reset_index(inplace=True)
    # ignore normal
    df_counts.drop(columns=['normal'], inplace=True)
    # binarize tumor counts to 0/1
    df_counts['tumor'] = df_counts['tumor'].apply(lambda x: 1 if x >= 1 else 0)
    # count cases where both WGS and RNA-seq have tumor == 1
    paired_counts = 0
    for donor_id, group in df_counts.groupby('donor_id'):
        if ('RNA-Seq' in group['data_type'].values) and ('WGS' in group['data_type'].values):
            if group['tumor'].sum() == 2:
                paired_counts += 1
        else:
            continue
    print('paired counts in {}: {}'.format(f, paired_counts))
"
python3 -c "$script"

# sort main file
sort -k1,1 -k2,2 "${outfile}" > tmp && mv tmp "${outfile}"

echo "# time elapsed: $(( $(date +%s) - t_0 )) seconds"
