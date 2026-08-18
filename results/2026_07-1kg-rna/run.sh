#!/usr/bin/env bash

HOST="ftp.sra.ebi.ac.uk"
DIROUT="output"
mkdir -p $DIROUT
# get 1000g rna ftp data
python3 -c 'from jkbiolib.datasets.loaders import thousg_rna_short_read_samples; df=thousg_rna_short_read_samples(); df.to_csv("thousg_rna_short_read_samples.tsv", index=False, sep="\t")'

# ftp paths
tail -n +2 thousg_rna_short_read_samples.tsv | \
    cut -f 1 | \
    cut -d "/" -f 4- > fastq_paths.txt

sed -i 's|^|\/|' fastq_paths.txt

# outpaths
tail -n +2 thousg_rna_short_read_samples.tsv | \
    cut -f 1 | \
    cut -d "/" -f 9 > outnames.txt
sed -i "s|^|${DIROUT}/|" outnames.txt

# stack columns into tsv
paste fastq_paths.txt outnames.txt <(tail -n +2 thousg_rna_short_read_samples.tsv | cut -f 2) > ftp_queries.tsv

# add host column
sed -i "s|^|${HOST}\t|" ftp_queries.tsv

# header
sed -i '1i\ftp_host\tftp_path\tout_path\tmd5' ftp_queries.tsv

# cleanup
rm fastq_paths.txt outnames.txt


ftp_reliable \
    --query ftp_queries.tsv \
    --logfile ftp_queries.log \
    --cpus 5 \
    --algo md5 \
    --retries 3
