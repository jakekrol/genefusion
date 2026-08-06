#!/usr/bin/env python3
from jkbiolib.datasets.loaders import *
from ftplib import FTP
import multiprocessing as mp
import os
import sys
import hashlib

CPUS=10
OUT_SR='1kg_short_read'
OUT_LR='1kg_long_read'
HOSTNAME='ftp.sra.ebi.ac.uk'
LOGFILE='download.log'
VERIFIED_DOWNLOADS='verified_downloads.tsv'

def downloader(ftp_host, ftp_path, out_path):
    ftp = FTP(ftp_host)
    ftp.login()
    print("# downloading {} from {} to {}".format(ftp_path, ftp_host, out_path))
    with open(out_path, 'wb') as f:
        try:
            ftp.retrbinary(f'RETR {ftp_path}', f.write)
        except Exception as e:
            print("# error downloading {} from {} to {}: {}".format(ftp_path, ftp_host, out_path, e))

def get_checksum(filepath, algo='md5'):
    h = hashlib.new(algo)
    with open (filepath, 'rb') as f:
        # 4096 is a common buffer size for reading files in chunks
        # could be any power of 2
        for chunk in iter(lambda: f.read(4096), b""):
            # hash bytes
            h.update(chunk)
    # get string representation of the hash
    return h.hexdigest()

def main():
    os.makedirs(OUT_SR, exist_ok=True)
    os.makedirs(OUT_LR, exist_ok=True)
    df_1kg_sr_rna = thousg_rna_short_read_samples()
    df_1kg_sr_rna['hostname'] = df_1kg_sr_rna['url'].apply(lambda x: x.split('//')[1].split('/')[0])
    df_1kg_sr_rna['ftp_path'] = df_1kg_sr_rna['url'].apply(lambda x: '/' + '/'.join(x.split('//')[1].split('/')[1:]))
    df_1kg_sr_rna['read_index'] = df_1kg_sr_rna['url'].apply(lambda x: os.path.basename(x).split('_')[1].split('.')[0])
    df_1kg_sr_rna['outfile'] = df_1kg_sr_rna.apply(lambda row: os.path.join(OUT_SR, row['Sample'] + '_' +  row['read_index'] + '.fastq.gz'), axis=1)
    # df_1kg_lr_rna = thousg_rna_long_read_samples()
    # df_1kg_lr_rna['hostname'] = df_1kg_lr_rna['url'].apply(lambda x: x.split('//')[1].split('/')[0])
    # df_1kg_lr_rna['ftp_path'] = df_1kg_lr_rna['url'].apply(lambda x:  '/' + '/'.join(x.split('//')[1].split('/')[1:]))
    if os.path.exists(VERIFIED_DOWNLOADS):
        df_verified = pd.read_csv(VERIFIED_DOWNLOADS, sep='\t')
    else:
        df_verified = pd.DataFrame(columns=['outfile', 'md5'])
    # only bother downloading files that are not already downloaded and verified
    df_1kg_sr_rna = df_1kg_sr_rna[~df_1kg_sr_rna['outfile'].isin(df_verified['outfile'])]
    args = []
    for i, row in df_1kg_sr_rna.iterrows():
        outfile = row['outfile']
        md5_given = row['md5']
        download_verified=False
        if os.path.exists(outfile):
            md5_calculated = get_checksum(outfile)
            download_verified = (md5_given == md5_calculated)
        if download_verified:
            print("# skipping {} from {} to {}: file already exists and checksum matches".format(row['ftp_path'], row['hostname'], outfile))
            continue
        else:
            args.append((row['hostname'], row['ftp_path'], outfile))
    # for i, row in df_1kg_lr_rna.iterrows():
    #     args.append((row['hostname'], row['ftp_path'], os.path.join(OUT_LR, os.path.basename(row['ftp_path'])), LOGFILE))
    # for a in args:
    #     print("# downloading {} from {} to {}".format(a[1], a[0], a[2]))
    #     downloader(*a)
    try:
        with mp.Pool(processes=CPUS) as pool:
            pool.starmap(downloader, args)
    except Exception as e:
        print("# error in multiprocessing: {}".format(e))

    # verify downloads
    for i, row in df_1kg_sr_rna.iterrows():
        outfile = row['outfile']
        md5_given = row['md5']
        if os.path.exists(outfile):
            md5_calculated = get_checksum(outfile)
            if md5_given == md5_calculated:
                print("# verified {}: checksum matches".format(outfile))
                df_verified = pd.concat([df_verified, pd.DataFrame({'outfile': [outfile], 'md5': [md5_given]})], ignore_index=True)
            else:
                print("# error verifying {}: checksum does not match (expected {}, got {})".format(outfile, md5_given, md5_calculated))
        else:
            print("# error verifying {}: file does not exist".format(outfile))
    
    df_verified.to_csv(VERIFIED_DOWNLOADS, sep='\t', index=False)

if __name__ == "__main__":
    main()
