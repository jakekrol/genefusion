#!/usr/bin/env python3
from jkbiolib.datasets.loaders import *
from ftplib import FTP
from urllib.parse import urlparse

df = thousg_rna_short_read_samples()
fastq_urls = df.iloc[0:2,:]['url'].values
# ftp download
ftp = FTP('ftp.sra.ebi.ac.uk')
ftp.login()
for url in fastq_urls:
    print(f'Downloading {url}...')
    parsed = urlparse(url)
    remote_path = parsed.path  # e.g., /vol1/fastq/SRR197/...
    filename = parsed.path.split('/')[-1]
    with open(filename, 'wb') as f:
        ftp.retrbinary('RETR ' + remote_path, f.write)
ftp.quit()