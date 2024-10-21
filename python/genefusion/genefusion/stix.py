#!/usr/bin/env python
from genefusion.genefusion import stix
import os,sys
# example
# genefusion-stix index stix.ped.db 21 39751949 40033704 42836478 42903043 erg.tmprss2.stix.out DEL

# objs assigned here aren't used. for reference. 
# use unpack operator instead for simpler code
def main():
    args = sys.argv[1:]
    giggle_idx = args[0]
    stix_db = args[1]
    chr = args[2]
    left_start = args[3]
    left_end = args[4]
    right_start = args[5]
    right_end = args[6]
    outfile = args[7]
    type = args[8]
    #gene_file =args[10]
    #strand = args[11]

    file,*stix_args = sys.argv
    print(stix_args)
    stix(*stix_args)



