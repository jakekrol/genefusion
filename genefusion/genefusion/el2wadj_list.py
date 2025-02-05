#!/usr/bin/env python3

from genefusion.genefusion import el2wadj_list
import os,sys

def main():
    args = sys.argv[1:]
    print(args)
    el, out = args
    el2wadj_list(el, out)
