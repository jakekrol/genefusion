from genefusion.genefusion import index_els
import os,sys

def main():
    args = sys.argv[1:]
    print(args)
    el, out = args
    index_els(el, out)
