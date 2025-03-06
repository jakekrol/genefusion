from genefusion.genefusion import el2json
import os,sys

def main():
    args = sys.argv[1:]
    print(args)
    el, out = args
    el2json(el, out)
