from genefusion.genefusion import g2ewdist
import os,sys

def main():
    args = sys.argv[1:]
    print(args)
    adj, out = args
    g2ewdist(adj, out)
