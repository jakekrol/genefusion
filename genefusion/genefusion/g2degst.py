from genefusion.genefusion import g2degst
import os,sys

def main():
    args = sys.argv[1:]
    print(args)
    adj, out = args
    g2degst(adj, out)
