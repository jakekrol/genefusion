#!/usr/bin/env python3

import random
import argparse

parser = argparse.ArgumentParser(description='Tissue-wise score hists')
parser.add_argument('--score', required=True)
parser.add_argument('--seed', type=int, default=0)
parser.add_argument('--sample_size', type=int, default=100000)
parser.add_argument('--out', required=True)
