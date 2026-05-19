#!/usr/bin/env python
from polymerization.star import fusiontsv2bedpe
import argparse

parser = argparse.ArgumentParser(description="Convert STAR-Fusion output TSV to BEDPE format")
parser.add_argument("--input", help="Path to STAR-Fusion output")
parser.add_argument("--output", help="Output BEDPE file path")
parser.add_argument("--bgzip", action="store_true", help="Compress output with bgzip")
parser.add_argument("--pad", action="store_true", help="Pad output with zeros for compatibility with downstream tools")
args = parser.parse_args()

fusiontsv2bedpe(args.input, args.output, bgzip=args.bgzip, pad=args.pad)
