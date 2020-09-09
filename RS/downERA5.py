#!/usr/bin/env python3
from libera5 import *
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("north", help="north")
parser.add_argument("south", help="south")
parser.add_argument("east", help="east")
parser.add_argument("west", help="west")
parser.add_argument("year", help="year (YYYY)")
parser.add_argument("outdir", help="output directory")
parser.add_argument("prefix", help="name of this specific download")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args=parser.parse_args()

#Feed the beast from args parser
north, south = args.north, args.south
east, west = args.east, args.west
outdir, prefix = args.outdir, args.prefix
i_year = args.year

download_ERA5_bbox(i_year, outdir, prefix, north, west, south, east)

