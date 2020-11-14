#!/usr/bin/env python3
from libmkMeteo import mkmeteo4lingrars
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("netcdf", help="name of NETCDF ERA5 file to access")
parser.add_argument("rsdir", help="name of RS LAI/ET/etc directory")
parser.add_argument("becsmosdir", help="name of RS BEC SMOS 1km directory")
parser.add_argument("longitude", help="value of longitude")
parser.add_argument("latitude", help="value of latitude")
parser.add_argument("output", help="name of output text file")
args = parser.parse_args()
print(args.netcdf)
print(args.rsdir)
print(args.becsmosdir)
print(args.longitude)
print(args.latitude)
print(args.output)
#################################
meteolist = mkmeteo4lingrars(args.netcdf, args.rsdir, args.becsmosdir, args.longitude, args.latitude)
#################################
with open(args.output, 'w') as f:
    for item in meteolist:
        f.write(str(item))
        f.write("\n")
