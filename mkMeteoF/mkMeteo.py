#!/usr/bin/env python3
import libmkMeteo
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("netcdf", help="name of NETCDF ERA5 file to access")
parser.add_argument("rsdir", help="name of RS LAI/ET/etc directory")
parser.add_argument("becsmosdir", help="name of RS BEC SMOS 1km directory")
parser.add_argument("longitude", help="value of longitude")
parser.add_argument("latitude", help="value of latitude")
parser.add_argument("output", help="name of output text file")
parser.add_argument("pickle", help="output to pickle (True/False)")
args = parser.parse_args()
print(args.netcdf)
print(args.rsdir)
print(args.becsmosdir)
print(args.longitude)
print(args.latitude)
print(args.output)
print(args.pickle)
#################################
meteolist = libmkMeteo.mkmeteo4lingrars(args.netcdf, args.rsdir, args.becsmosdir, args.longitude, args.latitude)
#################################
if not args.pickle:
    with open(args.output, 'w') as f:
        for item in meteolist:
            f.write(str(item))
            f.write("\n")
    f.close()
else:
    with open(args.output, 'wb') as f:
        pickle.dump(meteolist, f)
    f.close()