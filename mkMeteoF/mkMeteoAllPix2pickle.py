#!/usr/bin/env python3
import argparse
import sys, os
from subprocess import check_output
import multiprocessing
from time import sleep
from osgeo import gdal

# Import local libraries
import libera5
import libmkMeteo

def get_pid(name):
    return list(map(int,check_output(["pidof",name]).split()))

Requirements = """
*---------------------------------------------------------------------------*
* Requirements for Linux OS: datamash, paste, sed
* Requirements for local pylibs: libera5, libgetgeo
* Requirements for python libs: os, sys, netCDF4, numpy, datetime, argparse, gdal
*---------------------------------------------------------------------------*
"""

Header = """
*---------------------------------------------------------------------------*
* Station name: ERA5 Download
* Author: Yann Chemin
* Source: ECMWF CDSAPI Download
*
* Column  Daily value
* 1       year
* 2       day
* 3       irradiation                   (kJ m-2 d-1)
* 4       minimum temperature           (degrees Celsius)
* 5       maximum temperature           (degrees Celsius)
* 6       early morning vapour pressure (kPa)
* 7       mean wind speed (height: 10m) (m s-1)
* 8       precipitation                 (mm d-1)
* 9       RS evaporation                (mm d-1)
* 10      RS transpiration              (mm d-1)
* 11      RS LAI                        (m2.m-2)
* 12      RS cut event                  (0/1)
* 13      RS soil moisture              (cm3/cm3)
*---------------------------------------------------------------------------*
"""

parser = argparse.ArgumentParser()
parser.add_argument("grassland", help="name of Copernicus Grassland file to access (INT32!)")
parser.add_argument("netcdf", help="name of NETCDF ERA5 file to access")
parser.add_argument("rsdir", help="name of RS LAI/ET/etc directory")
parser.add_argument("smdir", help="name of RS BEC-SMOS 1Km directory")
parser.add_argument("output", help="name of output filename base for pickled lists")
if len(sys.argv) < 6:
    parser.print_help(sys.stderr)
    print(Requirements)
    print(Header)
    sys.exit(1)
args = parser.parse_args()

homedirout = "/Users/dnd/Downloads/data_"
# Extract time array in python datetime
time_convert = libera5.getTimeNetcdf(args.netcdf)

# -------------------------------------------------
# Grassland Copernicus pixels location extraction
# -------------------------------------------------
driver = gdal.GetDriverByName('GTiff')
filename = args.grassland  # path to raster
dataset = gdal.Open(filename)
transform = dataset.GetGeoTransform()
xO = transform[0]
yO = transform[3]
pixelWidth = transform[1]
pixelHeight = -transform[5]
band = dataset.GetRasterBand(1)
cols = dataset.RasterXSize
rows = dataset.RasterYSize
data = band.ReadAsArray(0, 0, cols, rows)

#Try to get the number of cores available
try:
    ncores = multiprocessing.cpu_count()
except:
    ncores = 8

# TODO Temporary force not to overwhelm poor lil'Mac
ncores=4

# Define slave name to search in ps aux
slavename = "mkMeteo.py"

# Start processing the model for each grass pixel
for c in range(cols):
    print("data[%d/%d]" % (c, cols))
    for r in range(rows):
        # Hold on more processing as long as we have too many slave PIDs
        while len(get_pid(slavename)) >= ncores:
            sleep(10)
        # Process the pixel
        if data[c][r] == 1:
            filename = homedirout + str(c) + "_" + str(r) + "_" + args.output
            try:
                f = open(filename)
                f.close()
            except FileNotFoundError:
                print("data[%d/%d][%d/%d]" % (c, cols, r, rows))
                longitude = c * pixelWidth + xO
                latitude = yO - r * pixelHeight
                # print(col, row, longitude, latitude, data)
                pickle = True
                os.system("/bin/bash /Users/dnd/Documents/GitHub/LINGRA_RS/mkMeteoF/mkMeteo.sh %s %s %s %f %f %s %s"
                          % (args.netcdf, args.rsdir, args.smdir, longitude, latitude, args.output, pickle))
                # meteolist = libmkMeteo.mkmeteo4lingrars(netcdf, rsdir, smdir, longitude, latitude)
                # with open(filename, 'wb') as f:
                #    pickle.dump(meteolist, f)
                # f.close()

# Complete the process by signalling the user
print('\033[32m', "Processing DONE", '\033[0m', sep='')
