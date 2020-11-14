#!/usr/bin/env python3
import argparse
import sys
from time import sleep

from osgeo import gdal
# Parallelization of runs, has to be defined before local libs!
import ray
# Import local libraries
import libera5
from libprocessLingraRow import processlingrarow

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
parser.add_argument("RSdir", help="name of RS LAI/ET/etc directory")
parser.add_argument("SMdir", help="name of RS BEC-SMOS 1Km directory")
if len(sys.argv) < 5:
    parser.print_help(sys.stderr)
    print(Requirements)
    print(Header)
    sys.exit(1)
args = parser.parse_args()

# Extract time array in python datetime
time_convert = libera5.getTimeNetcdf(args.netcdf)

# -------------------------------------------------
# Grassland Copernicus pixels location extraction
# -------------------------------------------------
driver = gdal.GetDriverByName('GTiff')
filename = args.grassland  # path to raster
dataset = gdal.Open(filename)
band = dataset.GetRasterBand(1)

cols = dataset.RasterXSize
rows = dataset.RasterYSize
transform = dataset.GetGeoTransform()
xOrigin = transform[0]
yOrigin = transform[3]
pixelWidth = transform[1]
pixelHeight = -transform[5]

data = band.ReadAsArray(0, 0, cols, rows)
# Create output files
# Set no_data value
no_data = 7777777
# Set pyramid building option to Erdas
gdal.SetConfigOption('HFA_USE_RRD', 'YES')
# Creation options
opts = ['TILED=YES', 'COMPRESS=DEFLATE', 'NUM_THREADS=ALL_CPUS']
# Build output files
# tiller, yielD, wlvg, wlvd1, parcu, grass, tracu, evacu
out0 = driver.CreateCopy("0_tiller.tif", dataset, 0, opts)
out1 = driver.CreateCopy("1_yield.tif", dataset, 0, opts)
out2 = driver.CreateCopy("2_wlvg.tif", dataset, 0, opts)
out3 = driver.CreateCopy("3_wlvd1.tif", dataset, 0, opts)
out4 = driver.CreateCopy("4_parcu.tif", dataset, 0, opts)
out5 = driver.CreateCopy("5_grass.tif", dataset, 0, opts)
out6 = driver.CreateCopy("6_tracu.tif", dataset, 0, opts)
out7 = driver.CreateCopy("7_evacu.tif", dataset, 0, opts)
# Access band handles for output files
b0 = out0.GetRasterBand(1)
b1 = out1.GetRasterBand(1)
b2 = out2.GetRasterBand(1)
b3 = out3.GetRasterBand(1)
b4 = out4.GetRasterBand(1)
b5 = out5.GetRasterBand(1)
b6 = out6.GetRasterBand(1)
b7 = out7.GetRasterBand(1)
# Set no_data value
b0.SetNoDataValue(no_data)
b1.SetNoDataValue(no_data)
b2.SetNoDataValue(no_data)
b3.SetNoDataValue(no_data)
b4.SetNoDataValue(no_data)
b5.SetNoDataValue(no_data)
b6.SetNoDataValue(no_data)
b7.SetNoDataValue(no_data)
# Access data from band handles as numpy arrays
d0 = b0.ReadAsArray(0, 0, cols, rows)
d1 = b1.ReadAsArray(0, 0, cols, rows)
d2 = b2.ReadAsArray(0, 0, cols, rows)
d3 = b3.ReadAsArray(0, 0, cols, rows)
d4 = b4.ReadAsArray(0, 0, cols, rows)
d5 = b5.ReadAsArray(0, 0, cols, rows)
d6 = b6.ReadAsArray(0, 0, cols, rows)
d7 = b7.ReadAsArray(0, 0, cols, rows)
# Fill with nan, so that only grass pixels get to be taken care of.
d0.fill(no_data)
d1.fill(no_data)
d2.fill(no_data)
d3.fill(no_data)
d4.fill(no_data)
d5.fill(no_data)
d6.fill(no_data)
d7.fill(no_data)

# DEBUG ONLY: True
plot = False

# Parallel code
ray.init(local_mode=True, address="localhost:6379")
#, num_cpus=16, ignore_reinit_error=True)
#ray.init(num_cpus=7, num_gpus=1) #, memory= < bytes >, object_store_memory = < bytes >)

# Start processing the model for each grass pixel
for c in range(cols):
    # for r in range(rows):
    # (c, r, d0[c][r], d1[c][r], d2[c][r], d3[c][r], d4[c][r], d5[c][r], d6[c][r], d7[c][r]) = \
    # processlingrapixel(r, data[c][r], pixelWidth, pixelHeight, xOrigin, yOrigin, plot, args.netcdf, args.RSdir)
    lingrarow = [processlingrarow.remote(c, rows, data[c], pixelWidth, pixelHeight, xOrigin, yOrigin,
                                           plot, args.netcdf, args.RSdir, args.SMdir) for r in range(rows)]
    for row in range(rows):
        ready_ids, not_ready_ids = ray.wait(lingrarow[row])
        while len(not_ready_ids) > 1:
            sleep(10)
            ready_ids, not_ready_ids = ray.wait(lingrarow[row])
            print(not_ready_ids)

    for row in range(row):
        output = ray.get(lingrarow[row])
        row = output[0]
        print(output[0])
        print(output[1])
        print(output[2])
        d0[c][row] = output[1]
        d1[c][row] = output[2]
        d2[c][row] = output[3]
        d3[c][row] = output[4]
        d4[c][row] = output[5]
        d5[c][row] = output[6]
        d6[c][row] = output[7]
        d7[c][row] = output[8]

# Write arrays to files
b0.WriteArray(d0)
b1.WriteArray(d1)
b2.WriteArray(d2)
b3.WriteArray(d3)
b4.WriteArray(d4)
b5.WriteArray(d5)
b6.WriteArray(d6)
b7.WriteArray(d7)

# Build overviews for all new files
out0.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out1.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out2.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out3.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out4.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out5.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out6.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])
out7.BuildOverviews(overviewlist=[2, 4, 8, 16, 32, 64, 128])

# Commit the changes to the files
out0.FlushCache()
out1.FlushCache()
out2.FlushCache()
out3.FlushCache()
out4.FlushCache()
out5.FlushCache()
out6.FlushCache()
out7.FlushCache()

# Close properly the dataset
out0 = None
out1 = None
out2 = None
out3 = None
out4 = None
out5 = None
out6 = None
out7 = None

# Complete the process by signalling the user
print('\033[32m', "Processing DONE", '\033[0m', sep='')
