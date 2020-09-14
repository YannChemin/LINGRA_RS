#!/usr/bin/env python3
import argparse
import datetime
import sys

import numpy as np

# Import local libraries
import libera5
import libgetgeo
from osgeo import gdal

# import main lingraRS library
from liblingraRS import lingrars

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
*---------------------------------------------------------------------------*
"""

parser = argparse.ArgumentParser()
parser.add_argument("grassland", help="name of Copernicus Grassland file to access (INT32!)")
parser.add_argument("netcdf", help="name of NETCDF ERA5 file to access")
parser.add_argument("RSdir", help="name of RS LAI/ET/etc directory")
if len(sys.argv) < 4:
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
pixelWidth = int(transform[1])
pixelHeight = int(-transform[5])

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

# Start processing the model for each grass pixel
for col in range(cols):
    for row in range(rows):
        if data[col][row] == 1:
            longitude = col * pixelWidth + xOrigin
            latitude = yOrigin - row * pixelHeight
            print(longitude, latitude, data[col][row])
            # Set the basic array generation loop
            layers = ['ssrd', 'u10', 'v10', 't2m', 'sp', 'tp', 'd2m']
            for lay in range(len(layers)):
                array = libera5.getPixelVals(args.netcdf, layers[lay], longitude, latitude)
                if layers[lay] == 'ssrd':
                    ssrd = np.copy(array)
                elif layers[lay] == 'u10':
                    u10 = np.copy(array)
                elif layers[lay] == 'v10':
                    v10 = np.copy(array)
                elif layers[lay] == 't2m':
                    t2m = np.copy(array)
                elif layers[lay] == 'sp':
                    sp = np.copy(array)
                elif layers[lay] == 'tp':
                    tp = np.copy(array)
                elif layers[lay] == 'd2m':
                    d2m = np.copy(array)
                else:
                    print("name of this dataset is unknown")

            # Make Wind speed (u10 x v10 components)
            windspeed = np.abs(np.multiply(u10, v10))
            # Solar Radiation J.m-2 -> kJ.m-2
            ssrd = np.divide(ssrd, 1000)
            for t in range(ssrd.shape[0]):
                if ssrd[t] < 1.0:
                    ssrd[t] = np.nan
            # Surface Pressure Pa -> hPa
            sp = np.divide(sp, 100)
            # Temperature K -> C
            t2m = np.subtract(t2m, 273.15)
            # Precipitation m -> mm
            tp = np.abs(np.multiply(tp, 1000))
            # Temperature K -> C
            d2m = np.subtract(d2m, 273.15)

            # Before daily conversion, remove night values for some
            # Replace nodata with NAN
            [np.nan if x == 0.0000 else x for x in ssrd]

            # Convert from hour to day
            import xarray as xr

            dsas = xr.Dataset(
                {
                    "ssrd": ("time", ssrd),
                },
                {"time": time_convert},
            )
            # Make mean on non NAN values
            # dsasD = dsas.resample(time='1D').mean(skipna=True)
            # Make mean on non NAN values
            dsasD = dsas.resample(time='1D').sum(skipna=True)

            dsa = xr.Dataset(
                {
                    "wspd": ("time", windspeed),
                    "srfp": ("time", sp),
                },
                {"time": time_convert},
            )
            dsaD = dsa.resample(time='1D').mean()

            # Precipitation summed by day
            dsap = xr.Dataset(
                {
                    "prmm": ("time", tp),
                },
                {"time": time_convert},
            )
            dsapD = dsap.resample(time='1D').sum()

            # Minimum temperature per day
            dstmin = xr.Dataset(
                {
                    "tmin": ("time", t2m),
                },
                {"time": time_convert},
            )
            dstminD = dstmin.resample(time='1D').min()

            # Maximum temperature per day
            dstmax = xr.Dataset(
                {
                    "tmax": ("time", t2m),
                },
                {"time": time_convert},
            )
            dstmaxD = dstmax.resample(time='1D').max()

            # Minimum dew point temperature per day
            dsdmin = xr.Dataset(
                {
                    "dmin": ("time", d2m),
                },
                {"time": time_convert},
            )
            dsdminD = dsdmin.resample(time='1D').min()

            # Convert Dewpoint T Min (C) to Eact (hPa)
            eactt = np.array(dsdminD.dmin)
            eact = []
            for d in range(len(eactt)):
                eact.append(libera5.d2m2eact(eactt[d]))

            # Generate year value for every day
            year = np.datetime_as_string(dsaD.time, unit='Y')
            # Generate DOY value for every day
            doyt = np.datetime_as_string(dsaD.time, unit='D')
            doy = []
            for t in range(len(doyt)):
                doy.append(datetime.datetime.strptime(doyt[t], '%Y-%m-%d').timetuple().tm_yday)

            # Join all data
            meteolist = [year.tolist(), doy, np.array(dsasD.ssrd).tolist(), np.array(dstminD.tmin).tolist(),
                         np.array(dstmaxD.tmax).tolist(), eact, np.array(dsaD.wspd).tolist(),
                         np.array(dsapD.prmm).tolist()]
            # meteolist.append(np.array(dsaD.srfp).tolist())
            # INSERT RS DATA HERE
            ####################
            # Create E & T from RS
            # Get the FV files in array
            FV = libgetgeo.getgeoydoy(args.RSdir, "MOD13Q1*_FV.tif", longitude, latitude, year.tolist(), doy, i=True)
            print("START FV", FV, "END FV")
            # Get ETa from both MOD and MYD platforms in arrays
            MODET = libgetgeo.getgeoydoy(args.RSdir, "MOD16A*.vrt", longitude, latitude, year.tolist(), doy)
            MYDET = libgetgeo.getgeoydoy(args.RSdir, "MYD16A*.vrt", longitude, latitude, year.tolist(), doy)
            print("START MODET", MODET, "END MODET")
            print("START MYDET", MYDET, "END MYDET")
            # Compute average ETa. Take NANs into account to maximise data count out
            ET = np.nanmean(MODET, MYDET)
            print("START ET", ET, "END ET")
            # Transpiration
            T = np.multiply(ET, FV)
            print("START T", T, "END T")
            # Evaporation is the complementary of Transpiration vis-a-vis ET
            E = np.subtract(ET, T)
            print("START E", E, "END E")
            # Add Evaporation data to the Meteo file
            meteolist.append(E)
            # Add Transpiration data to the Meteo file
            meteolist.append(T)
            # Add LAI data to the Meteo file
            LAI = libgetgeo.getgeoydoy(args.RSdir, "MCD15A2H*.vrt", longitude, latitude, year.tolist(), doy)
            print("START LAI", LAI, "END LAI")
            meteolist.append(libgetgeo.getgeoydoy(args.RSdir, "MCD15A2H*.vrt",
                                                  longitude, latitude, year.tolist(), doy))
            # Create artificial Cut data from a previous array
            Cut = np.zeros(np.array(dsapD.prmm))
            # Fill with NANs
            CutNans = [np.nan if x == 0 else x for x in Cut]
            # Fill with 2 cuts
            CutNans[127] = 1
            CutNans[235] = 1
            # Add Cut data to the Meteo file
            meteolist.append(CutNans)
            print("START LAI", LAI, "END LAI")
            # Launch LINGRA_RS and retrieve the yield into a pixel
            plot = False
            (tiller, yielD, wlvg, wlvd1, parcu, grass, tracu, evacu) = lingrars(latitude, meteolist, plot)
            # Let the pixels fit into each map (*10000 bc INT32 maps)
            d0[col][row] = tiller * 10000
            d1[col][row] = yielD * 10000
            d2[col][row] = wlvg * 10000
            d3[col][row] = wlvd1 * 10000
            d4[col][row] = parcu * 10000
            d5[col][row] = grass * 10000
            d6[col][row] = tracu * 10000
            d7[col][row] = evacu * 10000

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
