#!/usr/bin/env python3
import datetime
import numpy as np
import xarray as xr
# Import local libraries
import libera5
import libgetgeo

parser.add_argument("netcdf", help="name of NETCDF file to access")
parser.add_argument("RSdir", help="name of RS LAI/ET/etc directory")
parser.add_argument("longitude", help="Longitude")
parser.add_argument("latitude", help="Latitude")
parser.add_argument("output", help="name of output Meteo txt file")

def libmkmeteo(netcdf, rsdir, longitude, latitude):
    """
    Prepare meteo and RS constraints for lingraRS model
    Requirements for local pylibs: libera5, libgetgeo
    Requirements for python libs: os, sys, netCDF4, numpy, datetime, argparse
    :param netcdf: input ERA5 map weather data
    :param rsdir: input directory where is RS data
    :param longitude: longitude (dd.ddd)
    :param latitude: latitude (dd.ddd)
    :return: list of lists of input parameters
    --------------------------------------------------------
    Column  Daily value
    1       year
    2       day
    3       irradiation                   (kJ m-2 d-1)
    4       minimum temperature           (degrees Celsius)
    5       maximum temperature           (degrees Celsius)
    6       early morning vapour pressure (kPa)
    7       mean wind speed (height: 10m) (m s-1)
    8       precipitation                 (mm d-1)
    9       RS evaporation                (mm d-1)
    10      RS transpiration              (mm d-1)
    11      RS LAI                        (m2.m-2)
    12      RS cut event                  (0/1)
    --------------------------------------------------------
    """
# Extract time array in python datetime
time_convert = libera5.getTimeNetcdf(args.netcdf)

# Set the basic array generation loop
layers = ['ssrd', 'u10', 'v10', 't2m', 'sp', 'tp', 'd2m']
for ll in range(len(layers)):
    array = libera5.getPixelVals(args.netcdf, layers[ll], args.longitude, args.latitude)
    if layers[ll] == 'ssrd':
        ssrd = np.copy(array)
    elif layers[ll] == 'u10':
        u10 = np.copy(array)
    elif layers[ll] == 'v10':
        v10 = np.copy(array)
    elif layers[ll] == 't2m':
        t2m = np.copy(array)
    elif layers[ll] == 'sp':
        sp = np.copy(array)
    elif layers[ll] == 'tp':
        tp = np.copy(array)
    elif layers[ll] == 'd2m':
        d2m = np.copy(array)
    else:
        print("name of this dataset is unknown")

# Make Wind speed (u10 x v10 components)
windspeed = np.abs(np.multiply(u10, v10))
# Solar Radiation J.m-2 -> kJ.m-2
ssrd = np.divide(ssrd, 1000)
for t in range(ssrd.shape[0]):
    if ssrd[t] < 0.001:
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
ssrd_agg = np.zeros_like(time_convert, dtype=np.float)
for i in range(0, ssrd.shape[0], 2):
    ssrd_agg[int(np.floor(i / 2))] = (ssrd[i] + ssrd[i + 1]) * 0.5

ssrd_agg[ssrd_agg == 0.0000] = np.nan
# [np.nan if x == 0.0000 else x for x in ssrd_agg]

# Convert from hour to day
dsas = xr.Dataset(
    {
        "ssrd": ("time", ssrd_agg),
    },
    {"time": time_convert},
)
# Make mean on non NAN values
# dsasD=dsas.resample(time='1D').mean(skipna=True)
# Make mean on non NAN values
dsasD = dsas.resample(time='1D').reduce(np.nansum)

wspd_agg = np.zeros_like(time_convert, dtype=np.float)
for i in range(0, windspeed.shape[0], 2):
    wspd_agg[int(np.floor(i / 2))] = (windspeed[i] + windspeed[i + 1]) * 0.5

dsa = xr.Dataset(
    {
        "wspd": ("time", wspd_agg),
    },
    {"time": time_convert},
)
dsaD = dsa.resample(time='1D').mean()

# Precipitation summed by day
pmm_agg = np.zeros_like(time_convert, dtype=np.float)
for i in range(0, tp.shape[0], 2):
    pmm_agg[int(np.floor(i / 2))] = (tp[i] + tp[i + 1]) * 0.5

dsap = xr.Dataset(
    {
        "prmm": ("time", pmm_agg),
    },
    {"time": time_convert},
)
dsapD = dsap.resample(time='1D').sum()

# Minimum temperature per day
t2m_agg = np.zeros_like(time_convert, dtype=np.float)
for i in range(0, t2m.shape[0], 2):
    t2m_agg[int(np.floor(i / 2))] = (t2m[i] + t2m[i + 1]) * 0.5

dstmin = xr.Dataset(
    {
        "tmin": ("time", t2m_agg),
    },
    {"time": time_convert},
)
dstminD = dstmin.resample(time='1D').min()

# Maximum temperature per day
dstmax = xr.Dataset(
    {
        "tmax": ("time", t2m_agg),
    },
    {"time": time_convert},
)
dstmaxD = dstmax.resample(time='1D').max()

# Minimum dew point temperature per day
d2m_agg = np.zeros_like(time_convert, dtype=np.float)
for i in range(0, d2m.shape[0], 2):
    d2m_agg[int(np.floor(i / 2))] = (d2m[i] + d2m[i + 1]) * 0.5

dsdmin = xr.Dataset(
    {
        "dmin": ("time", d2m_agg),
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
             np.array(dstmaxD.tmax).tolist(), eact, np.array(dsaD.wspd).tolist(), np.array(dsapD.prmm).tolist()]
# INSERT RS DATA HERE
####################
# Create E & T from RS
# Get the FC files in numpy array converted to [0.0-1.0]
FC = libgetgeo.getgeoydoy(args.RSdir, "MOD13Q1*_FC.tif", args.longitude, args.latitude, year.tolist(), doy, i=True)
FC = np.asarray(FC, dtype=np.float)
FC = np.divide(FC, 100)
# print("START FC", FC, "END FC")
# Get ETa from both MOD and MYD platforms in arrays
MODET = libgetgeo.getgeoydoy(args.RSdir, "MOD16A*.vrt", args.longitude, args.latitude, year.tolist(), doy)
MYDET = libgetgeo.getgeoydoy(args.RSdir, "MYD16A*.vrt", args.longitude, args.latitude, year.tolist(), doy)
# TODO make 32765 -> NAN in MODET/MYDET
MODET = np.asarray(MODET, dtype=np.float)
MYDET = np.asarray(MYDET, dtype=np.float)
MODET[MODET == 32765] = np.nan
MYDET[MYDET == 32765] = np.nan
MODET = np.divide(MODET, 10)
MYDET = np.divide(MYDET, 10)
# print("START MODET", MODET, "END MODET")
# print("START MYDET", MYDET, "END MYDET")
# Compute average ETa. Take NANs into account to maximise data count out
ET = (MODET + MYDET) / 2.0
# print("START ET", ET, "END ET")
# Transpiration
T = np.multiply(ET, FC)
# print("START T", T, "END T")
# Evaporation is the complementary of Transpiration vis-a-vis ET
E = np.subtract(ET, T)
# print("START E", E, "END E")
# TODO check the output format of the arrays to be compatible to the csv format
# Add Evaporation data to the Meteo file
meteolist.append(E.tolist())
# Add Transpiration data to the Meteo file
meteolist.append(T.tolist())
# Add LAI data to the Meteo file
LAI = libgetgeo.getgeoydoy(args.RSdir, "MCD15A2H*.vrt", args.longitude, args.latitude, year.tolist(), doy)
LAI = np.array(LAI, dtype=np.float)
LAI = np.divide(LAI, 10)
# print("START LAI", LAI, "END LAI")
meteolist.append(LAI.tolist())
# Create artificial Cut data from a previous array
CutNans = np.zeros(LAI.shape[0])
# Fill with NANs
CutNans.fill(np.nan)
# Fill with 2 cuts
CutNans[127] = 1
CutNans[235] = 1
# Add Cut data to the Meteo file
meteolist.append(CutNans.tolist())
# print("START CUTNANS", CutNans, "END CUTNANS")

#################################
with open(args.output, 'w') as f:
    for item in meteolist:
        f.write(str(item))
        f.write("\n")
