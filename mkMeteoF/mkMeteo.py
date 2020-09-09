#!/usr/bin/env python

Requirements=
"""
*---------------------------------------------------------------------------*
* Requirements for Linux OS: datamash, paste, sed
* Requirements for local pylibs: libera5, libgetgeo
* Requirements for python libs: os, sys, netCDF4, numpy, datetime, argparse
*---------------------------------------------------------------------------*
"""

Header=
"""
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

import os, sys
import netCDF4
import numpy as np
import datetime
#Import local libraries
import libera5
import libgetgeo

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("netcdf", help="name of NETCDF file to access")
parser.add_argument("RSdir", help="name of RS LAI/ET/etc directory")
parser.add_argument("longitude", help="Longitude")
parser.add_argument("latitude", help="Latitude")
parser.add_argument("output", help="name of output Meteo txt file")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    print(Requirements)
    print(Header)
    sys.exit(1)
args = parser.parse_args()

#Extract time array in python datetime
time_convert=libera5.getTimeNetcdf(args.netcdf)

#Set the basic array generation loop
layers=['ssrd','u10','v10','t2m','sp','tp','d2m']
for l in range(len(layers)):
    array=libera5.getPixelVals(args.netcdf, layers[l], args.longitude, args.latitude)
    if(layers[l]=='ssrd'):
        ssrd=np.copy(array)
    elif(layers[l]=='u10'):
        u10=np.copy(array)
    elif(layers[l]=='v10'):
        v10=np.copy(array)
    elif(layers[l]=='t2m'):
        t2m=np.copy(array)
    elif(layers[l]=='sp'):
        sp=np.copy(array)
    elif(layers[l]=='tp'):
        tp=np.copy(array)
    elif(layers[l]=='d2m'):
        d2m=np.copy(array)
    else:
        print("name of this dataset is unknown")


#Make Wind speed (u10 x v10 components) 
windspeed=np.abs(np.multiply(u10,v10))
#Solar Radiation J.m-2 -> kJ.m-2
ssrd=np.divide(ssrd,1000)
for t in range(ssrd.shape[0]):
    if(ssrd[t]<1.0):
        ssrd[t]=np.nan
#Surface Pressure Pa -> hPa
sp=np.divide(sp,100)
#Temperature K -> C
t2m=np.subtract(t2m,273.15)
#Precipitation m -> mm
tp=np.abs(np.multiply(tp,1000))
#Temperature K -> C
d2m=np.subtract(d2m,273.15)

#Before daily conversion, remove night values for some
#Replace nodata with NAN
[np.nan if x==0.0000 in x else x for x in ssrd]

#Convert from hour to day
import xarray as xr
dsas = xr.Dataset(
    {
        "ssrd": (("time"), ssrd),
    },
    {"time": time_convert},
)
#Make mean on non NAN values
#dsasD=dsas.resample(time='1D').mean(skipna=True)
#Make mean on non NAN values
dsasD=dsas.resample(time='1D').sum(skipna=True)

dsa = xr.Dataset(
    {
        "wspd": (("time"), windspeed),
        "srfp": (("time"), sp),
    },
    {"time": time_convert},
)
dsaD=dsa.resample(time='1D').mean()

#Precipitation summed by day
dsap = xr.Dataset(
    {
        "prmm": (("time"), tp),
    },
    {"time": time_convert},
)
dsapD=dsap.resample(time='1D').sum()

#Minimum temperature per day
dstmin = xr.Dataset(
    {
        "tmin": (("time"), t2m),
    },
    {"time": time_convert},
)
dstminD=dstmin.resample(time='1D').min()

#Maximum temperature per day
dstmax = xr.Dataset(
    {
        "tmax": (("time"), t2m),
    },
    {"time": time_convert},
)
dstmaxD=dstmax.resample(time='1D').max()

#Minimum dew point temperature per day
dsdmin = xr.Dataset(
    {
        "dmin": (("time"), d2m),
    },
    {"time": time_convert},
)
dsdminD=dsdmin.resample(time='1D').min()

#Convert Dewpoint T Min (C) to Eact (hPa) 
eactt=np.array(dsdminD.dmin)
eact=[]
for d in range(len(eactt)):
    eact.append(libera5.d2m2eact(eactt[d]))

import pandas as pd
#Generate year value for every day
year=np.datetime_as_string(dsaD.time,unit='Y')
#Generate DOY value for every day
doyt=np.datetime_as_string(dsaD.time,unit='D')
doy=[]
for t in range(len(doyt)):
    doy.append(datetime.datetime.strptime(doyt[t],'%Y-%m-%d').timetuple().tm_yday)

#Join all data
listA=[]
listA.append(year.tolist())
listA.append(doy)
listA.append(np.array(dsasD.ssrd).tolist())
listA.append(np.array(dstminD.tmin).tolist())
listA.append(np.array(dstmaxD.tmax).tolist())
listA.append(eact)
#listA.append(np.array(dsaD.srfp).tolist())
listA.append(np.array(dsaD.wspd).tolist())
listA.append(np.array(dsapD.prmm).tolist())
#INSERT RS DATA HERE
####################
#Create E & T from RS
#Get the FV files in array 
FV=getgeoydoy(args.RSdir, "MOD13Q1*_FV.tif", args.longitude, args.latitude, year.tolist(), doy)
print("START FV",FV,"END FV")
#TODO expand the FV data to all days, include filling NANs
#Get ETa from both MOD and MYD platforms in arrays
MODET=getgeoydoy(args.RSdir, "MOD16A*.vrt", args.longitude, args.latitude, year.tolist(), doy)
MYDET=getgeoydoy(args.RSdir, "MYD16A*.vrt", args.longitude, args.latitude, year.tolist(), doy)
print("START MODET",MODET,"END MODET")
print("START MYDET",MYDET,"END MYDET")
#Compute average ETa. Take NANs into account to maximise data count out
ET=np.nanmean(MODET,MYDET)
print("START ET",ET,"END ET")
#Transpiration TODO Check how NANs are handled, maybe FV filling is easier
T=np.multiply(ET,FV)
print("START T",T,"END T")
#Evaporation is the complementary of Transpiration vis-a-vis ET
E=np.subtract(ET,T)
print("START E",E,"END E")
#TODO check the output format of the arrays to be compatible to the csv format
# Add Evaporation data to the Meteo file
listA.append(E)
# Add Transpiration data to the Meteo file
listA.append(T)
# Add LAI data to the Meteo file
LAI=getgeoydoy(args.RSdir, "MCD15A2H*.vrt", args.longitude, args.latitude, year.tolist(), doy)
print("START LAI",LAI,"END LAI")
listA.append(getgeoydoy(args.RSdir, "MCD15A2H*.vrt", args.longitude, args.latitude, year.tolist(), doy))
# Create artificial Cut data from a previous array
Cut=np.zeros(np.array(dsapD.prmm)) 
# TODO Fill with NANs
CutNans=CutNans[ np.nan for x in Cut ]
#Fill with 2 cuts
CutNans[127]=1
CutNans[235]=1
# Add Cut data to the Meteo file
listA.append(CutNans)
print("START LAI",LAI,"END LAI")

#################################
with open(args.output, 'w') as f:
    for item in listA:
        f.write(f"{item}\n")

########################################################
# Housekeeping for csv format
os.system('sed -i "s/\'//g" %s' % (args.output))
os.system('sed -i "s/,\\ /\t/g" %s' % (args.output))
os.system('sed -i "s/\[//g" %s' % (args.output))
os.system('sed -i "s/\]//g" %s' % (args.output))
os.system('datamash transpose < %s > t%s' % (args.output, args.output))
os.system('paste t%s rs.csv > %s' % (args.output, args.output))
#########################################################
