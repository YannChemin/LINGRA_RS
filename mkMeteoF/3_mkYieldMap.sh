#!/bin/bash
#For gdal tools
export PATH=$PATH:/Library/Frameworks/GDAL.framework/Versions/3.1/Programs/

####################
#parser.add_argument("grassland", help="name of Copernicus Grassland file to access")
#parser.add_argument("netcdf", help="name of NETCDF file to access")
#parser.add_argument("RSdir", help="name of RS LAI/ET/etc directory")
####################

grassland="$HOME/Documents/GRA_2018_10m/GRA_2018_BE_WA_BRABANT.tif"
netcdf="$HOME/Documents/ERA5/ERA5_EU_2020.nc"
RSdir="$HOME/Documents/MODIS/"
smosdir="$HOME/Documents/SMOS/"


ray start --head --port=6379
python3 ./mkMeteoAllPix.py "$grassland" "$netcdf" "$RSdir" "$smosdir"
ray stop