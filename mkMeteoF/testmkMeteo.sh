#!/bin/bash
#For gdal tools
export PATH=$PATH:/Library/Frameworks/GDAL.framework/Versions/3.1/Programs/

#ERA5 file access handle
netcdf="$HOME/Documents/ERA5/ERA5_EU_2020.nc"	#name of NETCDF file to access
#RS directory access handle
RSdir="$HOME/Documents/MODIS"	#name of RS LAI/ET/etc directory
#BEC SMOS 1km data dir
becsmosdir="$HOME/Documents/SMOS"

#TEST long lat
longitude="0.0"
latitude="47.0"
output="mkMeteoF/rs.csv" 	#name of output Meteo txt file
./mkMeteoF/mkMeteo.py "$netcdf" "$RSdir" "$becsmosdir" $longitude $latitude $output

