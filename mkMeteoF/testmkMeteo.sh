#!/bin/bash
#For gdal tools
export PATH=$PATH:/Library/Frameworks/GDAL.framework/Versions/3.1/Programs/

#ERA5 file access handle
netcdf="$HOME/Documents/ERA5/ERA5_EU_2020.nc"	#name of NETCDF file to access
#RS directory access handle
RSdir="$HOME/Documents/MODIS"	#name of RS LAI/ET/etc directory

longitude=-10.0
latitude=60.0
output="rs.csv" 	#name of output Meteo txt file
./mkMeteo.py "$netcdf" "$RSdir" $longitude $latitude $output
