#!/bin/bash
#For gdal tools
export PATH=$PATH:/Library/Frameworks/GDAL.framework/Versions/3.1/Programs/
#Copernicus GRASSLAND tif file (or subset)
grassland="$HOME/Documents/GRA_2018_10m/GRA_2018_BE_WA_BRABANT.tif"
#ERA5 file access handle
netcdf="$HOME/Documents/ERA5/ERA5_EU_2020.nc"	#name of NETCDF file to access
#RS directory access handle
RSdir="$HOME/Documents/MODIS"	#name of RS LAI/ET/etc directory
#BEC SMOS 1km data dir
becsmosdir="$HOME/Documents/SMOS"
#base name of output pickled list
output="Wallonie_Brabant" 	#name of output Meteo txt file
python3 /Users/dnd/Documents/GitHub/LINGRA_RS/mkMeteoF/mkMeteoAllPix2pickle.py "$grassland" "$netcdf" "$RSdir" "$becsmosdir" "$output"

