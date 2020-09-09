#!/bin/bash

#TODO Replace this by a Lat/Lon iterator based on 1x1Km grassland pixels from Copernicus 2015 

netcdf="/media/yann/KINGSTON/ERA5/ERA5_EU_2019.nc"	#name of NETCDF file to access
RSdir="/media/yann/KINGSTON/MODIS"	#name of RS LAI/ET/etc directory
longitude=47.5
latitude=-0.2 
output="rs.csv" 	#name of output Meteo txt file
./mkMeteo.py $netcdf $RSdir $longitude $latitude $output
