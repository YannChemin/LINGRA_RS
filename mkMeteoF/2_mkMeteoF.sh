#!/bin/bash
#ERA5 file access handle
netcdf="$HOME/Documents/ERA5/ERA5_EU_2020.nc"	#name of NETCDF file to access
#RS directory access handle
RSdir="$HOME/Documents/MODIS"	#name of RS LAI/ET/etc directory
#GRASSLAND master RS file original handle
grasslandF="$HOME/Documents/GRA_2018_10m/GRA_2018_100m_eu_03035_v010/DATA/GRA_2018_100m_eu_03035_V1_0.tif"


#############################################################
#Dealing with the grassland pixels location from Copernicus
#############################################################
#Cut to ERA5_EU_2020 boundaries and reproj to EPSG:4326
#For gdalwarp
export PATH=$PATH:/Library/Frameworks/GDAL.framework/Versions/3.1/Programs/
#Output file handle
gF="$HOME/Documents/GRA_2018_10m/GRA_2018_EU.tif"
#rm -f $gF
if [ -f $gF ]
then
  echo "$gF already exists, not overwriting."
else
  # EU boundaries from ERA5
  north=60.0
  south=35.0
  east=30.0
  west=-10.0
  #Cut and reproj (and set to 0/1)
  # -co "GDAL_DISABLE_READDIR_ON_OPEN=TRUE"
  # -co "GDAL_CACHEMAX=1000"
  gdalwarp -q -srcnodata "0 254 255" -dstnodata "0" -ot "INT32" -co "NUM_THREADS=ALL_CPUS" -co "COMPRESS=DEFLATE" -co "TILED=YES" -t_srs "EPSG:4326" -te $west $south $east $north $grasslandF $gF
fi

# For DEBUG
#from osgeo import gdal
#grass=gdalopen($gF)
#ForLoop
  #longitude=47.5
  #latitude=-0.2
  #output="rs.csv" 	#name of output Meteo txt file
  #./mkMeteo.py "$netcdf" "$RSdir" $longitude $latitude $output
#end for loop