#!/bin/bash

#Specify location of the data drive
#ROOTDIR=/media/yann/KINGSTON/
ROOTDIR=/Users/dnd/Documents/

##############################################
# ERA5 DATA Download                         #
#This code requires local downERA5.py to run #
##############################################

# EU
north=60.0
south=35.0
east=30.0
west=-10.0
# Year of Interest
year=2021

outdir=$ROOTDIR/ERA5
mkdir -p $outdir

prefix=EU

# Fix the "IDontKnowYou" bug
echo "url: https://cds.climate.copernicus.eu/api/v2" > $HOME/.cdsapirc
echo "key: 46111:c87052f5-b0c1-4560-b17a-b216348961fe" >> $HOME/.cdsapirc

# Download ERA5 Data for EU
./downERA5.py $north $south $east $west $year $outdir $prefix

####################################
# RS DATA Download                 #
#This code requires pymodis to run #
####################################
#GO TO RSDATA DIR
mkdir -p $ROOTDIR/MODIS
cd $ROOTDIR/MODIS
echo $PWD

export PATH=$PATH:/Library/Frameworks/GDAL.framework/Versions/3.1/Programs/
# Try to fix ssl certif error in pymodis
# pip3 install --upgrade certifi --force-reinstall --no-cache-dir
# Add in downmodis.py line 396 397
#  import ssl
#  ssl._create_default_https_context = ssl._create_unverified_context

machine='e4ftl01.cr.usgs.gov'
user='dr.yann.chemin'
password='Gipe-74321'

# Enforce authorization in host (MacOSX is picky)
echo "machine $machine login $user password $password" > $HOME/.netrc
echo "machine urs.earthdata.nasa.gov login $user password $password" >> $HOME/.netrc

startdate=$year'-01-01'
enddate=$year'-09-01'

for hv in h15v05 h16v02 h16v05 h16v06 h17v02 h17v03 h17v04 h17v05 h17v06 h18v02 h18v03 h18v04 h18v05 h19v02 h19v03 h19v04 h19v05 h20v03 h20v04 h20v05
#for hv in h18v02 h18v03 h18v04 h18v05 h19v02 h19v03 h19v04 h19v05 h20v03 h20v04 h20v05
do
	#LAI/fAPAR 8-days 500m
	#https://e4ftl01.cr.usgs.gov/MOTA/MCD15A2H.006/
	#2020.01.01/
	modis_download.py -U $user -P $password -r -t $hv -s MOTA -p MCD15A2H.006 -f $startdate -e $enddate .
	#ET MODIS TERRA 8-days 500m
	modis_download.py -U $user -P $password -r -t $hv -p MOD16A2.006 -f $startdate -e $enddate .
	#ET MODIS AQUA 8-days 500m
	modis_download.py -U $user -P $password -r -t $hv -s MOLA -p MYD16A2.006 -f $startdate -e $enddate .
	#NDVI MODIS TERRA 16-days 250m
	modis_download.py -U $user -P $password -r -t $hv -p MOD13Q1.006 -f $startdate -e $enddate .
done

#######################################################
#TODO Build meteo data variables for input to LINGRA  #
#./getPixValNcdf.py ERA5_EU_2019.nc t2m 0.5 47.5 -p

