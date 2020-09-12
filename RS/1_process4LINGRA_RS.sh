#!/bin/bash

export PATH=$PATH:/Applications/QGIS3.14.app/Contents/MacOS/bin
export PATH=/usr/local/bin:$PATH
# Needed to debug on Mac
# ln -s /usr/local/Cellar/gcc/9.3.0_1/lib/gcc/9 /usr/local/opt/gcc/lib/gcc/9

#Run this from the RSDATA/ dir
RSDTRoot=/Users/dnd/Documents/MODIS

#Clean & Make
progRoot=$HOME/Documents/GitHub/distRS/prog
cd $progRoot/prog_fc
make clean
make
cd $progRoot/prog_NDVI
make clean
make
cd $progRoot/prog_LAI
make clean
make
cd $progRoot/prog_ETa_MODIS
make clean
make

#GO TO PROCESSING DIRECTORY
cd $RSDTRoot

#Clean up previous jobs
mkdir -p xml
mv *.xml xml/

#PROCESS NDVI for VegFraction (split E&T after .vrt are made)
for file in MOD13Q1*.hdf
do
  NDVIF=$(gdalinfo "$file" | grep SUBDATASET_1_NAME | sed "s/\ \ SUBDATASET_1_NAME=//")
  NDVIOUT=$(echo "$NDVIF" | sed "s/HDF4_EOS:EOS_GRID:\"\(.*\)\":\(.*\)/\1/" | sed "s/.hdf/_NDVI.tif/")
  FEXISTS=$(test -f "$NDVIOUT" && echo 1 || echo 0)
  if [ -f "$NDVIOUT" ];
  then
    echo "$NDVIF"
    echo "$NDVIOUT exists, Not Overwriting.";
  else
        echo "$NDVIF"
        #Offset=$(gdalinfo $NDVIF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\1/")
        #Scale=$(gdalinfo $NDVIF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\2/")
        NDVIFQC=$(gdalinfo "$file" | grep SUBDATASET_12_NAME | sed "s/\ \ SUBDATASET_12_NAME=//")
        $progRoot/prog_NDVI/ndvi "$NDVIF" "$NDVIFQC" "$NDVIOUT" #$Offset $Scale
        #TODO integrate -co into distRS code directly
        gdal_translate -q -co "COMPRESS=DEFLATE" -co "TILED=YES" "$NDVIOUT" tmp"$NDVIOUT"
        mv tmp"$NDVIOUT" "$NDVIOUT";
        rm -f *.aux.xml &
  fi
done

#PROCESS LAI
for file in MCD15A2H*.hdf
do
  LAIF=$(gdalinfo "$file" | grep SUBDATASET_2_NAME | sed "s/\ \ SUBDATASET_2_NAME=//")
  LAIOUT=$(echo "$LAIF" | sed "s/\(.*\)\"\(.*\)\"\(.*\)/\2/" | sed "s/.hdf/_LAI.tif/")
  if [ -f "$LAIOUT" ]
  then
    echo "$LAIOUT exists, Not Overwriting."
  else
    echo "$LAIF"
    #Offset=$(gdalinfo $LAIF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\1/")
    #Scale=$(gdalinfo $LAIF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\2/")
    LAIFQC=$(gdalinfo "$file" | grep SUBDATASET_3_NAME | sed "s/\ \ SUBDATASET_3_NAME=//")
    $progRoot/prog_LAI/lai "$LAIF" "$LAIFQC" "$LAIOUT" #$Offset $Scale
    #TODO integrate -co into distRS code directly
    gdal_translate -q -co "COMPRESS=DEFLATE" -co "TILED=YES" "$LAIOUT" tmp"$LAIOUT"
    mv tmp"$LAIOUT" "$LAIOUT"
    rm -f *.aux.xml &
  fi
done

#PROCESS MOD ETA
for file in MOD16A2*.hdf
do
  ETAF=$(gdalinfo "$file" | grep SUBDATASET_1_NAME | sed "s/\ \ SUBDATASET_1_NAME=//")
  ETAOUT=$(echo "$ETAF" | sed "s/\(.*\)\"\(.*\)\"\(.*\)/\2/" | sed "s/.hdf/_ETA.tif/")
  if [ -f "$ETAOUT" ]
  then
    echo "$ETAOUT exists, Not Overwriting."
  else
    echo "$ETAF"
    #Offset=$(gdalinfo $ETAF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\1/")
    #Scale=$(gdalinfo $ETAF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\2/")
    ETAFQC=$(gdalinfo "$file" | grep SUBDATASET_5_NAME | sed "s/\ \ SUBDATASET_5_NAME=//")
    $progRoot/prog_ETa_MODIS/eta "$ETAF" "$ETAFQC" "$ETAOUT" #$Offset $Scale
    #TODO integrate -co into distRS code directly
    gdal_translate -q -co "COMPRESS=DEFLATE" -co "TILED=YES" "$ETAOUT" tmp"$ETAOUT"
    mv tmp"$ETAOUT" "$ETAOUT"
  fi
done

#PROCESS MYD ETA
for file in MYD16A2*.hdf
do
  ETAF=$(gdalinfo "$file" | grep SUBDATASET_1_NAME | sed "s/\ \ SUBDATASET_1_NAME=//")
  ETAOUT=$(echo "$ETAF" | sed "s/\(.*\)\"\(.*\)\"\(.*\)/\2/" | sed "s/.hdf/_ETA.tif/")
  if [ -f "$ETAOUT" ]
  then
    echo "$ETAOUT exists, Not Overwriting."
  else
    echo $ETAF
    #Offset=$(gdalinfo $ETAF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\1/")
    #Scale=$(gdalinfo $ETAF | grep Offset:\ | sed "s/\ \ Offset:\ \(.*\),\ \ \ Scale:\(.*\)/\2/")
    ETAFQC=$(gdalinfo "$file" | grep SUBDATASET_5_NAME | sed "s/\ \ SUBDATASET_5_NAME=//")
    $progRoot/prog_ETa_MODIS/eta $ETAF $ETAFQC $ETAOUT #$Offset $Scale
    #TODO integrate -co into distRS code directly
    gdal_translate -q -co "COMPRESS=DEFLATE" -co "TILED=YES" $ETAOUT tmp$ETAOUT
    mv tmp$ETAOUT $ETAOUT
  fi
done

#####################################
#Build VRT by MODIS product by date #
#####################################

#########################################################################
#TODO make an auto detect of years available and define startYYYY endYYYY
#startd=$startYYYY"001"
#endd=$endYYYY"366"
###################

startd=2020001
endd=2020366

# Overwrite by default : Clean slate
OVR='-overwrite'

for ((date = $startd; date <= $endd; date++))
do
  MODIS='MOD13Q1'
  l=$(find $MODIS.A$date*_NDVI.tif 2>/dev/null | wc -l)
  if [ $l -gt 2 ]; then
    gdalbuildvrt -q $OVR -srcnodata 255 -vrtnodata 255 $MODIS\_$date.vrt $(ls $MODIS*$date*_NDVI.tif)
    #Create Fraction of Vegetation Cover
    FCOUT=$(echo $MODIS\_$date.vrt | sed 's/.vrt/_FC.tif/')
    gdal_translate -q -co "COMPRESS=DEFLATE" -co "TILED=YES" $MODIS\_$date.vrt $FCOUT
    $progRoot/prog_fc/fc $FCOUT tmp$FCOUT
    mv tmp$FCOUT $FCOUT
  fi
  MODIS='MCD15A2H'
  l=$(find $MODIS.A$date*_LAI.tif 2>/dev/null | wc -l)
  if [ $l -gt 2 ]; then
    gdalbuildvrt -q $OVR -srcnodata 255 -vrtnodata 255 $MODIS\_$date.vrt $(ls $MODIS*$date*_LAI.tif) &
  fi
  MODIS='MOD16A2'
  l=$(find $MODIS.A$date*_ETA.tif 2>/dev/null | wc -l)
  if [ $l -gt 2 ]; then
    gdalbuildvrt -q $OVR -srcnodata 32767 -vrtnodata 32767 $MODIS\_$date.vrt $(ls $MODIS*$date*_ETA.tif) &
  fi
  MODIS='MYD16A2'
  l=$(find $MODIS.A$date*_ETA.tif 2>/dev/null | wc -l)
  if [ $l -gt 2 ]; then
    gdalbuildvrt -q $OVR -srcnodata 32767 -vrtnodata 32767 $MODIS\_$date.vrt $(ls $MODIS*$date*_ETA.tif) &
  fi
done

mv xml/*.xml .
