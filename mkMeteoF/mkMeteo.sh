#!/bin/bash

netcdf=$1
rsdir=$2
becsmosdir=$3
longitude=$4
latitude=$5
output=$6
pickle=$7

#Define number of (virtual) cores in Linux
ncores=8 # For %$#@ Mac
#ncores=`grep -c 'processor' /proc/cpuinfo`
echo "ncores=" $ncores

#Define number of mkMeteo.py running
npid=$(echo "$(ps aux | grep mkMeteo.py | wc -l) - 1" | bc)
while [ $npid -ge $ncores ]
do
	sleep 10
	#Update number of mkMeteo.py running
	npid=$(echo "$(ps aux | grep mkMeteo.py | wc -l) - 1" | bc)
	#Update number of (virtual) cores in Linux
	#ncores=`grep -c 'processor' /proc/cpuinfo`
done

python3 /Users/dnd/Documents/GitHub/LINGRA_RS/mkMeteoF/mkMeteo.py $netcdf $rsdir $becsmosdir $longitude $latitude $output $pickle &


