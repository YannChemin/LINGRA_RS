# LINGRA_RS
 [LINGRA][1] grass growth and production model translated in Python bypassed by Satellite Remote Sensing data and input from [ERA5][2] live weather dataset

![output](https://github.com/YannChemin/LINGRA_RS/blob/master/Figure_1.png "LINGRA_RS output")

## Keywords
grass growth, grass production, water use by soil evaporation and transpiration, sink versus source limited growth

## Description
[LINGRA][1] is simple model to calculate grass growth and production under potential and water limited conditions. 

Simulated key processes are:

```
1 - Light utilization
2 - Leaf formation
3 - Leaf elongation
4 - Tillering
5 - Carbon partitioning to roots and shoots
```

Source and sink limited growth are simulated independently. 

Sink-limited growth is characterized by temperature-dependent leaf expansion and tiller development, whereas source-limited growth is determined by photosynthetic light-use efficiency of the canopy and th remobilization of stored carbohydrates in the stubble. 

At each integration step of 1 day, the available amount of carbon from the source is compared with the carbon required by the sink. 

The actual growth is determined by the minimum value of either sink or source.

## Run

```
> python3 lingraRS.py 52.0 NL1.987.rs.csv
```

## Input file structure

```
0    year      # year in weather file
1    doy       # doy of the year
2    rdd       # solar radiation (kj m-2 day-1)
3    tmmn      # minimum temperature (degrees celsius)
4    tmmx      # maximum temperature (degrees celsius)
5    vp        # water vapour pressure (kpa)
6    wn        # average wind speed (m s-1)
7    rain      # daily rainfall (mm day-1)
8    RSevp     # daily RS Evaporation actual (mm day-1)
9    RStrn     # daily RS Transpiration actual (mm day-1)
10   RSlai     # daily RS LAI (-)
11   RScut     # daily RS cutting event (0/1)
```

## Input file example

```

```

[1]: https://models.pps.wur.nl/lingra-model-simple-grass-model-potential-and-water-limited-conditions 
[2]: https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
