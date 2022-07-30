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
1987		105		15140			4.1	  15	   0.95	 2.3	 0	   2.08	 1.8	 nan	 nan
1987		106		8520				7.5	  15.1	 1.1	  1.4	 0	   nan  	nan	 nan	 nan
1987		107		12270	  9.2	  17.1	 1.07	 1.2	 0	   nan	  nan	 nan	 nan
1987		108		18330	  8.6	  21.4	 1.14	 2.2	 0	   nan	  nan	 nan	 nan
1987		109		7260	   8.7	  16.7	 1.21	 2.6	 3.6	 nan	  nan	 nan	 nan
1987		110		9750	   5.2	  13.8	 0.94	 5.6	 0.9	 nan	  nan	 nan	 nan
1987		111		17220	  2.8	  13.4	 0.8	  2.1	 0	   nan	  nan	 4.5	 nan
1987		112		21460	  -0.5	 15.9	 0.78	 1.6	 0	   nan	  nan	 nan	 nan
1987		113		21250	  5.6	  19.4	 0.89	 1.7	 0	   nan	  nan	 nan	 nan
1987		114		21820	  3.3	  22.1	 0.87	 1.9	 0	   nan	  nan	 nan	 nan
1987		115		19030	  5.3	  23.5	 1.07	 1.5	 0	   nan	  nan	 nan	 nan
1987		116		18970	  5.4	  19.7	 0.92	 1.6	 0	   nan	  nan	 nan	 nan
1987		117		23270	  2.9	  18.7	 0.73	 2.3	 0	   nan	  nan	 nan	 nan
1987		118		22820	  8.3	  22.5	 0.74	 3	   0	   nan	  nan	 nan	 nan
1987		119		16850	  12.5	 23.9	 1.17	 2.6	 0	   nan	  nan	 4.8	 1
1987		120		13470	  11.5	 21.2	 1.34	 3.1	 0.2	 nan	  nan	 nan	 nan
1987		121		3400	   11.6	 14.8	 1.3	  3.6	 0.1	 nan	  nan	 nan	 nan
1987		122		16810	  3.8	  13.7	 0.84	 3.5	 0.3	 2.25	 4.8	 nan	 nan
1987		123		9280	   3.6	  9.4	  0.77	 2.8	 3.8	 nan	  nan	 nan	 nan
```

[1]: https://models.pps.wur.nl/lingra-model-simple-grass-model-potential-and-water-limited-conditions 
[2]: https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
