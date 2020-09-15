#/1usr/bin/env python3
Requirements=
"""
*---------------------------------------------------------------------------*
* Requirements: Input file grassland from Copernicus (modified 0/1)
* Requirements for python libs: os, sys, gdal, argparse
*---------------------------------------------------------------------------*
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("grassland", help="name of Copernicus Grassland file to access")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    print(Requirements)
    sys.exit(1)
args = parser.parse_args()

from osgeo import gdal

driver = gdal.GetDriverByName('GTiff')
filename = args.grassland #path to raster
dataset = gdal.Open(filename)
band = dataset.GetRasterBand(1)

cols = dataset.RasterXSize
rows = dataset.RasterYSize

transform = dataset.GetGeoTransform()

xOrigin = transform[0]
yOrigin = transform[3]
pixelWidth = transform[1]
pixelHeight = -transform[5]

data = band.ReadAsArray(0, 0, cols, rows)

#points_list = [ (355278.165927, 4473095.13829), (355978.319525, 4472871.11636) ] #list of X,Y coordinates

for col in pixelWidth:
    for row in pixelHeight:
        if data[col][row]==1:
            longitude=col*pixelWidth+xOrigin
            latitude=yOrigin-row*pixelHeight
            print(longitude, latitude, data[row][col])

#for point in points_list:
#    col = int((point[0] - xOrigin) / pixelWidth)
#    row = int((yOrigin - point[1] ) / pixelHeight)
#    print row,col, data[row][col]