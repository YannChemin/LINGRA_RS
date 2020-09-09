import os
import netCDF4
import cdsapi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def download_ERA5_bbox(i_year, out_dir, prefix, north, west, south, east):
    c = cdsapi.Client()
    c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            'surface_solar_radiation_downwards','10m_u_component_of_wind',
            '10m_v_component_of_wind','2m_temperature','surface_pressure',
            'total_precipitation', '2m_dewpoint_temperature'
            # 'evaporation', 'potential_evaporation',
            # 'snow_density', 'snow_depth'
            # '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            # '2m_temperature', 'angle_of_sub_gridscale_orography', 
            # 'anisotropy_of_sub_gridscale_orography',
            # 'evaporation', 'land_sea_mask', 
            # 'maximum_total_precipitation_rate_since_previous_post_processing',
            # 'mean_evaporation_rate', 'mean_potential_evaporation_rate',
            # 'minimum_total_precipitation_rate_since_previous_post_processing',
            # 'orography', 'potential_evaporation', 'precipitation_type',
            # 'snow_density', 'snow_depth', 'soil_temperature_level_1',
            # 'soil_temperature_level_2', 'soil_temperature_level_3', 'soil_temperature_level_4',
            # 'soil_type', 'surface_pressure', 'surface_solar_radiation_downwards',
            # 'total_column_snow_water', 'total_precipitation', 'volumetric_soil_water_layer_1',
            # 'volumetric_soil_water_layer_2', 'volumetric_soil_water_layer_3', 
            # 'volumetric_soil_water_layer_4',
        ],
        'area': [
            north, west, south, east,
        ], # North, West, South, East. Default: global
        'year': str(i_year),
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'format': 'netcdf',
    },
    out_dir + '/ERA5_'+ prefix +'_' + str(i_year) +'.nc')

def getTimeNetcdf(netcdfFile):
    """
    Function to extract the time component of a netCDF file
    """
    #Get the time array from NETCDF file and convert to Python datetime
    ncfile = netCDF4.Dataset(netcdfFile, 'r')
    tim = ncfile.variables['time'] # do not cast to numpy array yet
    time_convert = netCDF4.num2date(tim[:], tim.units, tim.calendar,only_use_cftime_datetimes=False)
    return(time_convert)

def getPixelVals(netcdfFile, layer, longitude, latitude):
    """
    Function to query pixels values for a given layer at a given Lat/long
    """
    #Construct the Layer name for GDAL
    lyr = "NETCDF:"+netcdfFile+":"+layer
    #print(lyr)

    #Parse GDALINFO for the needed
    #Define UR + Resolution
    URx=os.popen('gdalinfo %s | grep Origin | sed "s/\(.*\)(\(.*\),\(.*\))/\\2/"' % (lyr)).read()
    URy=os.popen('gdalinfo %s | grep Origin | sed "s/\(.*\)(\(.*\),\(.*\))/\\3/"' % (lyr)).read()
    #print(URx,URy)
    Rsx=os.popen('gdalinfo %s | grep Pixel\ Size | sed "s/\(.*\)(\(.*\),\(.*\))/\\2/"' % (lyr)).read()
    Rsy=os.popen('gdalinfo %s | grep Pixel\ Size | sed "s/\(.*\)(\(.*\),-\(.*\))/\\3/"' % (lyr)).read()
    #print(Rsx,Rsy)
    #Extract the offset and scale from netcdf file
    offset=os.popen('gdalinfo %s | grep Offset | tail -n1 | sed "s/\(.*\):\ \(.*\),\(.*\):\(.*\)/\\2/"' % (lyr)).read()
    scale=os.popen('gdalinfo %s | grep Offset | tail -n1 | sed "s/\(.*\):\ \(.*\),\(.*\):\(.*\)/\\4/"' % (lyr)).read()
    #print(offset,scale)
    #Extract NoData Value from layer
    nodata=os.popen('gdalinfo %s | grep NoData\ Value | tail -n1 | sed "s/\ \ NoData\ Value=//"' % (lyr)).read()

    #Get row and column numbers (Not needed at this point)
    #nX=os.popen('gdalinfo %s | grep Size\ is\ | sed "s/Size\ is\ \(.*\),\(.*\)/\\1/"' % (lyr)).read()
    #nY=os.popen('gdalinfo %s | grep Size\ is\ | sed "s/Size\ is\ \(.*\),\(.*\)/\\2/"' % (lyr)).read()
    #print(nX,nY)
    
    #Clean vars
    URx=float(URx.strip())
    URy=float(URy.strip())
    Rsx=float(Rsx.strip())
    Rsy=float(Rsy.strip())
    offset=float(offset.strip())
    scale=float(scale.strip())
    #print(offset,scale)
    nodata=int(nodata.strip())

    #Convert from Lat/Long to X/Y
    lon=float(longitude)
    lat=float(latitude)
    #lon=URx+X*Rsx
    X=int((lon-URx)/Rsx)
    #lat=URy-Y*Rsy
    Y=int((URy-lat)/Rsy)
    #print(X,Y)

    #Create GDALLOCATION X Y Var
    loc = str(X)+" "+str(Y)

    #Query the temporal values at pixel location in image column and row
    result = os.popen('gdallocationinfo -valonly %s %s' % (lyr, loc)).read()
    #NETCDF driver does not read in projected lon lat (Not needed)
    #result = os.popen('gdallocationinfo -valonly -wgs84 %s %s' % (lyr, loc)).read()

    #cleanup the \n everywhere and remove empty elements
    result1=list(result.split("\n"))
    while '' in result1:
        result1.remove('')

    #Create and fill a Numpy array
    array=np.zeros(len(result1))
    for i in range (len(result1)):
        try:
            array[i]=float(result1[i])
        except:
            print("###%s###" % (result1[i]))

    #Replace nodata with NAN
    [np.nan if x==nodata in x else x for x in array]
    #Rescale the data
    array=offset+array*scale
    #Return the array
    return(array)

def plotLayer(layer, time_convert, array):
    """
    Function to plot the array data with time
    """
    fig, ax = plt.subplots()
    ax.plot_date(time_convert,array, 'g-')
    if(layer.strip() == "ssrd"):
        ylab='surface downwelling shortwave flux in air (J m-2)'
    elif(layer.strip() == "u10"):
        ylab='windspeed u component (m s-1)'
    elif(layer.strip() == "v10"):
        ylab='windspeed v component (m s-1)'
    elif(layer.strip() == "t2m"):
        ylab='T2m (K)'
    elif(layer.strip() == "d2m"):
        ylab='Dewpoint T2m (K)'
    elif(layer.strip() == "sp"):
        ylab='Surface air pressure (Pa)'
    elif(layer.strip() == "tp"):
        ylab='Total precipitation (m)'
    else:
        ylab='undefined product'
    ax.set(xlabel='time',ylabel=ylab,title='ERA5 '+ylab)
    ax.grid()
    #fig.savefig("test.png")
    plt.show()

def d2m2eact(d2m):
    """
    Converts dew point temperature to actual vapour pressure (hPa)
    """
    return(0.6108*np.exp((17.27*d2m)/(d2m+237.3)))

