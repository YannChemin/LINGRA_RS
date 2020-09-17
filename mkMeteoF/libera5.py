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
                'surface_solar_radiation_downwards', '10m_u_component_of_wind',
                '10m_v_component_of_wind', '2m_temperature', 'surface_pressure',
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
            ],  # North, West, South, East. Default: global
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
        out_dir + '/ERA5_' + prefix + '_' + str(i_year) + '.nc')


def getTimeNetcdf(netcdffile):
    """
    Function to extract the time component of a netCDF file
    """
    # Get the time array from NETCDF file and convert to Python datetime
    ncfile = netCDF4.Dataset(netcdffile, 'r')
    tim = ncfile.variables['time']  # do not cast to numpy array yet
    time_convert = netCDF4.num2date(tim[:], tim.units, tim.calendar, only_use_cftime_datetimes=False)
    return time_convert


def getPixelVals(netcdffile, layer, longitude, latitude):
    """
    Function to query pixels values for a given layer at a given Lat/long
    """
    nc = netCDF4.Dataset(netcdffile)
    dict_nc = nc[layer].__dict__
    # Extract the offset and scale from netcdf file
    offset = float(dict_nc['add_offset'])
    scale = float(dict_nc['scale_factor'])
    # Extract NoData Value from layer
    nodata = int(dict_nc['missing_value'])
    fillva = int(dict_nc['_FillValue'])
    print(offset, scale, nodata, fillva)

    lat, lon = nc.variables['latitude'], nc.variables['longitude']
    # extract lat/lon values (in degrees) to numpy arrays
    latvals = lat[:]
    lonvals = lon[:]
    xl = np.linspace(lonvals.min(), lonvals.max(), lonvals.shape[0])
    yl = np.linspace(latvals.min(), latvals.max(), latvals.shape[0])
    xx, yy = np.meshgrid(xl, yl, sparse=True)

    def distance_2d(longitude, latitude, xgrid, ygrid):
        return np.hypot(xgrid - longitude, ygrid - latitude)

    dist_grid = distance_2d(longitude, latitude, xx, yy)
    minval = dist_grid.min()
    for x in range(dist_grid.shape[1]):
        for y in range(dist_grid.shape[0]):
            if dist_grid[y][x] == minval:
                ix_min = x
                iy_min = y

    # Select all of the temporal instances of the pixel
    arr = nc.variables[layer][:, 0, iy_min, ix_min]
    if layer == 'ssrd':
        arr[arr == 0] = np.nan
    # Replace nodata (often -32767) with NAN
    arr[arr == nodata] = np.nan
    arr[arr == -nodata] = np.nan
    # Fill in value is sometimes -32766, but not documented...
    arr[arr == 32766] = np.nan
    # Fill values to NAN
    arr[arr == fillva] = np.nan
    arr[arr == -fillva] = np.nan
    if layer == 'ssrd':
        return arr
    # Rescale the data
    array = offset + arr * scale
    # Return the array
    return array


def plotLayer(layer, time_convert, array):
    """
    Function to plot the array data with time
    """
    fig, ax = plt.subplots()
    ax.plot_date(time_convert, array, 'g-')
    if (layer.strip() == "ssrd"):
        ylab = 'surface downwelling shortwave flux in air (J m-2)'
    elif (layer.strip() == "u10"):
        ylab = 'windspeed u component (m s-1)'
    elif (layer.strip() == "v10"):
        ylab = 'windspeed v component (m s-1)'
    elif (layer.strip() == "t2m"):
        ylab = 'T2m (K)'
    elif (layer.strip() == "d2m"):
        ylab = 'Dewpoint T2m (K)'
    #elif (layer.strip() == "sp"):
    #    ylab = 'Surface air pressure (Pa)'
    elif (layer.strip() == "tp"):
        ylab = 'Total precipitation (m)'
    else:
        ylab = 'undefined product'
    ax.set(xlabel='time', ylabel=ylab, title='ERA5 ' + ylab)
    ax.grid()
    # fig.savefig("test.png")
    plt.show()


def d2m2eact(d2m):
    """
    Converts dew point temperature to actual vapour pressure (hPa)
    """
    return (0.6108 * np.exp((17.27 * d2m) / (d2m + 237.3)))
