import os
import glob
from netCDF4 import Dataset
import numpy as np
from datetime import date


def getyear(fname):
    return fname[-36:-32]


def getmonth(fname):
    return fname[-32:-30]


def getday(fname):
    return fname[-30:-28]


def getdoy(fname):
    """Extract doy from the file name
    :param fname: ex. BEC_SM____SMOS__EUM_L4__A_20200826T035430_001km_3d_REP_v5.0.nc
    :return: doy
    """
    year = int(getyear(fname))
    month = int(getmonth(fname))
    day = int(getday(fname))
    #print(fname, year,month,day)
    d = date(year, month, day)
    return int(d.strftime("%j"))


def getgeo(fname, longitude, latitude):
    """Get geographically bound data from NC file
    :param fname: netcdf filename (.nc)
    :param longitude: longitude dd.ddd)
    :param latitude: latitude (dd.ddd)
    :return: pixel value at long,lat
    """
    f = Dataset(fname)
    sm = f.variables['SM']
    longitude_array = f.variables['lon'][:].data
    latitude_array = f.variables['lat'][:].data
    smraw = f.variables['SM'][:].data
    i = np.abs(longitude_array - np.float(longitude)).argmin()
    j = np.abs(latitude_array - np.float(latitude)).argmin()
    pixval = smraw[:, j, i]
    if pixval == sm._FillValue or pixval == sm.missing_value:
        pixval = np.nan
    else:
        pixval = pixval
    # pixval = sm.scale_factor * pixval + sm.add_offset
    if pixval < 0.0:
        pixval = 0.0
    if pixval > 0.6:
        pixval = 0.6
    return pixval


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.
    https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        #>>> # linear interpolation of NaNs
        #>>> nans, x= nan_helper(y)
        #>>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return np.isnan(y), lambda z: z.nonzero()[0]


def getgeoydoy(rsdir, wildcard, longitude, latitude, listyear, listdoy, **kwargs):
    interpolate = kwargs.get('i', None)
    py_glob = os.path.join(rsdir, wildcard)
    #print(rsdir, wildcard)
    #print(py_glob)
    RSyear = []
    RSdoy = []
    value = []
    for path in glob.iglob(py_glob):
        year = getyear(path)
        #print(year)
        doy = getdoy(path)
        #print(doy)
        RSyear.append(year)
        RSdoy.append(doy)
        value.append(getgeo(path, longitude, latitude))
        #print(value[-1])

    resarr = np.zeros((len(listdoy),), dtype=np.float)
    resarr.fill(np.nan)
    # Fill list with nan for the len(listdoy)
    for d in range(len(listdoy)):
        for l in range(len(RSdoy)):
            if int(listyear[d]) == int(RSyear[l]):
                if int(listdoy[d]) == int(RSdoy[l]):
                    resarr[d] = value[l]
    # If Kwargs 'i' is set to True (i=True), interpolate the content
    try:
        if interpolate is True:
            # interpolate nans to values
            nans, x = nan_helper(resarr)
            resarr[nans] = np.interp(x(nans), x(~nans), resarr[~nans])
    except:
        pass

    # print(result)
    result = resarr.tolist()
    return result

#	int16 SM(time, lat, lon)
#   long_name: Surface Soil Moisture
#   units: m^3/m^3
#   scale_factor: 1e-04
#   add_offset: 0.0
#   valid_min: 0.0
#   valid_max: 0.6
#   coordinates: time lat lon
#   missing_value: -999
#   _FillValue: -999
