import os
import glob
import numpy as np


def getgeo(fname, longitude, latitude):
    return os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (fname, longitude, latitude)).read()


def getdoy(fname):
    return fname[-7:-4]


def getyear(fname):
    return fname[-11:-7]


def getgeoWLDCRD(VRTdir, wildcard, longitude, latitude):
    py_glob = os.path.join(VRTdir, wildcard)
    result = []
    for path in glob.iglob(py_glob):
        result.append(getyear(path), getdoy(path), getgeo(path, longitude, latitude))
    return result


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

def getgeoydoy(VRTdir, wildcard, longitude, latitude, listyear, listdoy, **kwargs):
    interpolate = kwargs.get('i', None)
    py_glob = os.path.join(VRTdir, wildcard)
    RSyear = []
    RSdoy = []
    value = []
    for path in glob.iglob(py_glob):
        year = getyear(path)
        doy = getdoy(path)
        RSyear.append(year)
        RSdoy.append(doy)
        value.append(year, doy, getgeo(path, longitude, latitude))

    result = []
    # Fill list with nan for the len(listdoy)
    for d in range(len(listdoy)):
        for l in range(len(RSdoy)):
            if listyear[d] == RSyear[l]:
                if listdoy[d] == RSdoy[l]:
                    result.append(value[l])
            else:
                result.append('nan')

    # If Kwargs 'i' is set to True (i=True), interpolate the content
    if interpolate is True:
        resarr = np.array(result, dtype=np.float)
        # interpolate nans to values
        nans, x = nan_helper(resarr)
        resarr[nans] = np.interp(x(nans), x(~nans), resarr[~nans])
        result = resarr.tolist()

    return result
