import os
import glob
import numpy as np


def getgeo(fname, longitude, latitude):
    return os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (fname, longitude, latitude)).read()


def getdoy(fname):
    if fname[-7:-4] == "_FC":
        return fname[-10:-7]
    else:
        return fname[-7:-4]


def getyear(fname):
    if fname[-7:-4] == "_FC":
        return fname[-14:-10]
    else:
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
        value.append(getgeo(path, longitude, latitude).split('\n')[0])

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
