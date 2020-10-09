#!/usr/bin/env python3

def mkmeteo4lingrars(netcdf, rsdir, longitude, latitude):
    """
    Prepare meteo and RS constraints for lingraRS model
    Requirements for local pylibs: libera5, libgetgeo
    Requirements for python libs: os, sys, netCDF4, numpy, datetime, argparse
    :param netcdf: input ERA5 map weather data
    :param rsdir: input directory where is RS data
    :param longitude: longitude (dd.ddd)
    :param latitude: latitude (dd.ddd)
    :return: list of lists of input parameters
    --------------------------------------------------------
    Column  Daily value
    1       year
    2       day
    3       irradiation                   (kJ m-2 d-1)
    4       minimum temperature           (degrees Celsius)
    5       maximum temperature           (degrees Celsius)
    6       early morning vapour pressure (kPa)
    7       mean wind speed (height: 10m) (m s-1)
    8       precipitation                 (mm d-1)
    9       RS evaporation                (mm d-1)
    10      RS transpiration              (mm d-1)
    11      RS LAI                        (m2.m-2)
    12      RS cut event                  (0/1)
    --------------------------------------------------------
    """
    import datetime
    import numpy as np
    import xarray as xr
    # Import local libraries
    import libera5
    import libgetgeo

    # Extract time array in python datetime
    time_convert = libera5.getTimeNetcdf(netcdf)

    # Set the basic array generation loop
    # layers = ['ssrd', 'u10', 'v10', 't2m', 'sp', 'tp', 'd2m']
    layers = ['ssrd', 'u10', 'v10', 't2m', 'tp', 'd2m']
    for ll in range(len(layers)):
        arr = libera5.getPixelVals(netcdf, layers[ll], longitude, latitude)
        if layers[ll] == 'ssrd':
            ssrd = np.copy(arr).astype('float64')
        elif layers[ll] == 'u10':
            u10 = np.copy(arr).astype('float64')
        elif layers[ll] == 'v10':
            v10 = np.copy(arr).astype('float64')
        elif layers[ll] == 't2m':
            t2m = np.copy(arr).astype('float64')
        # elif layers[ll] == 'sp':
        #   sp = np.copy(arr).astype('float64')
        elif layers[ll] == 'tp':
            tp = np.copy(arr).astype('float64')
        elif layers[ll] == 'd2m':
            d2m = np.copy(arr).astype('float64')
        else:
            print("name of this dataset is unknown")

    # Make Wind speed (u10 x v10 components)
    wspd = np.abs(np.multiply(u10, v10))
    # Solar Radiation J.m-2.s
    # Assuming rad per s so /3600 for an hour
    ssrd = np.divide(ssrd, 3600)
    # Surface Pressure Pa -> hPa
    # sp = np.divide(sp, 100)
    # Temperature K -> C
    t2m = np.subtract(t2m, 273.15)
    # Precipitation m -> mm
    tp = np.abs(np.multiply(tp, 1000))
    # Temperature K -> C
    d2m = np.subtract(d2m, 273.15)

    # Before daily conversion, remove night values for some
    # Replace nodata with NAN
    ssrd[ssrd == 0.0000] = np.nan

    # Convert from hour to day
    dsas = xr.Dataset(
        {
            "ssrd": ("time", ssrd),
        },
        {"time": time_convert},
    )
    # Make mean on non NAN values
    # dsasD=dsas.resample(time='1D').mean(skipna=True)
    # Make mean on non NAN values
    dsasD = dsas.resample(time='1D').reduce(np.nanmean)

    dsa = xr.Dataset(
        {
            "wspd": ("time", wspd),
        },
        {"time": time_convert},
    )
    dsaD = dsa.resample(time='1D').mean()

    # Precipitation summed by day
    dsap = xr.Dataset(
        {
            "prmm": ("time", tp),
        },
        {"time": time_convert},
    )
    dsapD = dsap.resample(time='1D').sum()

    # Minimum temperature per day
    dstmin = xr.Dataset(
        {
            "tmin": ("time", t2m),
        },
        {"time": time_convert},
    )
    dstminD = dstmin.resample(time='1D').min()

    # Maximum temperature per day
    dstmax = xr.Dataset(
        {
            "tmax": ("time", t2m),
        },
        {"time": time_convert},
    )
    dstmaxD = dstmax.resample(time='1D').max()

    # Minimum dew point temperature per day
    dsdmin = xr.Dataset(
        {
            "dmin": ("time", d2m),
        },
        {"time": time_convert},
    )
    dsdminD = dsdmin.resample(time='1D').min()

    # Convert Dewpoint T Min (C) to Eact (hPa)
    eactt = np.array(dsdminD.dmin)
    eact = []
    for d in range(len(eactt)):
        eact.append(libera5.d2m2eact(eactt[d]))

    # Generate year value for every day
    year = np.datetime_as_string(dsaD.time, unit='Y')
    # Generate DOY value for every day
    doyt = np.datetime_as_string(dsaD.time, unit='D')
    doy = []
    for t in range(len(doyt)):
        doy.append(datetime.datetime.strptime(doyt[t], '%Y-%m-%d').timetuple().tm_yday)

    # Join all data
    meteolist = [year.tolist(), doy, np.array(dsasD.ssrd).tolist(), np.array(dstminD.tmin).tolist(),
                 np.array(dstmaxD.tmax).tolist(), eact, np.array(dsaD.wspd).tolist(), np.array(dsapD.prmm).tolist()]
    # INSERT RS DATA HERE
    ####################
    # Create E & T from RS
    # Get the FC files in numpy array converted to [0.0-1.0]
    FC = libgetgeo.getgeoydoy(rsdir, "MOD13Q1*_FC.tif", longitude, latitude, year.tolist(), doy, i=True)
    FC = np.asarray(FC, dtype=np.float)
    FC = np.divide(FC, 100)
    # print("START FC", FC, "END FC")
    # Get ETa from both MOD and MYD platforms in arrays
    MODET = libgetgeo.getgeoydoy(rsdir, "MOD16A*.vrt", longitude, latitude, year.tolist(), doy)
    MYDET = libgetgeo.getgeoydoy(rsdir, "MYD16A*.vrt", longitude, latitude, year.tolist(), doy)
    # make 32765 -> NAN in MODET/MYDET
    MODET = np.asarray(MODET, dtype=np.float)
    MYDET = np.asarray(MYDET, dtype=np.float)
    MODET[MODET == 32765] = np.nan
    MODET[MODET == 32766] = np.nan
    MODET[MODET == 32767] = np.nan
    MYDET[MYDET == 32765] = np.nan
    MYDET[MYDET == 32766] = np.nan
    MYDET[MYDET == 32767] = np.nan
    MODET = np.divide(MODET, 10)
    MYDET = np.divide(MYDET, 10)
    # print("START MODET", MODET, "END MODET")
    # print("START MYDET", MYDET, "END MYDET")
    # Compute average ETa. Take NANs into account to maximise data count out
    ET = (MODET + MYDET) / 2.0
    # print("START ET", ET, "END ET")
    # Transpiration
    T = np.multiply(ET, FC)
    # print("START T", T, "END T")
    # Evaporation is the complementary of Transpiration vis-a-vis ET
    E = np.subtract(ET, T)
    # print("START E", E, "END E")
    # Add Evaporation data to the Meteo file
    meteolist.append(E.tolist())
    # Add Transpiration data to the Meteo file
    meteolist.append(T.tolist())
    # Add LAI data to the Meteo file
    LAI = libgetgeo.getgeoydoy(rsdir, "MCD15A2H*.vrt", longitude, latitude, year.tolist(), doy)
    LAI = np.array(LAI, dtype=np.float)
    LAI = np.divide(LAI, 10)
    # print("START LAI", LAI, "END LAI")
    meteolist.append(LAI.tolist())
    # Create artificial Cut data from a previous array
    CutNans = np.zeros(LAI.shape[0])
    # Fill with NANs
    CutNans.fill(np.nan)
    # Fill with 2 cuts
    CutNans[127] = 1
    CutNans[235] = 1
    # Add Cut data to the Meteo file
    meteolist.append(CutNans.tolist())
    # print("START CUTNANS", CutNans, "END CUTNANS")
    return meteolist
