#!/usr/bin/env python3

def mkmeteo4lingrars(netcdf, rsdir, becsmosdir, longitude, latitude):
    """
    Prepare meteo and RS constraints for lingraRS model
    Requirements for local pylibs: libera5, libgetgeo
    Requirements for python libs: os, sys, netCDF4, numpy, datetime, argparse
    :param netcdf: input ERA5 map weather data
    :param rsdir: input directory where is RS data
    :param becsmosdir: input directory where is BEC-SMOS soil moisture RS data
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
    11      RS lai                        (m2.m-2)
    12      RS cut event                  (0/1)
    13      RS soil moisture              (cm3/cm3)
    --------------------------------------------------------
    """
    import datetime
    import numpy as np
    import xarray as xr
    # Import local libraries
    import libera5
    import libgetgeo
    import libgetgeobecsmos

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
    # dsasd=dsas.resample(time='1D').mean(skipna=True)
    # Make mean on non NAN values
    dsasd = dsas.resample(time='1D').reduce(np.nanmean)

    dsa = xr.Dataset(
        {
            "wspd": ("time", wspd),
        },
        {"time": time_convert},
    )
    dsad = dsa.resample(time='1D').mean()

    # Precipitation summed by day
    dsap = xr.Dataset(
        {
            "prmm": ("time", tp),
        },
        {"time": time_convert},
    )
    dsapd = dsap.resample(time='1D').sum()

    # Minimum temperature per day
    dstmin = xr.Dataset(
        {
            "tmin": ("time", t2m),
        },
        {"time": time_convert},
    )
    dstmind = dstmin.resample(time='1D').min()

    # Maximum temperature per day
    dstmax = xr.Dataset(
        {
            "tmax": ("time", t2m),
        },
        {"time": time_convert},
    )
    dstmaxd = dstmax.resample(time='1D').max()

    # Minimum dew point temperature per day
    dsdmin = xr.Dataset(
        {
            "dmin": ("time", d2m),
        },
        {"time": time_convert},
    )
    dsdmind = dsdmin.resample(time='1D').min()

    # Convert Dewpoint T Min (C) to Eact (hPa)
    eactt = np.array(dsdmind.dmin)
    eact = []
    for d in range(len(eactt)):
        eact.append(libera5.d2m2eact(eactt[d]))

    # Generate year value for every day
    year = np.datetime_as_string(dsad.time, unit='Y')
    # Generate DOY value for every day
    doyt = np.datetime_as_string(dsad.time, unit='D')
    doy = []
    for t in range(len(doyt)):
        doy.append(datetime.datetime.strptime(doyt[t], '%Y-%m-%d').timetuple().tm_yday)

    # Join all data
    meteolist = [year.tolist(), doy, np.array(dsasd.ssrd).tolist(), np.array(dstmind.tmin).tolist(),
                 np.array(dstmaxd.tmax).tolist(), eact, np.array(dsad.wspd).tolist(), np.array(dsapd.prmm).tolist()]
    # INSERT RS DATA HERE
    ####################
    # Create E & T from RS
    # Get the fc files in numpy array converted to [0.0-1.0], fc can be interpolated
    fc = libgetgeo.getgeoydoy(rsdir, "MOD13Q1*_fc.tif", longitude, latitude, year.tolist(), doy, i=True)
    fc = np.asarray(fc, dtype=np.float)
    fc = np.divide(fc, 100)
    # print("START fc", fc, "END fc")
    # Get eta from both MOD and MYD platforms in arrays
    modet = libgetgeo.getgeoydoy(rsdir, "MOD16A*.vrt", longitude, latitude, year.tolist(), doy)
    mydet = libgetgeo.getgeoydoy(rsdir, "MYD16A*.vrt", longitude, latitude, year.tolist(), doy)
    # make 32765 -> NAN in modet/mydet
    modet = np.asarray(modet, dtype=np.float)
    mydet = np.asarray(mydet, dtype=np.float)
    modet[modet == 32765] = np.nan
    modet[modet == 32766] = np.nan
    modet[modet == 32767] = np.nan
    mydet[mydet == 32765] = np.nan
    mydet[mydet == 32766] = np.nan
    mydet[mydet == 32767] = np.nan
    modet = np.divide(modet, 10)
    mydet = np.divide(mydet, 10)
    # print("START modet", modet, "END modet")
    # print("START mydet", mydet, "END mydet")
    # Compute average eta. Take NANs into account to maximise data count out
    et = (modet + mydet) / 2.0
    # print("START et", et, "END et")
    # Transpiration
    t = np.multiply(et, fc)
    # print("START T", T, "END T")
    # Evaporation is the complementary of Transpiration vis-a-vis et
    e = np.subtract(et, t)
    # print("START E", E, "END E")
    # Add Evaporation data to the Meteo file
    meteolist.append(e.tolist())
    # Add Transpiration data to the Meteo file
    meteolist.append(t.tolist())
    # Add lai data to the Meteo file, LAI can be interpolated
    lai = libgetgeo.getgeoydoy(rsdir, "MCD15A2H*.vrt", longitude, latitude, year.tolist(), doy, i=True)
    lai = np.array(lai, dtype=np.float)
    lai = np.divide(lai, 10)
    # print("START lai", lai, "END lai")
    meteolist.append(lai.tolist())
    # TODO Write a function to get the data from imagery of temporal coherence
    s1temporalcoherence = None
    if s1temporalcoherence:
        # Skip the manual trick
        # Create artificial Cut data from a previous array
        cutnans = np.zeros(lai.shape[0])
        # Add Cut data to the Meteo file
        meteolist.append(cutnans.tolist())
    else:
        # Create artificial Cut data from a previous array
        cutnans = np.zeros(lai.shape[0])
        # Fill with NANs
        cutnans.fill(np.nan)
        # Fill with 2 cuts
        cutnans[127] = 1
        cutnans[235] = 1
        # Add Cut data to the Meteo file
        meteolist.append(cutnans.tolist())
    # print("START CUT NANS", cutnans, "END CUT NANS")
    # use soil moisture from the BEC-SMOS 1km dataset
    valuelist = libgetgeobecsmos.getgeoydoy(becsmosdir, "BEC*.nc", longitude, latitude, year.tolist(), doy)
    #print(valuelist)
    meteolist.append(valuelist)
    return meteolist
