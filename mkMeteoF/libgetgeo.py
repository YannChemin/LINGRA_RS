import os
import glob

def getgeo(fName, longitude, latitude):
    return(os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (fName, longitude, latitude)).read())

def getdoy(fName):
    return(fName[-7:-4])

def getyear(fName):
    return(fName[-11:-7])

def getgeoWLDCRD(VRTdir, wildcard, longitude, latitude):
    py_glob = os.path.join(VRTdir, wildcard)
    result=[]
    for path in glob.iglob(py_glob):
        result.append(getyear(path),getdoy(path),getgeo(path,longitude,latitude))

def getgeoydoy(VRTdir, wildcard, longitude, latitude, listyear, listdoy):
    py_glob = os.path.join(VRTdir, wildcard)
    RSyear=[]
    RSdoy=[]
    value=[]
    for path in glob.iglob(py_glob):
        year=getyear(path)
        doy=getdoy(path)
        RSyear.append(year)
        RSdoy.append(doy)
        value.append(year,doy,getgeo(path,longitude,latitude))
        
    result=[]
    #Fill list with nan for the len(listdoy)
    for d in range(len(listdoy)):
        for l in range(len(RSdoy)):
            if(listyear[d]==RSyear[l]):
                if(listdoy[d]==RSdoy[l]):
                    result.append(value[l])
            else:
                result.append('nan')



