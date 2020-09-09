#* version adapted for rye grass april 4, 2006, joost wolf for
#* fst modelling;
#* model is based on lingra model in cgms
#* for forage growth and production simulation which was written fortran
#* (grsim.pfo) and later in c; application of model is the
#* simulation of perennial ryegrass (l. perenne) growth under both
#* potential and water-limited growth conditions.
#*
#* model is different from lingra model in cgms with respect to:
#*  1) evaporation, transpiration, water balance, root depth growth
#*  and growth reduction by drought stress (tranrf) are derived from lingra
#*  model for thimothee (e.g. subroutine penman, evaptr and drunir)
#*  2) running average to calculate soil temperature (soitmp) is derived from
#*  approach in lingra model for thimothee
#*
#*
#*     name    type     description                               unit
#* variables:
#*     biomass in leaves, reserves, etc. in kg dm / ha
#*     terms of water balance in mm/day
#*

import numpy as np
import pandas as pd
from ast import literal_eval

#Parse input values
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("latitude", help="Latitude in d.ddd")
parser.add_argument("meteofile", help="Name of the meteo file")
args = parser.parse_args()

lat  = float(args.latitude) # latitude (used to calculate daylength;
                            # daylength not taken into account in 
                            # lingra, latitude does not affect yield)

df = pd.read_csv(args.meteofile,sep='\t')
year = np.asarray((df[df.columns[0]]), dtype=np.int)    # year in weather file
doy  = np.asarray((df[df.columns[1]]), dtype=np.int)    # doy of the year
rdd  = np.asarray((df[df.columns[2]]), dtype=np.float)  # solar radiation (kj m-2 day-1)
tmmn = np.asarray((df[df.columns[3]]), dtype=np.float)  # minimum temperature (degrees celsius)
tmmx = np.asarray((df[df.columns[4]]), dtype=np.float)  # maximum temperature (degrees celsius)
vp   = np.asarray((df[df.columns[5]]), dtype=np.float)  # water vapour pressure (kpa)
wn   = np.asarray((df[df.columns[6]]), dtype=np.float)  # average wind speed (m s-1)
rain = np.asarray((df[df.columns[7]]), dtype=np.float)  # daily rainfall (mm day-1)
RSevp= np.asarray((df[df.columns[8]]), dtype=np.float)  # daily RS Evaporation actual (mm day-1)
RStrn= np.asarray((df[df.columns[9]]), dtype=np.float)  # daily RS Transpiration actual (mm day-1)
RSlai= np.asarray((df[df.columns[10]]),dtype=np.float)  # daily RS LAI (-)
RScut= np.asarray((df[df.columns[11]]),dtype=np.float)  # daily RS cutting event (0/1)

print(year[0],doy[0],rdd[0],tmmn[0],tmmx[0],vp[0],wn[0],rain[0],RSevp[0],RStrn[0],RSlai[0],RScut[0])


# number of days in simulation (days)
fintim = len(year)

#initial
rootdi  = 0.4
LAIi    = 0.1  

tilli   = 7000.
wrei    = 200.
wrti    = 4.

#************************************************************************
#***   functions and parameters for grass
#************************************************************************
  
#* parameters
co2a   = 360.     # atmospheric co2 concentration (ppm)
kdif   = 0.60     # 
LAIcr  = 4.       # critical LAI 
tbase  = 0.        
luemax = 3.0      # light use efficiency (g dm mj-1 par intercepted)
sla    = 0.0025   # specific leaf area (m2 g-1 dm)
cLAI   = 0.8      # LAI after cutting (-)
nitmax = 3.34     # maximum nitrogen content (%)
nitr   = 3.34     # actual nitrogen content (%)
rdrd   = 0.01     # base death rate (fraction)
tmbas1 = 3.       # base temperature perennial ryegrass (degrees celsius)

#* parameters for water relations from lingra for thimothee
drate  = 50.      # drainage rate (mm day-1)
irrigf = 0.       # irrigation (0 = no irrigation till 1 = full irrigation)
rootdm = 0.4      # maximum root depth (m)
rrdmax = 0.012    # maximum root growth (m day-1)
wcad   = 0.005    # air dry water content (fraction) 
wcwp   = 0.12     # wilting point water content (fraction) 
wcfc   = 0.29     # field capacity water content (fraction)
wci    = 0.29     # initial water content (fraction)
wcwet  = 0.37     # minimum water content at water logging (fraction) 
wcst   = 0.41     # saturation water content (fraction)

pi     = 3.1415927
rad    = pi / 180.


#* harvest dates
mndat0 = RScut
#mndat = [135.,165.,200.,240.,280.]  # harvest days in the simulation 
#mndat0 = np.zeros(fintim)
#mndat0[int(mndat[0])] = 1
#mndat0[int(mndat[1])] = 1
#mndat0[int(mndat[2])] = 1
#mndat0[int(mndat[3])] = 1 
#mndat0[int(mndat[4])] = 1

#      initial available water (mm)
wai     = 1000. * rootdi * wci
#*     initial leaf weight is initialized as initial
#*     leaf area divided by initial specific leaf area, kg ha-1
wlvgi   = LAIi / sla
#*     remaining leaf weight after cutting is initialized at remaining
#*     leaf area after cutting divided by initial specific leaf area, kg ha-1
cwlvg = cLAI/sla
#*      maximum site filling new buds (fsmax) decreases due
#*      to low nitrogen contents, van loo and schapendonk (1992)
#*      theoretical maximum tillering size = 0.693
fsmax = nitr/nitmax*0.693

#************************************************************************
#***   11. data
#************************************************************************

# specifying variables:
davtmp       = np.mean(np.array([tmmn,tmmx]), axis=0)
photmp       = np.divide(np.add(tmmn,3.*tmmx),4.0)
dec          = np.zeros(fintim)
decc         = np.zeros(fintim)
dayl         = np.zeros(fintim)
dtr          = rdd / 1.e+3
parav        = np.zeros(fintim)
rsoitm       = np.zeros(fintim)
soitmp       = np.zeros(fintim+1)
soitmp[0]    = 5. 
efftmp       = np.maximum(davtmp, tbase)
incut        = np.zeros(fintim+1)

# rate variables
redtmp       = np.zeros(fintim)  # reduction in light use efficiency due to temperature
redrdd       = np.zeros(fintim)  # reduction in light use efficiency due to solar radiation
tmeff        = np.zeros(fintim)  # effective temperature (daily temperature above 0 degrees celsius)   
par          = np.zeros(fintim)  # photosynthetic active radiation (mj par m-2 day-1)
fint         = np.zeros(fintim)  # fraction light intercepted
lue1         = np.zeros(fintim)  # light use efficiency (g dm mj-1 par intercepted)
parint       = np.zeros(fintim)  # par intercepted (mj par m-2 day-1)
frt          = np.zeros(fintim)  # fraction assimilates to roots
flv          = np.zeros(fintim)  # fraction assimulates to leaves

# state variables
parcu        = np.zeros(fintim+1)   # cumulative intercepted par, mj par intercepted m-2 ha-1
tracu        = np.zeros(fintim+1)   # cumulative transpiration, mm
tramcu       = np.zeros(fintim+1)   # cumulative maximal transpiration, mm
evacu        = np.zeros(fintim+1)   # cumulative evaporation, mm
evamcu       = np.zeros(fintim+1)   # cumulative maximum evaporation, mm
irrcu        = np.zeros(fintim+1)   # cumulative irrigation, mm
tsum         = np.zeros(fintim+1)   # sum of temperatures above base temperature, gr. c.d
dvs          = np.zeros(fintim)     # hypothetical development stage, 600 gr. c.d taken from subroutine tilsub
LAI          = np.zeros(fintim+1)
LAI[0]       = LAIi                 # leaf area index, ha ha-1
daha         = np.zeros(fintim+1)   # days after harv, d
tiller       = np.zeros(fintim+1)   
tiller[0]    = tilli                # number of tillers, tillers m-2
wlvg         = np.zeros(fintim+1)   
wlvg[0]      = wlvgi                # dry weight of green leaves, kg ha-1
wlvd         = np.zeros(fintim+1)   # dry weight of dead leaves, kg ha-1 incl. harvests
wlvd1        = np.zeros(fintim)     # dry weight of dead leaves, kg ha-1
hrvbl        = np.zeros(fintim+1) 
hrvbl[0]     = wlvg[0] - cwlvg      # harvestable leaf weight
grass        = np.zeros(fintim+1)   # dry weight of cutted green leaves, kg ha-1
wre          = np.zeros(fintim+1)   
wre[0]       = wrei                 # dry weight of storage carbohydrates, kg ha-1
wrt          = np.zeros(fintim+1) 
wrt[0]       = wrti                 # dry weight of roots, kg ha-1
tadrw        = np.zeros(fintim)     # total above ground dry weight including harvests, kg ha-1
yielD        = np.zeros(fintim+1)   # harvestable part of total above ground dry weight and previous harvests, kg ha-1
length       = np.zeros(fintim+1)   # length of leaves, cm
sLAInt       = np.zeros(fintim+1)   # running specific leaf area in model, ha kg-1
sLAInt[0]    = LAI[0] / wlvg[0]
rootd        = np.zeros(fintim+1)  
rootd[0]     = rootdi               # rooting depth (from lingra for timothee)
wa           = np.zeros(fintim+1)   
wa[0]        = wai                  # soil water in rooted zone  (from lingra for timothee)
raincu       = np.zeros(fintim+1)   # cumulative rainfall
dracu        = np.zeros(fintim+1)   # cumulative drainage (avdl)
runcu        = np.zeros(fintim+1)     # cumulative runoff (avdl)
intLAI       = np.zeros(fintim+1)     # cumulative rainfall intercepted by leaves (avdl)

#* ---------------------------------------------------------------------*
#*     subroutine penman                                                *
#*     purpose: computation of the penman equation                      *
#* ---------------------------------------------------------------------*

dtrjm2       = rdd * 1.e3
boltzm       = 5.668e-8
lhvap        = 2.4e6
psych        = 0.067

bbrad        = boltzm * (davtmp+273.)**4 * 86400.
svp          = 0.611 * np.exp(17.4 * davtmp / (davtmp + 239.))  
slope        = 4158.6 * svp / np.power((davtmp + 239.),2)

rlwn         = np.zeros(fintim)
nrads        = np.zeros(fintim)
nradc        = np.zeros(fintim)

penmrs       = np.zeros(fintim)
penmrc       = np.zeros(fintim)

wdf          = np.zeros(fintim)
penmd        = np.zeros(fintim)

pevap        = np.zeros(fintim)
ptran        = np.zeros(fintim)
rnintc       = np.zeros(fintim)


#* ---------------------------------------------------------------------*
#*     subroutine evaptr                                                *
#*     purpose: to compute actual rates of evaporation and transpiration*
#* ---------------------------------------------------------------------*

wc           = np.zeros(fintim)
waad         = np.zeros(fintim)
wafc         = np.zeros(fintim)
evap1        = np.zeros(fintim)
wccr = wcwp + 0.5 * (wcfc-wcwp)
fr           = np.zeros(fintim)
tran         = np.zeros(fintim)
availf       = np.zeros(fintim)
evap         = np.zeros(fintim)
tran         = np.zeros(fintim)

#* ---------------------------------------------------------------------*
#*     subroutine drunir                                                *
#*     purpose: to compute rates of drainage, runoff and irrigation     *
#* ---------------------------------------------------------------------*

wast         = np.zeros(fintim)
drain        = np.zeros(fintim)
runoff       = np.zeros(fintim)
irrig        = np.zeros(fintim)

#************************************************************************
#***   5. water balance and root depth growth (from lingra for thymothee)
#************************************************************************
rdaha        = np.zeros(fintim)
harv         = np.zeros(fintim)
rrootd       = np.zeros(fintim)
explor       = np.zeros(fintim)
tranrf       = np.zeros(fintim)
rwa          = np.zeros(fintim)
leafn        = np.zeros(fintim)
lera         = np.zeros(fintim)
lera2        = np.zeros(fintim)

# subroutine tilsub
dtil         = np.zeros(fintim)
reftil       = np.zeros(fintim)
dtild        = np.zeros(fintim)

# subroutine sosub
lued         = np.zeros(fintim)
gtwso        = np.zeros(fintim)

gtwso2       = np.zeros(fintim)
dLAIs        = np.zeros(fintim)
dre          = np.zeros(fintim)
gtwsi        = np.zeros(fintim)
gre          = np.zeros(fintim)
gtw          = np.zeros(fintim)
rre          = np.zeros(fintim)

rdrsh        = np.zeros(fintim)
rdrsm        = np.zeros(fintim)
rdrs         = np.zeros(fintim)
rdr          = np.zeros(fintim)
grt          = np.zeros(fintim)
gLAI         = np.zeros(fintim)
dLAI         = np.zeros(fintim)
rLAI         = np.zeros(fintim)
dlv          = np.zeros(fintim)
glv          = np.zeros(fintim)
rlv          = np.zeros(fintim)

newbio       = np.zeros(fintim)
rratio       = np.zeros(fintim)
lueycu       = np.zeros(fintim)
wrtmin       = np.zeros(fintim)

time = list(range(0,fintim,1))

#dynamic
for i in range(fintim):
  vp[i] = min(vp[i], svp[i])
  rlwn[i]         = bbrad[i] * max(0.,0.55*(1.-vp[i]/svp[i]))
  nrads[i]        = dtrjm2[i] * (1.-0.15) - rlwn[i]
  nradc[i]        = dtrjm2[i] * (1.-0.25) - rlwn[i]
  penmrs[i]       = nrads[i] * slope[i]/(slope[i]+psych)
  penmrc[i]       = nradc[i] * slope[i]/(slope[i]+psych)
  wdf[i]          = 2.63 * (1.0 + 0.54 * wn[i])
  penmd[i]        = lhvap * wdf[i] * (svp[i]-vp[i]) * psych/(slope[i]+psych)
  
  # northern / southern hemisphere! warnings from this line!
  dec[i]    = -np.arcsin(np.sin(23.45*rad)*np.cos(2.*pi*(time[i]+10.)/365.))

  if(np.arctan(-1./np.tan(rad*lat)) > dec[i]):
      decc[i] = np.arctan(-1./np.tan(rad*lat))
  elif(np.arctan(1./np.tan(rad*lat)) < dec[i]):
      decc[i] = np.arctan(-1./np.tan(rad*lat)) 
  else: 
      decc[i] = dec[i] 
  
  dayl[i]   = 0.5 * ( 1. + 2. * np.arcsin(np.tan(rad*lat)*np.tan(decc[i])) / pi )
  
  parav[i]  = 0.5 * dtr[i] / dayl[i]
  #soil temperature changes
  rsoitm[i] = (davtmp[i]-soitmp[i]) / 10.
  soitmp[i+1] = soitmp[i] + rsoitm[i]
  tmeff[i]    = max(davtmp[i] - tmbas1, 0.)
  
  # rate variables
  if(soitmp[i] <3):
      redtmp[i] =0 
  elif(soitmp[i]>8):
      redtmp[i] = 1 
  else:
      redtmp[i] = (soitmp[i]-3)*0.2 # luerd1

  if(dtr[i] <10):   
      redrdd[i] =1 
  else: 
      redrdd[i] = ((1+(10/(40-10)*0.67)) - dtr[i] * (0.67/(40-10)))  # luerd2
  
  # daily photosynthetically active radiation, mj m-2 d-1
  par[i] = dtr[i] * 0.50
  # fraction of light interception
  fint[i] = (1. - np.exp(-kdif*LAI[i]))
  # light use efficiency, g mj par-1
  lue1[i] = luemax * redtmp[i] * redrdd[i]
  # total intercepted photosynthetically active radiation, mj m-2 d-1       
  parint[i] = fint[i] * par[i]
  # fraction of dry matter allocated to roots, kg kg-1
  
  #* ---------------------------------------------------------------------*
  #*     subroutine penman                                                *
  #*     purpose: computation of the penman equation                      *
  #* ---------------------------------------------------------------------*
  pevap[i]  = max(0.0, np.exp(-0.5*LAI[i])  * (penmrs[i] + penmd[i]) / lhvap)
  ptran[i]  = (1.-np.exp(-0.5*LAI[i])) * (penmrc[i] + penmd[i]) / lhvap
  rnintc[i] = min(rain[i], 0.25*LAI[i])
  ptran[i]  = max(0.0, ptran[i]-0.5*rnintc[i])
  
  #* ---------------------------------------------------------------------*
  #*     subroutine evaptr                                                *
  #*     purpose: to compute actual rates of evaporation and transpiration*
  #* ---------------------------------------------------------------------*
  wc[i]   = 0.001 * wa[i]   / rootd[i]
  waad[i] = 1000. * wcad * rootd[i]
  wafc[i] = 1000. * wcfc * rootd[i]
  evap1[i]  = pevap[i] * max(0., min(1.,(wc[i]-wcad)/(wcfc-wcad)))
  if(wc[i]>wccr):
      fr[i] = max(0.,min(1.,(wcst-wc[i])/(wcst-wcwet))) 
  else:
      fr[i] = max(0., min(1.,(wc[i]-wcwp)/(wccr-wcwp))) 

  ####### INSERT_RStrn_HERE ##########
  if(np.isnan(RStrn[i])):
    tran[i] = ptran[i] * fr[i]
  else:
    tran[i] = RStrn[i]
  ####################################

  if(evap1[i]+tran[i] == 0.0):
      availf[i] = min(1., (wa[i]-waad[i])/0.01) # notnul removed!
  else:
      availf[i] = min(1., (wa[i]-waad[i])/(evap1[i]+tran[i])) # notnul removed!
  ####### INSERT_RSevp_HERE ##########
  if(np.isnan(RSevp[i])):
    evap[i] = evap1[i] * availf[i]
  else:
    evap[i] = RSevp[i]
  ####################################

  tran[i] = tran[i] * availf[i]
  
  #* ---------------------------------------------------------------------*
  #*     subroutine drunir                                                *
  #*     purpose: to compute rates of drainage, runoff and irrigation     *
  #* ---------------------------------------------------------------------*
  wast[i] = 1000. * wcst * rootd[i]
  drain[i]  = max(0, min(drate,(wa[i]-wafc[i] + rain[i] - rnintc[i] - evap[i] - tran[i])))   
  runoff[i] = max(0, wa[i]-wast[i] + rain[i] - rnintc[i] - evap[i] - tran[i] - drain[i])
  irrig[i]  = irrigf * (wafc[i]-wa[i] - (rain[i] - rnintc[i] - evap[i] - tran[i] - drain[i] - runoff[i]))
                                            
  #************************************************************************
  #***   5. water balance and root depth growth (from lingra for thymothee)
  #************************************************************************
  if(rootdm-rootd[i] <= 0 or wc[i]-wcwp <=0):
      rrootd[i] = rrdmax * 0 
  else:
      rrootd[i] = rrdmax * 1  
  
  explor[i] = 1000. * rrootd[i] * wcfc # assumption that explored layers are at fc!
  
  if(ptran[i] <= 0):
      tranrf[i] = 1 
  else: 
      tranrf[i] = min(1, tran[i] / ptran[i]) # insw removed!
                        
  rwa[i] = (rain[i]+explor[i]+irrig[i]) - (rnintc[i]+runoff[i]+tran[i]+evap[i]+drain[i])
  frt[i] = min(0.263, 0.263 - tranrf[i]*(0.263-0.165))        # frrttb 
  flv[i] = 1.-frt[i]
  
  # subroutine mowing 
  rdaha[i] = 1  
  harv[i] = 0
  incut[i+1] = incut[i]
  
  #     mowing at observation dates
  #     reset days after harv
  if(RScut[i] == 1 and wlvg[i] > cwlvg):
      harv[i]  = wlvg[i]-cwlvg
  if(RScut[i] == 1 and wlvg[i] > cwlvg):
      rdaha[i] = -daha[i]
  if(RScut[i] == 1 and wlvg[i] > cwlvg):
      incut[i+1] = incut[i] + 1.
  
  #     no mowing in current season, do not increase rate
  #     of days after harv
   
  if(incut[i] == 0.):
      rdaha[i] = 0
  
  #     mowing in current season, increase rate of days
  #     after harvests
  
  # temperature dependent leaf appearance rate, according to
  # (davies and thomas, 1983), soil temperature (soitmp)is used as
  # driving force which is estimated from a 10 day running average
         
  if(redtmp[i] > 0):
      leafn[i] = soitmp[i]*0.01 
  else:
      leafn[i] = 0
  
  # leaf elongation rate affected by temperature
  # cm day-1 tiller-1
  
  if(davtmp[i]-tmbas1 > 0):
      lera[i] = 0.83*np.log(max(davtmp[i], 2.))-0.8924 
  else:
      lera[i] = 0   # log10 or log?
  
  if((harv[i]-0.1)<0):
      lera2[i] = lera[i] 
  else:
      lera2[i] = -1*length[i] 
  
  # subroutine tilsub
  dtil[i] = 0.
  if(daha[i] < 8):
      reftil[i] = max(0.,0.335-0.067*LAI[i]) * redtmp[i] 
  else:
      reftil[i] = max(0., min(fsmax, 0.867-0.183*LAI[i])) * redtmp[i]
  
  # relative rate of tiller formation when defoliation less
  # than 8 days ago, tiller tiller-1 d-1
  # relative rate of tiller formation when defoliation is more
  # than 8 days ago, tiller tiller-1 d-1
  # relative death rate of tillers due to self-shading (dtild),
  # tiller tiller-1 d-1
  
  dtild[i] = max(0.01*(1.+tsum[i]/600.), 0.05 * (LAI[i]-LAIcr)/LAIcr)
  
  if(tiller[i] < 14000):
      dtil[i] = (reftil[i]-dtild[i]) * leafn[i] * tiller[i] 
  else:
      dtil[i] = -dtild[i] * leafn[i] * tiller[i] 
  
  # subroutine sosub 
  lued[i] = min(lue1[i]*(0.336+0.224*nitr)/(0.336+0.224*nitmax), lue1[i]*tranrf[i])
   
  # start of growing season
  if(harv[i] == 0):
      gtwso[i] = lued[i] * parint[i] * (1.+0.8*np.log(co2a/360.)) * 10. 
  else:
      gtwso[i] = 0. 
  
  # rate of sink limited leaf growth, unit of tiller is tillers m-2 (!),
  # 1.0e-8 is conversion from cm-2 to ha-1, ha leaf ha ground-1 d-1
  dLAIs[i] = (tiller[i] * 1.0e4 * (lera[i] * 0.3)) * 1.0e-8
  gtwso2[i] = gtwso[i]+wre[i]
  dre[i]   = wre[i]
  
  # conversion to total sink limited carbon demand,
  # kg leaf ha ground-1 d-1
  if(harv[i] <= 0):
      gtwsi[i] = dLAIs[i] * (1./sla) * (1./flv[i]) 
  else:
      gtwsi[i] = 0 
  
  # actual growth switches between sink- and source limitation
  # (more or less dry matter formed than can be stored)
  if((gtwso2[i] - gtwsi[i]) > 0):
      gre[i] = gtwso2[i] - gtwsi[i] 
  else:
      gre[i] = 0
  if((gtwso2[i] - gtwsi[i]) <= 0):
      gtw[i] = gtwso2[i] 
  else:
      gtw[i] = gtwsi[i] 
  
  # change in reserves
  rre[i] = gre[i]-dre[i]
  # relative death rate of leaves due to self-shading, d-1
  rdrsh[i] = max(0., min(0.03, 0.03 * (LAI[i]-LAIcr) /LAIcr))
  
  # relative death rate of leaves due to drought stress, d-1
  rdrsm[i] = max(0., min(0.05, 0.05 * (1.-tranrf[i])))
  
  # maximum of relative death rate of leaves due to
  # and drought stres, d-1
  rdrs[i] = max(rdrsh[i], rdrsm[i])
  
  # actual relative death rate of leaves is sum of base death
  # rate plus maximum of death rates rdrsm and rdrsh, d-1
  rdr[i] = rdrd + rdrs[i]
  
  # actual growth rate of roots, kg ha-1 d-1
  grt[i] = gtw[i] * frt[i]
  
  # actual growth rate of leaf area, ha ha-1 d-1
  gLAI[i] = gtw[i] * flv[i] * sla
  
  # actual death rate of leaf area, due to relative death
  # rate of leaf area or rate of change due to cutting, ha ha-1 d-1
  if(harv[i] <= 0):
      dLAI[i] = LAI[i] * (1. - np.exp(-rdr[i])) 
  else:
      dLAI[i] = harv[i] * sLAInt[i]  
  
  # change in LAI
  rLAI[i]= gLAI[i]-dLAI[i]
  # actual death rate of leaves, kg ha-1 d-1 incl. harvested leaves
  dlv[i] = dLAI[i] / sLAInt[i]    # notnul removed!
  
  # rate of change of dry weight of green leaves due to
  # growth and senescence of leaves or periodical harvest, kg ha-1 d-1
  if(harv[i] <= 0):
      glv[i] = gtw[i]*flv[i] 
  else:
      glv[i] = 0  
  
  # change in green leaf weight
  rlv[i] = glv[i]-dlv[i]
  
  ####### INSERT_RSlai_HERE ##########
  if((i+1) < len(RSlai)):
    if(np.isnan(RSlai[i+1])):
        LAI[i+1]    = LAI[i]    + rLAI[i]             #leaf area index, ha ha-1
    else:
        LAI[i+1]    = RSlai[i+1]

  ####### INSERT_RStrn_HERE ##########
  if((i+1) < len(RStrn)):
    if(np.isnan(RStrn[i+1])):
        tracu[i+1]    = tracu[i] + tran[i]             #cumulative transpiration, mm
    else:
        tracu[i+1]    = tracu[i] + RStrn[i+1]

  ####### INSERT_RSevp_HERE ##########
  if((i+1) < len(RSevp)):
    if(np.isnan(RSevp[i+1])):
        evacu[i+1]      = evacu[i] + evap[i]             #cumulative evaporation, mm
    else:
        evacu[i+1]      = evacu[i] + RSevp[i]

  # state variables
  parcu[i+1]      = parcu[i]  + parint[i]           #cumulative intercepted par, mj par intercepted m-2 ha-1
  tramcu[i+1]     = tramcu[i] + ptran[i]            #cumulative maximal transpiration, mm
  evamcu[i+1]     = evamcu[i] + pevap[i]            #cumulative maximal evaporation, mm
  irrcu[i+1]      = irrcu[i]  + irrig[i]            #cumulative irrigation, mm
  tsum[i+1]       = tsum[i]   + tmeff[i]            #sum of temperatures above base temperature, gr. c.d
  dvs[i]          = tsum[i] / 600.                  #hypot development stage, 600 gr. c.d taken fr subr tilsub
  daha[i+1]       = daha[i]   + rdaha[i]            #days after harv, d
  tiller[i+1]     = tiller[i] + dtil[i]             #number of tillers, tillers m-2
  wlvg[i+1]       = wlvg[i]   + rlv[i]              #dry weight of green leaves, kg ha-1
  wlvd[i+1]       = wlvd[i]   + dlv[i]              #dry weight of dead leaves, kg ha-1 incl. harvests
  wlvd1[i]        = wlvd[i]   - grass[i]            #dry weight of dead leaves, kg ha-1
  hrvbl[i+1]      = wlvg[i+1] - cwlvg               #harvestable leaf weight
  grass[i+1]      = grass[i]  + harv[i]             #dry weight of cutted green leaves, kg ha-1
  wre[i+1]        = wre[i]    + rre[i]              #dry weight of storage carbohydrates, kg ha-1
  wrt[i+1]        = wrt[i]    + grt[i]              #dry weight of roots, kg ha-1
  tadrw[i]        = grass[i]  + wlvg[i]             #total above ground dry weight including harvests,kg ha-1
  yielD[i+1]      = grass[i+1]+ max(0., hrvbl[i+1]) #harvestable part of total above ground dry weight and previous harvests, kg ha-1
  length[i+1]     = length[i] + lera2[i]            #length of leaves, cm
  sLAInt[i+1]     = LAI[i+1]  / wlvg[i+1]           #running specific leaf area in model,hakg-1 notnul rmved!
  rootd[i+1]      = rootd[i]  + rrootd[i]           #rooting depth (from lingra for timothee)  
  wa[i+1]         = wa[i]     + rwa[i]              #soil water in rooted zone  (from lingra for timothee)
  raincu[i+1]     = raincu[i] + rain[i]             #cumulative rainfall  
  dracu[i+1]      = dracu[i]  + drain[i]            #cumulative drainage (avdl)
  runcu[i+1]      = runcu[i]  + runoff[i]           #cumulative runoff (avdl)
  intLAI[i+1]     = intLAI[i] + rnintc[i]           #cumulative rainfall intercepted by leaves (avdl)

  #************************************************************************
  #  ***   9. additional variables and parameters for output
  #************************************************************************
  
  #newbio[i] = wlvg[i]+wrt[i]+wre[i] - (wlvgi+wrti+wrei) + grass[i]
  #if(newbio[i] == 0.):
  #    rratio[i] = 0 
  #else:
  #    rratio[i] = max(0., min(1., (wrt[i]-wrti) / newbio[i])) 
  #lueycu[i] = yielD[i]  / parcu[i]
  #wrtmin[i] = -wrt[i]
  
#ENDLOOP}
print("tiller[365], yielD[365], wlvg[365], wlvd1[365], parcu[365], grass[365], tracu[365], evacu[365]")
output = [tiller[-1], yielD[-1], wlvg[-1], wlvd1[-1], parcu[-1], grass[-1], tracu[-1], evacu[-1]]
print(output)

#PLOT
import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(15,7))

plt.subplot(231)
plt.plot(time, np.divide(yielD[:len(time)],1000),'g-')
plt.title('Yield accumulated (T/ha)')
plt.grid()

plt.subplot(232)
plt.plot(time, tracu[:len(time)],'b-')
plt.title('Transpiration accumulated (mm)')
plt.grid()

plt.subplot(236)
plt.plot(time, np.divide(tiller[:len(time)],100),'r-')
plt.title('Tiller (M#/ha)')
plt.grid()

plt.subplot(234)
plt.plot(time, LAI[:len(time)],'g-')
plt.title('LAI (-)')
plt.grid()

plt.subplot(235)
plt.plot(time, evacu[:len(time)],'b-')
plt.title('Evaporation accumulated (mm)')
plt.grid()

plt.subplot(233)
plt.plot(time, np.divide(grass[:len(time)],1000),'r-')
plt.title('Grass accumulated (T/ha)')
plt.grid()

plt.tight_layout()
plt.show()

# data processing
#waterinput = raincu + irrcu  
#wateroutput = tracu + evacu + intLAI + dracu + runcu    

# output in table and .csv file
#output = [year, davtmp, dtr, dayl, LAI, tiller, yielD, wlvg, wlvd1, wre, wrt, wrtmin, 
#                rratio, sla, sLAInt, parcu, lueycu, grass, tadrw, vp, newbio, tranrf, 
#                tran, ptran, evap, pevap, wc, wa, rain, raincu, wn,tramcu, tracu, evamcu,
#                evacu, irrcu, gtw, gtwsi, length]
#write.csv(output, file= "m:/r/output/lingra.csv") # write output to a .csv file

  
