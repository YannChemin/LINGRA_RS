// version adapted for rye grass april 4, 2006, joost wolf for fst modelling;
// model is based on lingra model in cgms
// for forage growth and production simulation which was written fortran
// (grsim.pfo) and later in c; application of model is the
// simulation of perennial ryegrass (l. perenne) growth under both
// potential and water-limited growth conditions.
// model is different from lingra model in cgms with respect to:
//  1) evaporation, transpiration, water balance, root depth growth
//  and growth reduction by drought stress (tranrf) are derived from lingra
//  model for thimothee (e.g. subroutine penman, evaptr and drunir)
//  2) running average to calculate soil temperature (soitmp) is derived from
//  approach in lingra model for thimothee
//  variables:
//     biomass in leaves, reserves, etc. in kg dm / ha
//     terms of water balance in mm/day

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "liblingraRSconst.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

int lingrars(double *results,
            double latitude,
            int lengthsimdays,
            int *year,
            int *doy,
            double *rdd,
            double *tmmn,
            double *tmmx,
            double *vp,
            double *wn,
            double *rain,
            double *RSevp,
            double *RStrn,
            double *RSlai,
            double *RScut,
            double *RSsms){
    //"""
    //Function for Lingra_RS
    //:param latitude: required for daylength calculation
    //:param meteolist: input for weather and RS params
    //:param plot: False (process mode) or True (debug only)
    //:return:[tiller[-1], yielD[-1], wlvg[-1], wlvd1[-1], wa[-1], grass[-1], tracu[-1], evacu[-1]]
    //"""
    //import numpy as np
    //year = meteolist[0]  // year in weather file
    //doy = meteolist[1]  // doy of the year
    //rdd = meteolist[2]  // solar radiation (kj m-2 day-1)
    //tmmn = meteolist[3]  // minimum temperature (degrees celsius)
    //tmmx = meteolist[4]  // maximum temperature (degrees celsius)
    //vp = meteolist[5]  // water vapour pressure (kpa)
    //wn = meteolist[6]  // average wind speed (m s-1)
    //rain = meteolist[7]  // daily rainfall (mm day-1)
    //RSevp = meteolist[8]  // daily RS Evaporation actual (mm day-1)
    //RStrn = meteolist[9]  // daily RS Transpiration actual (mm day-1)
    //RSlai = meteolist[10]  // daily RS LAI (-)
    //RScut = meteolist[11]  // daily RS cutting event (0/1)
    //RSsms = meteolist[11]  // daily RS soil moisture (cm3/cm3)
    // print(year[0], doy[0], rdd[0], tmmn[0], tmmx[0], vp[0], wn[0], rain[0], RSevp[0], RStrn[0], RSlai[0], RScut[0])
    // print(year[1], doy[1], rdd[1], tmmn[1], tmmx[1], vp[1], wn[1], rain[1], RSevp[1], RStrn[1], RSlai[1], RScut[1])
    // print(year[2], doy[2], rdd[2], tmmn[2], tmmx[2], vp[2], wn[2], rain[2], RSevp[2], RStrn[2], RSlai[2], RScut[2])

    // number of days in simulation (days)
    int fintim = lengthsimdays;

    // specifying variables:
    double davtmp = np.mean(np.array([tmmn, tmmx]), axis=0); //TODO
    double *dec = ad1d(int fintim);
    double *decc = ad1d(int fintim);
    double *dayl = ad1d(int fintim);
    double dtr = rdd / 1000.0;
    double parav = ad1d(int fintim);
    double *rsoitm = ad1d(int fintim);
    double *soitmp = ad1d(int fintim + 1);
    soitmp[0] = 5.0;
    double *incut = ad1d(int fintim + 1);

    // rate variables
    double *redtmp = ad1d(int fintim);  // reduction in light use efficiency due to temperature
    double *redrdd = ad1d(int fintim);  // reduction in light use efficiency due to solar radiation
    double *tmeff = ad1d(int fintim);  // effective temperature (daily temperature above 0 degrees celsius)
    double *par = ad1d(int fintim);  // photosynthetic active radiation (mj par m-2 day-1)
    double *fint = ad1d(int fintim);  // fraction light intercepted
    double *lue1 = ad1d(int fintim);  // light use efficiency (g dm mj-1 par intercepted)
    double *parint = ad1d(int fintim);  // par intercepted (mj par m-2 day-1)
    double *frt = ad1d(int fintim);  // fraction assimilates to roots
    double *flv = ad1d(int fintim);  // fraction assimulates to leaves

    // state variables
    double *parcu = ad1d(int fintim + 1);  // cumulative intercepted par, mj par intercepted m-2 ha-1
    double *tracu = ad1d(int fintim + 1);  // cumulative transpiration, mm
    double *tramcu = ad1d(int fintim + 1);  // cumulative maximal transpiration, mm
    double *evacu = ad1d(int fintim + 1);  // cumulative evaporation, mm
    double *evamcu = ad1d(int fintim + 1);  // cumulative maximum evaporation, mm
    double *irrcu = ad1d(int fintim + 1);  // cumulative irrigation, mm
    double *tsum = ad1d(int fintim + 1);  // sum of temperatures above base temperature, gr. c.d
    double *dvs = ad1d(int fintim);  // hypothetical development stage, 600 gr. c.d taken from subroutine tilsub
    double *LAI = ad1d(int fintim + 1);
    LAI[0] = LAIi  // leaf area index, ha ha-1
    double *daha = ad1d(int fintim + 1);  // days after harv, d
    double *tiller = ad1d(int fintim + 1);
    tiller[0] = tilli  // number of tillers, tillers m-2
    double *wlvg = ad1d(int fintim + 1);
    wlvg[0] = wlvgi  // dry weight of green leaves, kg ha-1
    double *wlvd = ad1d(int fintim + 1);  // dry weight of dead leaves, kg ha-1 incl. harvests
    double *wlvd1 = ad1d(int fintim);  // dry weight of dead leaves, kg ha-1
    double *hrvbl = ad1d(int fintim + 1);
    hrvbl[0] = wlvg[0] - cwlvg  // harvestable leaf weight
    double *grass = ad1d(int fintim + 1);  // dry weight of cutted green leaves, kg ha-1
    double *wre = ad1d(int fintim + 1);
    wre[0] = wrei  // dry weight of storage carbohydrates, kg ha-1
    double *wrt = ad1d(int fintim + 1);
    wrt[0] = wrti  // dry weight of roots, kg ha-1
    double *tadrw = ad1d(int fintim);  // total above ground dry weight including harvests, kg ha-1
    double *yielD = ad1d(int fintim + 1);  // harvestable part of total above ground dry weight and previous harvests, kg ha-1
    double *length = ad1d(int fintim + 1);  // length of leaves, cm
    double *sLAInt = ad1d(int fintim + 1);  // running specific leaf area in model, ha kg-1
    sLAInt[0] = LAI[0] / wlvg[0]
    double *rootd = ad1d(int fintim + 1);
    rootd[0] = rootdi  // rooting depth (from lingra for timothee)
    double *wa = ad1d(int fintim + 1);
    wa[0] = wai  // soil water in rooted zone  (from lingra for timothee)
    double *raincu = ad1d(int fintim + 1);  // cumulative rainfall
    double *dracu = ad1d(int fintim + 1);  // cumulative drainage (avdl)
    double *runcu = ad1d(int fintim + 1);  // cumulative runoff (avdl)
    double *intLAI = ad1d(int fintim + 1);  // cumulative rainfall intercepted by leaves (avdl)

    // ---------------------------------------------------------------------*
    //     subroutine penman                                                *
    //     purpose: computation of the penman equation                      *
    // ---------------------------------------------------------------------*

    double dtrjm2 = rdd * 1000.0;
    double boltzm = pow(5.668,-8);
    double lhvap = pow(2.4,6);
    double psych = 0.067;

    double bbrad = boltzm * pow((davtmp + 273.15),4) * 86400.0;
    double svp = 0.611 * exp(17.4 * davtmp / (davtmp + 239.));
    double slope = 4158.6 * svp / pow((davtmp + 239.), 2);

    double *rlwn = ad1d(int fintim);
    double *nrads = ad1d(int fintim);
    double *nradc = ad1d(int fintim);

    double *penmrs = ad1d(int fintim);
    double *penmrc = ad1d(int fintim);

    double *wdf = ad1d(int fintim);
    double *penmd = ad1d(int fintim);

    double *pevap = ad1d(int fintim);
    double *ptran = ad1d(int fintim);
    double *rnintc = ad1d(int fintim);

    // ---------------------------------------------------------------------*
    //     subroutine evaptr                                                *
    //     purpose: to compute actual rates of evaporation and transpiration*
    // ---------------------------------------------------------------------*

    double *wc = ad1d(int fintim);
    double *waad = ad1d(int fintim);
    double *wafc = ad1d(int fintim);
    double *evap1 = ad1d(int fintim);
    double wccr = wcwp + 0.5 * (wcfc - wcwp);
    double *fr = ad1d(int fintim);
    double *availf = ad1d(int fintim);
    double *evap = ad1d(int fintim);
    double *tran = ad1d(int fintim);

    // ---------------------------------------------------------------------*
    //     subroutine drunir                                                *
    //     purpose: to compute rates of drainage, runoff and irrigation     *
    // ---------------------------------------------------------------------*

    double *wast = ad1d(int fintim);
    double *drain = ad1d(int fintim);
    double *runoff = ad1d(int fintim);
    double *irrig = ad1d(int fintim);

    // 5. water balance and root depth growth (from lingra for thymothee)
    double *rdaha = ad1d(int fintim);
    double *harv = ad1d(int fintim);
    double *rrootd = ad1d(int fintim);
    double *explor = ad1d(int fintim);
    double *tranrf = ad1d(int fintim);
    double *rwa = ad1d(int fintim);
    double *leafn = ad1d(int fintim);
    double *lera = ad1d(int fintim);
    double *lera2 = ad1d(int fintim);

    // subroutine tilsub
    double *dtil = ad1d(int fintim);
    double *reftil = ad1d(int fintim);
    double *dtild = ad1d(int fintim);

    // subroutine sosub
    double *lued = ad1d(int fintim);
    double *gtwso = ad1d(int fintim);

    double *gtwso2 = ad1d(int fintim);
    double *dLAIs = ad1d(int fintim);
    double *dre = ad1d(int fintim);
    double *gtwsi = ad1d(int fintim);
    double *gre = ad1d(int fintim);
    double *gtw = ad1d(int fintim);
    double *rre = ad1d(int fintim);

    double *rdrsh = ad1d(int fintim);
    double *rdrsm = ad1d(int fintim);
    double *rdrs = ad1d(int fintim);
    double *rdr = ad1d(int fintim);
    double *grt = ad1d(int fintim);
    double *gLAI = ad1d(int fintim);
    double *dLAI = ad1d(int fintim);
    double *rLAI = ad1d(int fintim);
    double *dlv = ad1d(int fintim);
    double *glv = ad1d(int fintim);
    double *rlv = ad1d(int fintim);

    // dynamic
    for ( i = 0 ; i < fintim ; i++ ){
        vp[i] = MIN(vp[i], svp[i]);
        rlwn[i] = bbrad[i] * MAX(0., 0.55 * (1. - vp[i] / svp[i]));
        nrads[i] = dtrjm2[i] * (1. - 0.15) - rlwn[i];
        nradc[i] = dtrjm2[i] * (1. - 0.25) - rlwn[i];
        penmrs[i] = nrads[i] * slope[i] / (slope[i] + psych);
        penmrc[i] = nradc[i] * slope[i] / (slope[i] + psych);
        wdf[i] = 2.63 * (1.0 + 0.54 * wn[i]);
        penmd[i] = lhvap * wdf[i] * (svp[i] - vp[i]) * psych / (slope[i] + psych);

        // northern / southern hemisphere! warnings from this line!
        dec[i] = -asin(sin(23.45 * rad) * cos(2. * pi * (doy[i] + 10.) / 365.));

        if (atan2(-1. / tan(rad * latitude)) > dec[i]){
            decc[i] = atan2(-1. / tan(rad * latitude));
        } else if (atan2(1. / tan(rad * latitude)) < dec[i]){
            decc[i] = atan2(-1. / tan(rad * latitude));
        } else {
            decc[i] = dec[i];
        }
        dayl[i] = 0.5 * (1. + 2. * asin(tan(rad * latitude) * tan(decc[i])) / pi);

        parav[i] = 0.5 * dtr[i] / dayl[i];
        // soil temperature changes
        rsoitm[i] = (davtmp[i] - soitmp[i]) / 10.0;
        soitmp[i + 1] = soitmp[i] + rsoitm[i];
        tmeff[i] = max(davtmp[i] - tmbas1, 0.);

        // rate variables
        if (soitmp[i] < 3.0){
            redtmp[i] = 0;
        } else if (soitmp[i] > 8){
            redtmp[i] = 1;
        } else {
            redtmp[i] = (soitmp[i] - 3) * 0.2;  // luerd1
        }
        if (dtr[i] < 10.0){
            redrdd[i] = 1;
        } else {
            redrdd[i] = ((1 + (10 / (40 - 10) * 0.67)) - dtr[i] * (0.67 / (40 - 10)));  // luerd2
        }
        // daily photosynthetically active radiation, mj m-2 d-1
        par[i] = dtr[i] * 0.50;
        // fraction of light interception
        fint[i] = (1. - exp(-kdif * LAI[i]));
        // light use efficiency, g mj par-1
        lue1[i] = luemax * redtmp[i] * redrdd[i];
        // total intercepted photosynthetically active radiation, mj m-2 d-1
        parint[i] = fint[i] * par[i];
        // fraction of dry matter allocated to roots, kg kg-1

        // ---------------------------------------------------------------------*
        //     subroutine penman                                                *
        //     purpose: computation of the penman equation                      *
        // ---------------------------------------------------------------------*
        pevap[i] = max(0.0, np.exp(-0.5 * LAI[i]) * (penmrs[i] + penmd[i]) / lhvap);
        ptran[i] = (1. - np.exp(-0.5 * LAI[i])) * (penmrc[i] + penmd[i]) / lhvap;
        rnintc[i] = min(rain[i], 0.25 * LAI[i]);
        ptran[i] = max(0.0, ptran[i] - 0.5 * rnintc[i]);

        // ---------------------------------------------------------------------*
        //     subroutine evaptr                                                *
        //     purpose: to compute actual rates of evaporation and transpiration*
        // ---------------------------------------------------------------------*
        //wc[i] = 0.001 * wa[i] / rootd[i]
        // INSERT_RSsoilmoisture_HERE ////////////////////
        if (np.isnan(RSsms[i]){
            wc[i] = 0.001 * wa[i] / rootd[i];
        } else {
            wc[i] = RSsms[i];
        }
        ////////////////////////////////////////////////////////////////////////
        waad[i] = 1000. * wcad * rootd[i];
        wafc[i] = 1000. * wcfc * rootd[i];
        evap1[i] = pevap[i] * MAX(0., MIN(1., (wc[i] - wcad) / (wcfc - wcad)));
        if (wc[i] > wccr){
            fr[i] = MAX(0., MIN(1., (wcst - wc[i]) / (wcst - wcwet)));
        }else{
            fr[i] = MAX(0., MIN(1., (wc[i] - wcwp) / (wccr - wcwp)));
        }
        // INSERT_RStrn_HERE ////////////////////
        if (np.isnan(RStrn[i])){
            tran[i] = ptran[i] * fr[i];
        } else {
            tran[i] = RStrn[i];
        }
        ////////////////////////////////////////////////////////////////////////

        if (evap1[i] + tran[i] == 0.0){
            availf[i] = MIN(1., (wa[i] - waad[i]) / 0.01);  // notnul removed!
        }else{
            availf[i] = MIN(1., (wa[i] - waad[i]) / (evap1[i] + tran[i]));  // notnul removed!
        }
        // INSERT_RSevp_HERE ////////////////////
        if (np.isnan(RSevp[i])){
            evap[i] = evap1[i] * availf[i];
        } else {
            evap[i] = RSevp[i];
        }
        ////////////////////////////////////////////////////////////////////////

        tran[i] = tran[i] * availf[i];

        // * ---------------------------------------------------------------------*
        // *     subroutine drunir                                                *
        // *     purpose: to compute rates of drainage, runoff and irrigation     *
        // * ---------------------------------------------------------------------*
        wast[i] = 1000. * wcst * rootd[i];
        drain[i] = MAX(0, MIN(drate, (wa[i] - wafc[i] + rain[i] - rnintc[i] - evap[i] - tran[i])));
        runoff[i] = MAX(0, wa[i] - wast[i] + rain[i] - rnintc[i] - evap[i] - tran[i] - drain[i]);
        irrig[i] = irrigf * (wafc[i] - wa[i] - (rain[i] - rnintc[i] - evap[i] - tran[i] - drain[i] - runoff[i]));

        // ************************************************************************
        // ***   5. water balance and root depth growth (from lingra for thymothee)
        // ************************************************************************
        if ((rootdm - rootd[i]) <= 0 || (wc[i] - wcwp) <= 0){
            rrootd[i] = rrdmax * 0;
        } else {
            rrootd[i] = rrdmax * 1;
        }
        explor[i] = 1000. * rrootd[i] * wcfc;  // assumption that explored layers are at fc!

        if (ptran[i] <= 0){
            tranrf[i] = 1;
        } else {
            tranrf[i] = MIN(1, tran[i] / ptran[i]);  // insw removed!
        }
        rwa[i] = (rain[i] + explor[i] + irrig[i]) - (rnintc[i] + runoff[i] + tran[i] + evap[i] + drain[i]);
        frt[i] = MIN(0.263, 0.263 - tranrf[i] * (0.263 - 0.165));  // frrttb
        flv[i] = 1. - frt[i];

        // subroutine mowing
        rdaha[i] = 1;
        harv[i] = 0;
        incut[i + 1] = incut[i];

        //     mowing at observation dates
        //     reset days after harv
        if (RScut[i] == 1 && wlvg[i] > cwlvg){
            harv[i] = wlvg[i] - cwlvg;
            rdaha[i] = -daha[i];
            incut[i + 1] = incut[i] + 1.0;
        } else {
            // no mowing in current season, do not increase rate
            // of days after harv
            if (incut[i] == 0.0){
                rdaha[i] = 0;
            }
        }
        // temperature dependent leaf appearance rate, according to
        // (davies and thomas, 1983), soil temperature (soitmp)is used as
        // driving force which is estimated from a 10 day running average

        if (redtmp[i] > 0){
            leafn[i] = soitmp[i] * 0.01;
        } else {
            leafn[i] = 0;
        }
        // leaf elongation rate affected by temperature
        // cm day-1 tiller-1

        if (davtmp[i] - tmbas1 > 0){
            lera[i] = 0.83 * log(MAX(davtmp[i], 2.0)) - 0.8924;
        } else {
            lera[i] = 0;  // log10 or log?
        }
        if ((harv[i] - 0.1) < 0){
            lera2[i] = lera[i];
        } else {
            lera2[i] = -1 * length[i];
        }
        // subroutine tilsub
        dtil[i] = 0.0;
        if (daha[i] < 8){
            reftil[i] = MAX(0., 0.335 - 0.067 * LAI[i]) * redtmp[i];
        } else {
            reftil[i] = MAX(0., MIN(fsmax, 0.867 - 0.183 * LAI[i])) * redtmp[i];
        }
        // relative rate of tiller formation when defoliation less
        // than 8 days ago, tiller tiller-1 d-1
        // relative rate of tiller formation when defoliation is more
        // than 8 days ago, tiller tiller-1 d-1
        // relative death rate of tillers due to self-shading (dtild),
        // tiller tiller-1 d-1
        dtild[i] = MAX(0.01 * (1. + tsum[i] / 600.), 0.05 * (LAI[i] - LAIcr) / LAIcr);

        if (tiller[i] < 14000){
            dtil[i] = (reftil[i] - dtild[i]) * leafn[i] * tiller[i];
        } else {
            dtil[i] = -dtild[i] * leafn[i] * tiller[i];
        }
        // subroutine sosub
        lued[i] = MIN(lue1[i] * (0.336 + 0.224 * nitr) / (0.336 + 0.224 * nitmax), lue1[i] * tranrf[i]);

        // start of growing season
        if (harv[i] == 0){
            gtwso[i] = lued[i] * parint[i] * (1. + 0.8 * log(co2a / 360.0)) * 10.0;
        } else {
            gtwso[i] = 0.0;
        }
        // rate of sink limited leaf growth, unit of tiller is tillers m-2 (!),
        // 1.0e-8 is conversion from cm-2 to ha-1, ha leaf ha ground-1 d-1
        dLAIs[i] = (tiller[i] * 10000 * (lera[i] * 0.3)) * pow(1,-8);
        gtwso2[i] = gtwso[i] + wre[i];
        dre[i] = wre[i];

        // conversion to total sink limited carbon demand,
        // kg leaf ha ground-1 d-1
        if (harv[i] <= 0){
            gtwsi[i] = dLAIs[i] * (1.0 / sla) * (1.0 / flv[i]);
        } else {
            gtwsi[i] = 0;
        }
        // actual growth switches between sink- and source limitation
        // (more or less dry matter formed than can be stored)
        if ((gtwso2[i] - gtwsi[i]) > 0){
            gre[i] = gtwso2[i] - gtwsi[i];
        } else {
            gre[i] = 0;
        }
        if ((gtwso2[i] - gtwsi[i]) <= 0){
            gtw[i] = gtwso2[i];
        } else {
            gtw[i] = gtwsi[i];
        }
        // change in reserves
        rre[i] = gre[i] - dre[i];
        // relative death rate of leaves due to self-shading, d-1
        rdrsh[i] = MAX(0., MIN(0.03, 0.03 * (LAI[i] - LAIcr) / LAIcr));

        // relative death rate of leaves due to drought stress, d-1
        rdrsm[i] = MAX(0., MIN(0.05, 0.05 * (1.0 - tranrf[i])));

        // maximum of relative death rate of leaves due to
        // and drought stres, d-1
        rdrs[i] = MAX(rdrsh[i], rdrsm[i]);

        // actual relative death rate of leaves is sum of base death
        // rate plus maximum of death rates rdrsm and rdrsh, d-1
        rdr[i] = rdrd + rdrs[i];

        // actual growth rate of roots, kg ha-1 d-1
        grt[i] = gtw[i] * frt[i];

        // actual growth rate of leaf area, ha ha-1 d-1
        gLAI[i] = gtw[i] * flv[i] * sla;

        // actual death rate of leaf area, due to relative death
        // rate of leaf area or rate of change due to cutting, ha ha-1 d-1
        if (harv[i] <= 0){
            dLAI[i] = LAI[i] * (1. - exp(-rdr[i]));
        } else {
            dLAI[i] = harv[i] * sLAInt[i];
        }
        // change in LAI
        rLAI[i] = gLAI[i] - dLAI[i];
        // actual death rate of leaves, kg ha-1 d-1 incl. harvested leaves
        dlv[i] = dLAI[i] / sLAInt[i];  // notnul removed!

        // rate of change of dry weight of green leaves due to
        // growth and senescence of leaves or periodical harvest, kg ha-1 d-1
        if (harv[i] <= 0){
            glv[i] = gtw[i] * flv[i];
        } else {
            glv[i] = 0;
        }
        // change in green leaf weight
        rlv[i] = glv[i] - dlv[i];

        // INSERT_RSlai_HERE ////////////////////
        if ((i + 1) < len(RSlai)){
            if (np.isnan(RSlai[i + 1])){
                LAI[i + 1] = LAI[i] + rLAI[i];  // leaf area index, ha ha-1
            } else {
                LAI[i + 1] = RSlai[i + 1];
            }
        }
        // INSERT_RStrn_HERE ////////////////////
        if ((i + 1) < len(RStrn)){
            if (np.isnan(RStrn[i + 1])){
                tracu[i + 1] = tracu[i] + tran[i];  // cumulative transpiration, mm
            } else {
                tracu[i + 1] = tracu[i] + RStrn[i + 1];
            }
        }
        // INSERT_RSevp_HERE ////////////////////
        if ((i + 1) < len(RSevp)){
            if (np.isnan(RSevp[i + 1])){
                evacu[i + 1] = evacu[i] + evap[i];  // cumulative evaporation, mm
            } else {
                evacu[i + 1] = evacu[i] + RSevp[i];
            }
        }
        // state variables
        parcu[i + 1] = parcu[i] + parint[i];  // cumulative intercepted par, mj par intercepted m-2 ha-1
        tramcu[i + 1] = tramcu[i] + ptran[i];  // cumulative maximal transpiration, mm
        evamcu[i + 1] = evamcu[i] + pevap[i];  // cumulative maximal evaporation, mm
        irrcu[i + 1] = irrcu[i] + irrig[i];  // cumulative irrigation, mm
        tsum[i + 1] = tsum[i] + tmeff[i];  // sum of temperatures above base temperature, gr. c.d
        dvs[i] = tsum[i] / 600.0;  // hypot development stage, 600 gr. c.d taken fr subr tilsub
        daha[i + 1] = daha[i] + rdaha[i];  // days after harv, d
        tiller[i + 1] = tiller[i] + dtil[i];  // number of tillers, tillers m-2
        wlvg[i + 1] = wlvg[i] + rlv[i];  // dry weight of green leaves, kg ha-1
        wlvd[i + 1] = wlvd[i] + dlv[i];  // dry weight of dead leaves, kg ha-1 incl. harvests
        wlvd1[i] = wlvd[i] - grass[i];  // dry weight of dead leaves, kg ha-1
        hrvbl[i + 1] = wlvg[i + 1] - cwlvg;  // harvestable leaf weight
        grass[i + 1] = grass[i] + harv[i];  // dry weight of cutted green leaves, kg ha-1
        wre[i + 1] = wre[i] + rre[i];  // dry weight of storage carbohydrates, kg ha-1
        wrt[i + 1] = wrt[i] + grt[i];  // dry weight of roots, kg ha-1
        tadrw[i] = grass[i] + wlvg[i];  // total above ground dry weight including harvests,kg ha-1
        yielD[i + 1] = grass[i + 1] + MAX(0., hrvbl[i + 1]);  // harvestable part of total
        // above ground dry weight and previous harvests, kg ha-1
        length[i + 1] = length[i] + lera2[i];  // length of leaves, cm
        sLAInt[i + 1] = LAI[i + 1] / wlvg[i + 1];  // running specific leaf area in model,hakg-1 notnul rmved!
        rootd[i + 1] = rootd[i] + rrootd[i];  // rooting depth (from lingra for timothee)
        wa[i + 1] = wa[i] + rwa[i];  // soil water in rooted zone  (from lingra for timothee)
        raincu[i + 1] = raincu[i] + rain[i];  // cumulative rainfall
        dracu[i + 1] = dracu[i] + drain[i];  // cumulative drainage (avdl)
        runcu[i + 1] = runcu[i] + runoff[i];  // cumulative runoff (avdl)
        intLAI[i + 1] = intLAI[i] + rnintc[i];  // cumulative rainfall intercepted by leaves (avdl)

        // ************************************************************************
        //  ***   9. additional variables and parameters for output
        // ************************************************************************
        // newbio[i] = wlvg[i]+wrt[i]+wre[i] - (wlvgi+wrti+wrei) + grass[i]
        // if(newbio[i] == 0.):
        //    rratio[i] = 0
        // else:
        //    rratio[i] = max(0., min(1., (wrt[i]-wrti) / newbio[i]))
        // lueycu[i] = yielD[i]  / parcu[i]
        // wrtmin[i] = -wrt[i]

    }
    // print("tiller[365], yielD[365], wlvg[365], wlvd1[365], parcu[365], grass[365], tracu[365], evacu[365]")

    // data processing
    // waterinput = raincu + irrcu
    // wateroutput = tracu + evacu + intLAI + dracu + runcu

    // output in table and .csv file
    // output = [year, davtmp, dtr, dayl, LAI, tiller, yielD, wlvg, wlvd1, wre, wrt, wrtmin,
    //                rratio, sla, sLAInt, parcu, lueycu, grass, tadrw, vp, newbio, tranrf,
    //                tran, ptran, evap, pevap, wc, wa, rain, raincu, wn,tramcu, tracu, evamcu,
    //                evacu, irrcu, gtw, gtwsi, length]
    // write.csv(output, file= "m:/r/output/lingra.csv") // write output to a .csv file

    // Return last day of processing values
    results[0] = tiller[-1];
    results[1] = yielD[-1];
    results[2] = wlvg[-1];
    results[3] = wlvd1[-1];
    results[4] = wa[-1];
    results[5] = grass[-1];
    results[6] = tracu[-1];
    results[7] = evacu[-1]];
    return (EXIT_SUCCESS)
}