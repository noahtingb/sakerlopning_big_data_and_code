

import datetime
import math
import numpy as np

def _PET(ta,RH,tmrt,v,mbody,age,ht,work,icl,sex):
    """
    Args:
        ta: air temperature
        RH: relative humidity
        tmrt: Mean Radiant temperature
        v: wind at pedestrian heigh
        mbody: body masss (kg)
        age: person's age (years)
        ht: height (meters)
        work: activity level (W)
        icl: clothing amount (0-5)
        sex: 1=male 2=female
    Returns:
    """

    # humidity conversion
    vps = 6.107 * (10. ** (7.5 * ta / (238. + ta)))
    vpa = RH * vps / 100  # water vapour presure, kPa

    po = 1013.25  # Pressure
    p = 1013.25  # Pressure
    rob = 1.06
    cb = 3.64 * 1000
    food = 0
    emsk = 0.99
    emcl = 0.95
    evap = 2.42e6
    sigma = 5.67e-8
    cair = 1.01 * 1000

    eta = 0  # No idea what eta is

    c_1 = 0.
    c_2 = 0.
    c_3 = 0.
    c_4 = 0.
    c_5 = 0.
    c_6 = 0.
    c_7 = 0.
    c_8 = 0.
    c_9 = 0.
    c_10 = 0.
    c_11 = 0.

    # INBODY
    metbf = 3.19 * mbody ** (3 / 4) * (1 + 0.004 * (30 - age) + 0.018 * ((ht * 100 / (mbody ** (1 / 3))) - 42.1))
    metbm = 3.45 * mbody ** (3 / 4) * (1 + 0.004 * (30 - age) + 0.010 * ((ht * 100 / (mbody ** (1 / 3))) - 43.4))
    if sex == 1:
        met = metbm + work
    else:
        met = metbf + work

    h = met * (1 - eta)
    rtv = 1.44e-6 * met

    # sensible respiration energy
    tex = 0.47 * ta + 21.0
    eres = cair * (ta - tex) * rtv

    # latent respiration energy
    vpex = 6.11 * 10 ** (7.45 * tex / (235 + tex))
    erel = 0.623 * evap / p * (vpa - vpex) * rtv
    # sum of the results
    ere = eres + erel

    # calcul constants
    feff = 0.725
    adu = 0.203 * mbody ** 0.425 * ht ** 0.725
    facl = (-2.36 + 173.51 * icl - 100.76 * icl * icl + 19.28 * (icl ** 3)) / 100
    if facl > 1:
        facl = 1
    rcl = (icl / 6.45) / facl
    y = 1

    # should these be else if statements?
    if icl < 2:
        y = (ht-0.2) / ht
    if icl <= 0.6:
        y = 0.5
    if icl <= 0.3:
        y = 0.1

    fcl = 1 + 0.15 * icl
    r2 = adu * (fcl - 1. + facl) / (2 * 3.14 * ht * y)
    r1 = facl * adu / (2 * 3.14 * ht * y)
    di = r2 - r1
    acl = adu * facl + adu * (fcl - 1)

    tcore = [0] * 8

    wetsk = 0
    hc = 2.67 + 6.5 * v ** 0.67
    hc = hc * (p / po) ** 0.55
    c_1 = h + ere
    he = 0.633 * hc / (p * cair)
    fec = 1 / (1 + 0.92 * hc * rcl)
    htcl = 6.28 * ht * y * di / (rcl * np.log(r2 / r1) * acl)
    aeff = adu * feff
    c_2 = adu * rob * cb
    c_5 = 0.0208 * c_2
    c_6 = 0.76075 * c_2
    rdsk = 0.79 * 10 ** 7
    rdcl = 0

    count2 = 0
    j = 1

    while count2 == 0 and j < 7:
        tsk = 34
        count1 = 0
        tcl = (ta + tmrt + tsk) / 3
        count3 = 1
        enbal2 = 0

        while count1 <= 3:
            enbal = 0
            while (enbal*enbal2) >= 0 and count3 < 200:
                enbal2 = enbal
                # 20
                rclo2 = emcl * sigma * ((tcl + 273.15) ** 4 - (tmrt + 273.15) ** 4) * feff
                tsk = 1 / htcl * (hc * (tcl - ta) + rclo2) + tcl

                # radiation balance
                rbare = aeff * (1 - facl) * emsk * sigma * ((tmrt + 273.15) ** 4 - (tsk + 273.15) ** 4)
                rclo = feff * acl * emcl * sigma * ((tmrt + 273.15) ** 4 - (tcl + 273.15) ** 4)
                rsum = rbare + rclo

                # convection
                cbare = hc * (ta - tsk) * adu * (1 - facl)
                cclo = hc * (ta - tcl) * acl
                csum = cbare + cclo

                # core temperature
                c_3 = 18 - 0.5 * tsk
                c_4 = 5.28 * adu * c_3
                c_7 = c_4 - c_6 - tsk * c_5
                c_8 = -c_1 * c_3 - tsk * c_4 + tsk * c_6
                c_9 = c_7 * c_7 - 4. * c_5 * c_8
                c_10 = 5.28 * adu - c_6 - c_5 * tsk
                c_11 = c_10 * c_10 - 4 * c_5 * (c_6 * tsk - c_1 - 5.28 * adu * tsk)
                # tsk[tsk==36]=36.01
                if tsk == 36:
                    tsk = 36.01

                tcore[7] = c_1 / (5.28 * adu + c_2 * 6.3 / 3600) + tsk
                tcore[3] = c_1 / (5.28 * adu + (c_2 * 6.3 / 3600) / (1 + 0.5 * (34 - tsk))) + tsk
                if c_11 >= 0:
                    tcore[6] = (-c_10-c_11 ** 0.5) / (2 * c_5)
                if c_11 >= 0:
                    tcore[1] = (-c_10+c_11 ** 0.5) / (2 * c_5)
                if c_9 >= 0:
                    tcore[2] = (-c_7+abs(c_9) ** 0.5) / (2 * c_5)
                if c_9 >= 0:
                    tcore[5] = (-c_7-abs(c_9) ** 0.5) / (2 * c_5)
                tcore[4] = c_1 / (5.28 * adu + c_2 * 1 / 40) + tsk

                # transpiration
                tbody = 0.1 * tsk + 0.9 * tcore[j]
                sw = 304.94 * (tbody - 36.6) * adu / 3600000
                vpts = 6.11 * 10 ** (7.45 * tsk / (235. + tsk))
                if tbody <= 36.6:
                    sw = 0
                if sex == 2:
                    sw = 0.7 * sw
                eswphy = -sw * evap

                eswpot = he * (vpa - vpts) * adu * evap * fec
                wetsk = eswphy / eswpot
                if wetsk > 1:
                    wetsk = 1
                eswdif = eswphy - eswpot
                if eswdif <= 0:
                    esw = eswpot
                else:
                    esw = eswphy
                if esw > 0:
                    esw = 0

                # diffusion
                ed = evap / (rdsk + rdcl) * adu * (1 - wetsk) * (vpa - vpts)

                # MAX VB
                vb1 = 34 - tsk
                vb2 = tcore[j] - 36.6
                if vb2 < 0:
                    vb2 = 0
                if vb1 < 0:
                    vb1 = 0
                vb = (6.3 + 75 * vb2) / (1 + 0.5 * vb1)

                # energy balance
                enbal = h + ed + ere + esw + csum + rsum + food

                # clothing's temperature
                if count1 == 0:
                    xx = 1
                if count1 == 1:
                    xx = 0.1
                if count1 == 2:
                    xx = 0.01
                if count1 == 3:
                    xx = 0.001
                if enbal > 0:
                    tcl = tcl + xx
                else:
                    tcl = tcl - xx

                count3 = count3 + 1
            count1 = count1 + 1
            enbal2 = 0

        if j == 2 or j == 5:
            if c_9 >= 0:
                if tcore[j] >= 36.6 and tsk <= 34.050:
                    if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                        pass
                    else:
                        if vb > 90:
                            vb = 90
                        count2 = 1

        if j == 6 or j == 1:
            if c_11 > 0:
                if tcore[j] >= 36.6 and tsk > 33.850:
                    if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                        pass
                    else:
                        if vb > 90:
                            vb = 90
                        count2 = 1

        if j == 3:
            if tcore[j] < 36.6 and tsk <= 34.000:
                if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                    pass
                else:
                    if vb > 90:
                        vb = 90
                    count2 = 1

        if j == 7:
            if tcore[j] < 36.6 and tsk > 34.000:
                if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                    pass
                else:
                    if vb > 90:
                        vb = 90
                    count2 = 1

        if j == 4:
            if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                pass
            else:
                if vb > 90:
                    vb = 90
                count2 = 1

        j = j + 1

    # PET_cal
    tx = ta
    enbal2 = 0
    count1 = 0
    enbal = 0

    hc = 2.67 + 6.5 * 0.1 ** 0.67
    hc = hc * (p / po) ** 0.55

    while count1 <= 3:
        while (enbal * enbal2) >= 0:
            enbal2 = enbal

            # radiation balance
            rbare = aeff * (1 - facl) * emsk * sigma * ((tx + 273.15) ** 4 - (tsk + 273.15) ** 4)
            rclo = feff * acl * emcl * sigma * ((tx + 273.15) ** 4 - (tcl + 273.15) ** 4)
            rsum = rbare + rclo

            # convection
            cbare = hc * (tx - tsk) * adu * (1 - facl)
            cclo = hc * (tx - tcl) * acl
            csum = cbare + cclo

            # diffusion
            ed = evap / (rdsk + rdcl) * adu * (1 - wetsk) * (12 - vpts)

            # respiration
            tex = 0.47 * tx + 21
            eres = cair * (tx - tex) * rtv
            vpex = 6.11 * 10 ** (7.45 * tex / (235 + tex))
            erel = 0.623 * evap / p * (12 - vpex) * rtv
            ere = eres + erel

            # energy balance
            enbal = h + ed + ere + esw + csum + rsum

            # iteration concerning Tx
            if count1 == 0:
                xx = 1
            if count1 == 1:
                xx = 0.1
            if count1 == 2:
                xx = 0.01
            if count1 == 3:
                xx = 0.001
            if enbal > 0:
                tx = tx - xx
            if enbal < 0:
                tx = tx + xx
        count1 = count1 + 1
        enbal2 = 0

    return tx
def Perez_v3(zen, azimuth, radD, radI, jday):
    """
    This function calculates distribution of luminance on the skyvault based on
    Perez luminince distribution model.
    
    Created by:
    Fredrik Lindberg 20120527, fredrikl@gvc.gu.se
    Gothenburg University, Sweden
    Urban Climte Group
    
    Input parameters:
     - zen:     Zenith angle of the Sun (in degrees)
     - azimuth: Azimuth angle of the Sun (in degrees)
     - radD:    Horizontal diffuse radiation (W m-2)
     - radI:    Direct radiation perpendicular to the Sun beam (W m-2)
     - jday:    Day of year
    
    Output parameters:
     - lv:   Relative luminance map (same dimensions as theta. gamma)

    :param zen:
    :param azimuth:
    :param radD:
    :param radI:
    :param jday:
    :return:
    """
    acoeff=np.array([[1.3525,-0.2576,-0.269,-1.4366],[-1.2219,-0.773,1.4148,1.1016],[-1.1,-0.2515,0.8952,0.0156],[-0.5484,-0.6654,-0.2672,0.7117],[-0.6,-0.3566,-2.5,2.325],[-1.0156,-0.367,1.0078,1.4051],[-1.,0.0211,0.5025,-0.5119],[-1.05,0.0289,0.426,0.359]])    
    bcoeff=np.array([[-7.6700e-01,7.0000e-04,1.2734e+00,-1.2330e-01],[-2.0540e-01,3.6700e-02,-3.9128e+00,9.1560e-01],[2.7820e-01,-1.8120e-01,-4.5000e+00,1.1766e+00],[7.2340e-01,-6.2190e-01,-5.6812e+00,2.6297e+00],[2.9370e-01,4.9600e-02,-5.6812e+00,1.8415e+00],[2.8750e-01,-5.3280e-01,-3.8500e+00,3.3750e+00],[-3.0000e-01,1.9220e-01,7.0230e-01,-1.6317e+00],[-3.2500e-01,1.1560e-01,7.7810e-01,2.5000e-03]])
    ccoeff=np.array([[2.8,0.6004,1.2375,1.],[6.975,0.1774,6.4477,-0.1239],[24.7219,-13.0812,-37.7,34.8438],[33.3389,-18.3,-62.25,52.0781],[21.,-4.7656,-21.5906,7.2492],[14.,-0.9999,-7.1406,7.5469],[19.,-5.,1.2438,-1.9094],[31.0625,-14.5,-46.1148,55.375]])    
    dcoeff=np.array([[1.8734e+00,6.2970e-01,9.7380e-01,2.8090e-01],[-1.5798e+00,-5.0810e-01,-1.7812e+00,1.0800e-01],[-5.0000e+00,1.5218e+00,3.9229e+00,-2.6204e+00],[-3.5000e+00,1.6000e-03,1.1477e+00,1.0620e-01],[-3.5000e+00,-1.5540e-01,1.4062e+00,3.9880e-01],[-3.4000e+00,-1.0780e-01,-1.0750e+00,1.5702e+00],[-4.0000e+00,2.5000e-02,3.8440e-01,2.6560e-01],[-7.2312e+00,4.0500e-01,1.3350e+01,6.2340e-01]])    
    ecoeff=np.array([[0.0356,-0.1246,-0.5718,0.9938],[0.2624,0.0672,-0.219,-0.4285],[-0.0156,0.1597,0.4199,-0.5562],[0.4659,-0.3296,-0.0876,-0.0329],[0.0032,0.0766,-0.0656,-0.1294],[-0.0672,0.4016,0.3017,-0.4844],[1.0468,-0.3788,-2.4517,1.4656],[1.5,-0.6426,1.8564,0.5636]])
    
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    altitude = 90-zen
    zen = zen * deg2rad
    azimuth = azimuth * deg2rad
    altitude = altitude * deg2rad
    Idh = radD
    # Ibh = radI/sin(altitude)
    Ibn = radI

    # Skyclearness
    PerezClearness = ((Idh+Ibn)/(Idh+1.041*np.power(zen, 3)))/(1+1.041*np.power(zen, 3))
    # Extra terrestrial radiation
    day_angle = jday*2*np.pi/365
    #I0=1367*(1+0.033*np.cos((2*np.pi*jday)/365))
    I0 = 1367*(1.00011+0.034221*np.cos(day_angle) + 0.00128*np.sin(day_angle)+0.000719 *
               np.cos(2*day_angle)+0.000077*np.sin(2*day_angle))    # New from robinsson

    # Optical air mass
    # m=1/altitude; old
    if altitude >= 10*deg2rad:
        AirMass = 1/np.sin(altitude)
    elif altitude < 0:   # below equation becomes complex
        AirMass = 1/np.sin(altitude)+0.50572*np.power(180*complex(altitude)/np.pi+6.07995, -1.6364)
    else:
        AirMass = 1/np.sin(altitude)+0.50572*np.power(180*altitude/np.pi+6.07995, -1.6364)

    # Skybrightness
    # if altitude*rad2deg+6.07995>=0
    PerezBrightness = (AirMass*radD)/I0
    if Idh <= 10:
        # m_a=0;m_b=0;m_c=0;m_d=0;m_e=0;
        PerezBrightness = 0
    #if altitude < 0:
    #    print("Airmass")
    #    print(AirMass)
    #    print(PerezBrightness)
    # sky clearness bins
    if PerezClearness < 1.065:
        intClearness = 0
    elif PerezClearness < 1.230:
        intClearness = 1
    elif PerezClearness < 1.500:
        intClearness = 2
    elif PerezClearness < 1.950:
        intClearness = 3
    elif PerezClearness < 2.800:
        intClearness = 4
    elif PerezClearness < 4.500:
        intClearness = 5
    elif PerezClearness < 6.200:
        intClearness = 6
    elif PerezClearness > 6.200:
        intClearness = 7
    else:
        raise ValueError('No valid PerezClearness, are inputs NaN?')

    m_a = acoeff[intClearness,  0] + acoeff[intClearness,  1] * zen + PerezBrightness * (acoeff[intClearness,  2] + acoeff[intClearness,  3] * zen)
    m_b = bcoeff[intClearness,  0] + bcoeff[intClearness,  1] * zen + PerezBrightness * (bcoeff[intClearness,  2] + bcoeff[intClearness,  3] * zen)
    m_e = ecoeff[intClearness,  0] + ecoeff[intClearness,  1] * zen + PerezBrightness * (ecoeff[intClearness,  2] + ecoeff[intClearness,  3] * zen)

    if intClearness > 0:
        m_c = ccoeff[intClearness, 0] + ccoeff[intClearness, 1] * zen + PerezBrightness * (ccoeff[intClearness, 2] + ccoeff[intClearness, 3] * zen)
        m_d = dcoeff[intClearness, 0] + dcoeff[intClearness, 1] * zen + PerezBrightness * (dcoeff[intClearness, 2] + dcoeff[intClearness, 3] * zen)
    else:
        # different equations for c & d in clearness bin no. 1,  from Robinsson
        m_c = np.exp(np.power(PerezBrightness * (ccoeff[intClearness, 0] + ccoeff[intClearness, 1] * zen), ccoeff[intClearness, 2]))-1
        m_d = -np.exp(PerezBrightness * (dcoeff[intClearness, 0] + dcoeff[intClearness, 1] * zen)) + dcoeff[intClearness, 2] + \
            PerezBrightness * dcoeff[intClearness, 3] * PerezBrightness

    skyvaultalt = np.atleast_2d([])
    skyvaultazi = np.atleast_2d([])

    if True:
        # Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
        skyvaultaltint = [6, 18, 30, 42, 54, 66, 78]
        skyvaultaziint = [12, 12, 15, 15, 20, 30, 60]
        for j in range(7):
            for k in range(1, int(360/skyvaultaziint[j]) + 1):
                skyvaultalt = np.append(skyvaultalt, skyvaultaltint[j])
                skyvaultazi = np.append(skyvaultazi, k*skyvaultaziint[j])

        skyvaultalt = np.append(skyvaultalt, 90)
        skyvaultazi = np.append(skyvaultazi, 360)

    skyvaultzen = (90 - skyvaultalt) * deg2rad
    skyvaultalt = skyvaultalt * deg2rad
    skyvaultazi = skyvaultazi * deg2rad

    # Angular distance from the sun from Robinsson
    cosSkySunAngle = np.sin(skyvaultalt) * np.sin(altitude) + \
                     np.cos(altitude) * np.cos(skyvaultalt) * np.cos(np.abs(skyvaultazi-azimuth))

    # Main equation
    lv = (1 + m_a * np.exp(m_b / np.cos(skyvaultzen))) * ((1 + m_c * np.exp(m_d * np.arccos(cosSkySunAngle)) +
                                                           m_e * cosSkySunAngle * cosSkySunAngle))

    # Normalisation
    lv = lv / np.sum(lv)

    #x = np.atleast_2d([])
    #lv = np.transpose(np.append(np.append(np.append(x, skyvaultalt*rad2deg), skyvaultazi*rad2deg), lv))
    x = np.transpose(np.atleast_2d(skyvaultalt*rad2deg))
    y = np.transpose(np.atleast_2d(skyvaultazi*rad2deg))
    z = np.transpose(np.atleast_2d(lv))
    lv = np.append(np.append(x, y, axis=1), z, axis=1)
    return lv, PerezClearness, PerezBrightness

def Lside_veg_v2020a(azimuth,altitude,Ta,Tw,SBC,ewall,Ldown,esky,t,F_sh,CI,Lup):
    
    # This m-file is the current one that estimates L from the four cardinal points 20100414
    #Building height angle from svf
    svfalfa=np.arcsin(np.exp((np.log(1-0.6))/2))
    
    aziW=azimuth+t
    aziN=azimuth-90+t
    aziE=azimuth-180+t
    aziS=azimuth-270+t
    
    F_sh = 2*F_sh-1  #(cylindric_wedge scaled 0-1)
    
    c=1-CI
    Lsky_allsky = esky*SBC*((Ta+273.15)**4)*(1-c)+c*SBC*((Ta+273.15)**4)
    
    viktveg, viktwall, viktsky, viktrefl = -1.1102230246251565e-16, 0.719056214891864, 0.28094378510813617, 0.7190562148918639
    alfaB=np.arctan(svfalfa)
    betaB=np.arctan(np.tan((svfalfa)*F_sh))
    betasun=((alfaB-betaB)/2)+betaB
    Lsky = (0.6 * Lsky_allsky) * viktsky * 0.5
    Lrefl = (Ldown + Lup) * (viktrefl) * (1 - ewall) * 0.5
    Lground = Lup * 0.5
    Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5

    if altitude > 0:  # daytime
        # betasun = np.arctan(0.5*np.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
        if (azimuth > (180-t))  and  (azimuth <= (360-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziE*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    
    Lsum = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    if altitude>0: # daytime
        if (azimuth <= (90-t))  or  (azimuth > (270-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziS*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum+= Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        if (azimuth > (360-t))  or  (azimuth <= (180-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziW*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum += Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        if (azimuth > (90-t))  and  (azimuth <= (270-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziN*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum += Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        return Lsum
    return Lsum*4

def Kside_veg_v2019a(radI, radD, radG,altitude, albedo, F_sh, Kup, lv, diffsh):
    # New reflection equation 2012-05-25
    
    deg2rad = np.pi / 180
    KsideD = 0

    # Direct radiation ###
    # Kside with cylinder ###
    KsideI = radI * np.cos(altitude * deg2rad)

    ### Diffuse and reflected radiation ###
#    viktveg, viktwall = 4.440892098500626e-16, 0.7190562148918634
    svfviktbuveg = 0.7190562148918634# + (4.440892098500626e-16) * (1 - 1)

    ### Anisotropic Diffuse Radiation after Perez et al. 1993 ###

    aniAlt = lv[0][:, 0]
    aniLum = lv[0][:, 2]

    phiVar = np.zeros((145, 1))
    radTot = np.zeros(1)

    for ix in range(0, 145):  # Azimuth delta
        if ix < 60:
            aziDel = 12
        elif ix >= 60 and ix < 108:
            aziDel = 15
        elif ix >= 108 and ix < 126:
            aziDel = 20
        elif ix >= 126 and ix < 138:
            aziDel = 30
        elif ix >= 138 and ix < 144:
            aziDel = 60
        elif ix == 144:
            aziDel = 360

        phiVar[ix] = (aziDel * deg2rad) * (np.sin((aniAlt[ix] + 6) * deg2rad) - np.sin((aniAlt[ix] - 6) * deg2rad))  # Solid angle / Steradian
        radTot = radTot + (aniLum[ix] * phiVar[ix] * np.sin(aniAlt[ix] * deg2rad))  # Radiance fraction normalization
    lumChi = (aniLum * radD) / radTot  # Radiance fraction normalization
    for idx in range(0, 145):
            anglIncC = np.cos(aniAlt[idx] * deg2rad)  # Angle of incidence, np.cos(0) because cylinder - always perpendicular
            KsideD += diffsh[idx] * lumChi[idx] * anglIncC * phiVar[idx]  # Diffuse vertical radiation
    Keast = (albedo * (svfviktbuveg * (radG * (1 - F_sh) + radD * F_sh)) + Kup) * 0.5
    return Keast, KsideI, KsideD

def cylindric_wedge(zen, svfalfa):

    # Fraction of sunlit walls based on sun altitude and svf wieghted building angles
    # input: 
    # sun zenith angle "beta"
    # svf related angle "alfa"

    beta=zen
    alfa= svfalfa
    
    # measure the size of the image
    # sizex=size(svfalfa,2)
    # sizey=size(svfalfa,1)
    
    xa=1-2./(np.tan(alfa)*np.tan(beta))
    ha=2./(np.tan(alfa)*np.tan(beta))
    ba=(1./np.tan(alfa))
    hkil=2.*ba*ha
    if xa<0:
        qa=np.tan(beta)/2
        Za=((ba**2)-((qa**2)/4))**0.5
        phi=np.arctan(Za/qa)
        A=(np.sin(phi)-phi*np.cos(phi))/(1-np.cos(phi))
        ukil=2*ba*xa*A
    else:
        qa,Za,phi,A,ukil=0,0,0,0,0    
    
    Ssurf=hkil+ukil
    
    F_sh=(2*np.pi*ba-Ssurf)/(2*np.pi*ba)#Xa
    
    return F_sh
def diffusefraction(radG,altitude,Kt,Ta,RH):
    """
    This function estimates diffuse and directbeam radiation according to
    Reindl et al (1990), Solar Energy 45:1

    :param radG:
    :param altitude:
    :param Kt:
    :param Ta:
    :param RH:
    :return:
    """

    alfa = altitude*(np.pi/180)

    if Ta <= -999.00 or RH <= -999.00 or np.isnan(Ta) or np.isnan(RH):
        if Kt <= 0.3:
            radD = radG*(1.020-0.248*Kt)
        elif Kt > 0.3 and Kt < 0.78:
            radD = radG*(1.45-1.67*Kt)
        else:
            radD = radG*0.147
    else:
        RH = RH/100
        if Kt <= 0.3:
            radD = radG*(1 - 0.232 * Kt + 0.0239 * np.sin(alfa) - 0.000682 * Ta + 0.0195 * RH)
        elif Kt > 0.3 and Kt < 0.78:
            radD = radG*(1.329- 1.716 * Kt + 0.267 * np.sin(alfa) - 0.00357 * Ta + 0.106 * RH)
        else:
            radD = radG*(0.426 * Kt - 0.256 * np.sin(alfa) + 0.00349 * Ta + 0.0734 * RH)

    radI = (radG - radD)/(np.sin(alfa))

    # Corrections for low sun altitudes (20130307)
    if radI < 0:
        radI = 0
    if altitude < 1 and radI > radG:
        radI=radG
    if radD > radG:
        radD = radG
    return radI, radD

def daylen(doy, xlat,xlang=0):
    # Calculation of declination of sun (Eqn. 16). Amplitude= +/-23.45
    # deg. Minimum = DOY 355 (DEC 21), maximum = DOY 172.5 (JUN 21/22).
    # Sun angles.  SOC limited for latitudes above polar circles.
    # Calculate daylength, sunrise and sunset (Eqn. 17)
    #doy is day of year, xlat is the latitude in degrees
    rad=np.pi/180.0
    DEC = -23.45 * np.cos(2.0*np.pi*(doy+10.0)/365.0)
    SOC = np.tan(rad*DEC) * np.tan(rad*xlat)
    SOC = min(max(SOC,-1.0),1.0) #correction for polar circles
    # SOC=alt
    day_len = 12 + 24.0*np.arcsin(SOC)/np.pi
    sun_upp = 12.0-0/15 - day_len/2.0 
    sun_down = 12.0-0/15 + day_len/2.0
    return day_len, DEC, sun_down, sun_upp


def Solweig1D_2020a_calc(absK, absL, Fside, Fup, Fcyl,location, Ta, RH, radG, radD, radI, P,radD1,radI1,altitude, azimuth, zen, jday, dectime, altmax):
    # This is the core function of the SOLWEIG1D model, 2019-Jun-21
    # Fredrik Lindberg, fredrikl@gvc.gu.se, Goteborg Urban Climate Group, Gothenburg University, Sweden
    # Instrument offset in degrees
    
    # If metfile starts at night
    skyvaultalt = np.atleast_2d([])
    skyvaultaltint = [6, 18, 30, 42, 54, 66, 78]
    skyvaultaziint = [12, 12, 15, 15, 20, 30, 60]
    for j in range(7):
        for k in range(1, int(360/skyvaultaziint[j]) + 1):
            skyvaultalt = np.append(skyvaultalt, skyvaultaltint[j])

    skyvaultalt = np.append(skyvaultalt, 90)

    diffsh = np.zeros((145))
    svfalfadeg = 0.6847192 / (np.pi / 180.)
    for k in range(0, 145):
        if skyvaultalt[k] > svfalfadeg:
            diffsh[k] = 1

    t = 0.
    # Stefan Bolzmans Constant
    SBC = 5.67051e-8
    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(jday, location['latitude'],location["longitude"])
    # Vapor pressure
    ea = 6.107 * 10 ** ((7.5 * Ta) / (237.3 + Ta)) * (RH / 100.)
    # Determination of clear - sky emissivity from Prata (1996)
    msteg = 46.5 * (ea / (Ta + 273.15))
    esky = (1 - (1 + msteg) * np.exp(-((1.2 + 3.0 * msteg) ** 0.5))) + 0  # -0.04 old error from Jonsson et al.2006

    if altitude > 0: # # # # # # DAYTIME # # # # # #
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        #  factor for low sun elevations after Lindberg et al.(2008)

        # Estimation of radD and radI if not measured after Reindl et al.(1990)        
        I0, CI, Kt = clearnessindex_2013b(zen, jday, Ta, RH / 100., radG, location)
        if (CI > 1) or (CI == np.inf):
            CI = 1

        radI, radD = diffusefraction(radG, altitude, Kt, Ta, RH)
        #print("he",radI,radD)
        #print("he1",radI1,radD1)

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        zenDeg = zen*(180/np.pi)
        lv = Perez_v3(zenDeg, azimuth, radD, radI, jday)   # Relative luminance

        aniLum = 0.
        for idx in range(0, 145):
            aniLum = aniLum + diffsh[idx] * lv[0][idx][2]     # Total relative luminance from sky into each cell

        dRad = aniLum * radD   # Total diffuse radiation from sky into each cell
        #print("he1",dRad)


        # # # Surface temperature parameterisation during daytime # # # #
        # new using max sun alt.instead of  dfm

        const=np.sin(((dectime/1 - SNUP / 24) / (15 / 24 - (SNUP+0/15) / 24)) * np.pi / 2) -3.41
        Tg = max(0.37 * altmax  *const -3.41,0) # 2015 a, based on max sun altitude # temporary for removing low Tg during morning 20130205
        Tgwall = 0.58 * altmax * const -3.41 # 2015a, based on max sun altitude
        #print("TG",Tg,Tgwall)
        # New estimation of Tg reduction for non - clear situation based on Reindl et al.1990
        radI0, _ = diffusefraction(I0, altitude, 1., Ta, RH)

        corr = 0.1473 * np.log(90 - (zen / np.pi * 180)) + 0.3454  # 20070329 correction of lat, Lindberg et al. 2008
        CI_Tg = (radI / radI0) + (1 - corr)

        if (CI_Tg > 1) or (CI_Tg == np.inf):
            CI_Tg = 1
        Tg = Tg * CI_Tg  # new estimation
        Tgwall = Tgwall * CI_Tg

        if Tg < 0.:
            Tg = 0.
        #Tg[Tg < 0] = 0  # temporary for removing low Tg during morning 20130205

        gvf = 1

        Lup = SBC * 0.95 * ((gvf + Ta + Tg + 273.15) ** 4)

        # Building height angle from svfs
        F_sh = float(cylindric_wedge(zen, 0.6847192))  # Fraction shadow on building walls based on sun alt and svf
        #F_sh[np.isnan(F_sh)] = 0.5

        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        albedo_b=0.2
        Kdown = radI * np.sin(altitude * (np.pi / 180)) + dRad  # *sin(altitude(i) * (pi / 180))
        
        Kup =  0.15 * (radI * np.sin(altitude * (np.pi / 180.))) + dRad

        Kesnw, KsideI, KsideD = Kside_veg_v2019a(radI, radD, radG, altitude, albedo_b, F_sh, Kup, lv, diffsh)
            # # # # Ldown # # # # 
        ewall=0.9   
        Ldown = (0.6) * esky * SBC * ((Ta + 273.15) ** 4) + (0.4) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) + \
                (0.4) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)

        if CI < 0.95:  # non - clear conditions
            c = 1 - CI
            Ldown = Ldown * (1 - c) + c * (0.6 * SBC * ((Ta + 273.15) ** 4) + (0.4) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) +
                        (0.4) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4))

        # # # # Lside # # # # FIXAD
        Lsum = Lside_veg_v2020a(azimuth, altitude, Ta, Tgwall, SBC, ewall, Ldown,esky, t, F_sh, CI, Lup)

        # # # # Calculation of radiant flux density and Tmrt # # # #
        # Human body considered as a cylinder with Perez et al. (1993)
        Sstr = absK * ((KsideI + KsideD) * Fcyl + (Kdown + Kup) * Fup + (4*Kesnw) * Fside) + absL * ((Ldown + Lup) * Fup + (Lsum) * Fside) #4 directions

        Tmrt = float(np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.15)
        return Tmrt
    
    else:  # # # # # # # NIGHTTIME # # # # # # # #
        # # # # Lup # # # #
        Lup = SBC * 0.95 * ((0 + Ta + 273.15) ** 4)
        # # # # Ldown # # # # 
        ewall=0.9   
        Ldown = (0.6) * esky * SBC * ((Ta + 273.15) ** 4) + (0.4) * ewall * SBC * ((Ta + 273.15) ** 4) + \
                (0.4) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)
        if 0 < 0.95:  # non - clear conditions
            c = 1
            Ldown =  1 * (0.6 * SBC * ((Ta + 273.15) ** 4) + (0.4) * ewall * SBC * ((Ta + 273.15) ** 4) +
                        (0.4) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4))

        # # # # Lside # # # # FIXAD
        Lsum = Lside_veg_v2020a(azimuth, altitude, Ta, 0, SBC, ewall, Ldown,esky, 0, 0, 0, Lup)

        # # # # Calculation of radiant flux density and Tmrt # # # #
        # Human body considered as a cylinder with Perez et al. (1993)
        Sstr = absL * ((Ldown + Lup) * Fup + (Lsum) * Fside) #4 directions

        Tmrt = float(np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.15)
        return Tmrt
def sun_distance(jday):
    """
    #% Calculatesrelative earth sun distance
    #% with day of year as input.
    #% Partridge and Platt, 1975
    """
    b = 2.*np.pi*jday/365.
    return np.sqrt((1.00011+np.dot(0.034221, np.cos(b))+np.dot(0.001280, np.sin(b))+np.dot(0.000719,
                                        np.cos((2.*b)))+np.dot(0.000077, np.sin((2.*b)))))
def clearnessindex_2013b(zen, jday, Ta, RH, radG, location):

    """ Clearness Index at the Earth's surface calculated from Crawford and Duchon 1999

    :param zen: zenith angle in radians
    :param jday: day of year
    :param Ta: air temperature
    :param RH: relative humidity
    :param radG: global shortwave radiation
    :param location: distionary including lat, lon and alt
    :param P: pressure
    :return:
    """
    #P=1013

    Itoa = 1370.0  # Effective solar constant
    D = sun_distance(jday)  # irradiance differences due to Sun-Earth distances
    m = 35. * np.cos(zen) * ((1224. * (np.cos(zen)**2) + 1) ** (-1/2.))     # optical air mass at p=1013
    Trpg = 1.021-0.084*(m*(0.000949*1013+0.051))**0.5  # Transmission coefficient for Rayliegh scattering and permanent gases

    # empirical constant depending on latitude
    #print(location['latitude']//10)
    G=[[3.37, 2.85, 2.80, 2.64],[2.99, 3.02, 2.70, 2.93],[3.60, 3.00, 2.98, 2.93], [3.04, 3.11, 2.92, 2.94], [2.70, 2.95, 2.77, 2.71], [2.52, 3.07, 2.67, 2.93],[1.76, 2.69, 2.61, 2.61],[1.60, 1.67, 2.24, 2.63],[1.11, 1.44, 1.94, 2.02]][max(0,int(location['latitude']//10))]
        
    if jday > 335 or jday <= 60:
        G = G[0]
    elif jday > 60 and jday <= 152:
        G = G[1]
    elif jday > 152 and jday <= 244:
        G = G[2]
    elif jday > 244 and jday <= 335:
        G = G[3]

    # dewpoint calculation
    a2 = 17.27
    b2 = 237.7
    Td = (b2*(((a2*Ta)/(b2+Ta))+np.log(RH)))/(a2-(((a2*Ta)/(b2+Ta))+np.log(RH)))
    Td = (Td*1.8)+32  # Dewpoint (F)
    u = np.exp(0.1133-np.log(G+1)+0.0393*Td)  # Precipitable water
    Tw = 1-0.077*((u*m)**0.3)  # Transmission coefficient for water vapor
    Tar = 0.935**m  # Transmission coefficient for aerosols

    I0=Itoa*np.cos(zen)*Trpg*Tw*D*Tar
    if abs(zen)>np.pi/2:
        I0 = 0
    if not(np.isreal(I0)):
        I0 = 0

    corr=0.1473*np.log(90-(zen/np.pi*180))+0.3454  # 20070329

    CIuncorr = radG / I0
    CI = CIuncorr + (1-corr)
    I0et = Itoa*np.cos(zen)*D  # extra terrestial solar radiation 
    Kt = radG / I0et
    if math.isnan(CI):
        CI = float('Inf')
    return I0, CI, Kt

def Solweig_2015a_metdata_noload(year,doy,hour,minu, location, UTC=0):
    """
    This function is used to process the input meteorological file.
    It also calculates Sun position based on the time specified in the met-file

    :param inputdata:
    :param location:
    :param UTC:
    :return:
    """
    dectime = doy+hour / 24 + minu / (60*24.)
    halftimestepdec = 0
    sunmaximum = 0.
    
    sunmax = dict()

    YMD = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(doy) - 1)
    # Finding maximum altitude in 15 min intervals (20141027)
    fifteen = 0.
    sunmaximum = -90.
    sunmax={'zenith': 90.}
    while sunmaximum <= 90. - sunmax['zenith']:
                sunmaximum = 90. - sunmax['zenith']
                fifteen = fifteen + 15. / 1440.
                YMDHM = YMD+datetime.timedelta(days=(60*10)/1440.0 + fifteen)
                time = {"year":YMDHM.year,'month': YMDHM.month,'day':YMDHM.day,'hour': YMDHM.hour,'min':YMDHM.minute,'sec':0,"UTC":UTC}
                sunmax = sun_position(time,location)
                #print("sunmax",sunmax,time)
    altmax = sunmaximum 
    #YMDHM = YMD +datetime.timedelta(hours=12)
    #print(YMDHM,location["longitude"],location["latitude"],time)
    #print("sunmax",sunmax,sun_position({"year":YMDHM.year,'month': YMDHM.month,'day':YMDHM.day,'hour': YMDHM.hour,'min':YMDHM.minute,'sec':0,"UTC":UTC},location))

    time = {"year":YMDHM.year,'month': YMDHM.month,'day':YMDHM.day,'hour': YMDHM.hour,'min':YMDHM.minute,'sec':0,"UTC":UTC}
    sun = sun_position(time, location)
    altitude = 90. - sun['zenith']
    #print("altmax",altmax,altitude)
    azimuth = sun['azimuth']
    zen = sun['zenith'] * (3.141592653589793238/180.)
        # day of year and check for leap year
    if (time['year'] % 4) == 0 and ( ((time['year'] % 100)==0 and (time['year'] % 400) == 0) or ((time['year'] % 100)!=0)):
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    jday = sum(dayspermonth[0:time['month'] - 1]) + time['day']

    return None, altitude, azimuth, zen, jday, None, dectime, altmax

def indexflask(form):
        month = int(form["month"])
        day = int(form["day"])
        hour = int(form["hour"])
        year = int(form["year"])
        minu = 30
        Ta = float(form["Ta"])        
        RH = float(form["RH"])
        Ws = float(form["Ws"])
        radD1,radI1=float(form["radD"]),float(form["radI"])
        location = form["loc"]

        if month > 12 or month < 0:
            print("petresult.html","Incorrect month filled in")
        if day > 31 or day < 0:
            print("petresult.html","Incorrect day filled in")
        if hour > 23 or hour < 0:
            print("petresult.html","Incorrect hour filled in")
        if Ta > 60 or Ta < -75:
            print("petresult.html", "Unreasonable air temperature filled in",Ta)
        if RH > 100 or RH < 0:
            print("petresult.html", "Unreasonable relative humidity filled in")
        if Ws > 100 or Ws < 0:
            print("petresult.html", "Unreasonable Wind speed filled in")
        #day of year
        if (year % 4) == 0 and ( ((year % 100)==0 and (year % 400) == 0) or ((year % 100)!=0)):
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        doy= sum(dayspermonth[0:month - 1]) + day
        # Currently looked to England
        # Radiation
        _, altitude, azimuth, zen, jday, _, dectime, altmax = Solweig_2015a_metdata_noload(year,doy,hour,minu, location, 0)
        #print(altitude,zen,jday)
        
        
        if altitude > 0.:
            I0, _, _= clearnessindex_2013b(zen, jday, Ta, RH / 100., 40, location)
            radG = I0[0]*(1)#clear
        else:
            radG = 0.
        #print(radG)
        #print(radD1,radI1)

        # Calculate the global iradiation
        # Main calculation
        if Ta is not None and RH is not None and Ws is not None and radI1+radD1 is not None:
            Tmrt, resultPET, _ = petcalc(Ta, RH, Ws, radG,location,radD1,radI1,altitude, azimuth, zen, jday, dectime, altmax)
        return resultPET,Tmrt

def petcalc(Ta, RH, Ws, radG,location,radD1,radI1,altitude, azimuth, zen, jday, dectime, altmax):
#    sh = 1.  # 0 if shadowed by building
#    vegsh = 1.  # 0 if shadowed by tree

    # Location and time settings. Should be moved out later on

    # Human parameter data. Should maybe be move out later on
    pos=0

    if pos == 0:
        Fside = 0.22
        Fup = 0.06
        Fcyl = 0.28
    else:
        Fside = 0.166666
        Fup = 0.166666
        Fcyl = 0.2

    radD = -999.
    radI = -999.
    P   = -999.

    # main loop
    # Nocturnal cloudfraction from Offerle et al. 2003
    #cant have this funktion with this data

    Tmrt = Solweig1D_2020a_calc(0.70, 0.95,Fside, Fup, Fcyl,location,Ta, RH, radG, radD, radI, P,radD1,radI1,altitude, azimuth, zen, jday, dectime, altmax)

    # Recalculating wind speed based on pwerlaw
    WsPET = (1.1 / 10) ** 0.2 * Ws

    mbody,ht,clo,age,activity,sex=75.,1.8,0.9,35,80,1
    resultPET = _PET(Ta, RH, Tmrt, WsPET, mbody, age, ht, activity, clo, sex)

    return Tmrt, resultPET, None

def sun_position(time, location):
    """
    % sun = sun_position(time, location)
    %
    % This function compute the sun position (zenith and azimuth angle at the observer
    % location) as a function of the observer local time and position.
    %
    % It is an implementation of the algorithm presented by Reda et Andreas in:
    %   Reda, I., Andreas, A. (2003) Solar position algorithm for solar
    %   radiation application. National Renewable Energy Laboratory (NREL)
    %   Technical report NREL/TP-560-34302.
    % This document is avalaible at www.osti.gov/bridge
    %
    % This algorithm is based on numerical approximation of the exact equations.
    % The authors of the original paper state that this algorithm should be
    % precise at +/- 0.0003 degrees. I have compared it to NOAA solar table
    % (http://www.srrb.noaa.gov/highlights/sunrise/azel.html) and to USNO solar
    % table (http://aa.usno.navy.mil/data/docs/AltAz.html) and found very good
    % correspondance (up to the precision of those tables), except for large
    % zenith angle, where the refraction by the atmosphere is significant
    % (difference of about 1 degree). Note that in this code the correction
    % for refraction in the atmosphere as been implemented for a temperature
    % of 10C (283 kelvins) and a pressure of 1010 mbar. See the subfunction
    % �sun_topocentric_zenith_angle_calculation� for a possible modification
    % to explicitely model the effect of temperature and pressure as describe
    % in Reda & Andreas (2003).
    %
    % Input parameters:
    %   time: a structure that specify the time when the sun position is
    %   calculated.
    %       time.year: year. Valid for [-2000, 6000]
    %       time.month: month [1-12]
    %       time.day: calendar day [1-31]
    %       time.hour: local hour [0-23]
    %       time.min: minute [0-59]
    %       time.sec: second [0-59]
    %       time.UTC: offset hour from UTC. Local time = Greenwich time + time.UTC
    %   This input can also be passed using the Matlab time format ('dd-mmm-yyyy HH:MM:SS').
    %   In that case, the time has to be specified as UTC time (time.UTC = 0)
    %
    %   location: a structure that specify the location of the observer
    %       location.latitude: latitude (in degrees, north of equator is
    %       positive)
    %       location.longitude: longitude (in degrees, positive for east of
    %       Greenwich)
    %       location.altitude: altitude above mean sea level (in meters)
    %
    % Output parameters
    %   sun: a structure with the calculated sun position
    %       sun.zenith = zenith angle in degrees (angle from the vertical)
    %       sun.azimuth = azimuth angle in degrees, eastward from the north.
    % Only the sun zenith and azimuth angles are returned as output, but a lot
    % of other parameters are calculated that could also extracted as output of
    % this function.
    %
    % Exemple of use
    %
    % location.longitude = -105.1786;
    % location.latitude = 39.742476;
    % location.altitude = 1830.14;
    % time.year = 2005;
    % time.month = 10;
    % time.day = 17;
    % time.hour = 6;
    % time.min = 30;
    % time.sec = 30;
    % time.UTC = -7;
    % %
    % location.longitude = 11.94;
    % location.latitude = 57.70;
    % location.altitude = 3.0;
    % time.UTC = 1;
    % sun = sun_position(time, location);
    %
    % sun =
    %
    %      zenith: 50.1080438859849
    %      azimuth: 194.341174010338
    %
    % History
    %   09/03/2004  Original creation by Vincent Roy (vincent.roy@drdc-rddc.gc.ca)
    %   10/03/2004  Fixed a bug in julian_calculation subfunction (was
    %               incorrect for year 1582 only), Vincent Roy
    %   18/03/2004  Correction to the header (help display) only. No changes to
    %               the code. (changed the �elevation� field in �location� structure
    %               information to �altitude�), Vincent Roy
    %   13/04/2004  Following a suggestion from Jody Klymak (jklymak@ucsd.edu),
    %               allowed the 'time' input to be passed as a Matlab time string.
    %   22/08/2005  Following a bug report from Bruce Bowler
    %               (bbowler@bigelow.org), modified the julian_calculation function. Bug
    %               was 'MATLAB has allowed structure assignment  to a non-empty non-structure
    %               to overwrite the previous value.  This behavior will continue in this release,
    %               but will be an error in a  future version of MATLAB.  For advice on how to
    %               write code that  will both avoid this warning and work in future versions of
    %               MATLAB,  see R14SP2 Release Notes'. Script should now be
    %               compliant with futher release of Matlab...
    """

    # 1. Calculate the Julian Day, and Century. Julian Ephemeris day, century
    # and millenium are calculated using a mean delta_t of 33.184 seconds.
    julian = julian_calculation(time)
    #print(julian)

    # 2. Calculate the Earth heliocentric longitude, latitude, and radius
    # vector (L, B, and R)
    earth_heliocentric_position = earth_heliocentric_position_calculation(julian)

    # 3. Calculate the geocentric longitude and latitude
    sun_geocentric_position = sun_geocentric_position_calculation(earth_heliocentric_position)

    # 4. Calculate the nutation in longitude and obliquity (in degrees).
    nutation = nutation_calculation(julian)

    # 5. Calculate the true obliquity of the ecliptic (in degrees).
    true_obliquity = true_obliquity_calculation(julian, nutation)

    # 6. Calculate the aberration correction (in degrees)
    aberration_correction = abberation_correction_calculation(earth_heliocentric_position)

    # 7. Calculate the apparent sun longitude in degrees)
    apparent_sun_longitude = apparent_sun_longitude_calculation(sun_geocentric_position, nutation, aberration_correction)

    # 8. Calculate the apparent sideral time at Greenwich (in degrees)
    apparent_stime_at_greenwich = apparent_stime_at_greenwich_calculation(julian, nutation, true_obliquity)

    # 9. Calculate the sun rigth ascension (in degrees)
    sun_rigth_ascension = sun_rigth_ascension_calculation(apparent_sun_longitude, true_obliquity, sun_geocentric_position)

    # 10. Calculate the geocentric sun declination (in degrees). Positive or
    # negative if the sun is north or south of the celestial equator.
    sun_geocentric_declination = sun_geocentric_declination_calculation(apparent_sun_longitude, true_obliquity,
                                                                        sun_geocentric_position)

    # 11. Calculate the observer local hour angle (in degrees, westward from south).
    observer_local_hour = observer_local_hour_calculation(apparent_stime_at_greenwich, location, sun_rigth_ascension)

    # 12. Calculate the topocentric sun position (rigth ascension, declination and
    # rigth ascension parallax in degrees)
    topocentric_sun_position = topocentric_sun_position_calculate(earth_heliocentric_position, location,
                                                                  observer_local_hour, sun_rigth_ascension,
                                                                  sun_geocentric_declination)

    # 13. Calculate the topocentric local hour angle (in degrees)
    topocentric_local_hour = topocentric_local_hour_calculate(observer_local_hour, topocentric_sun_position)

    # 14. Calculate the topocentric zenith and azimuth angle (in degrees)
    sun = sun_topocentric_zenith_angle_calculate(location, topocentric_sun_position, topocentric_local_hour)

    return sun


def julian_calculation(t_input):
    """
    % This function compute the julian day and julian century from the local
    % time and timezone information. Ephemeris are calculated with a delta_t=0
    % seconds.

    % If time input is a Matlab time string, extract the information from
    % this string and create the structure as defined in the main header of
    % this script.
    """
    if not isinstance(t_input, dict):
        # tt = datetime.datetime.strptime(t_input, "%Y-%m-%d %H:%M:%S.%f")    # if t_input is a string of this format
        # t_input should be a datetime object
        time = dict()
        time['UTC'] = 0
        time['year'] = t_input.year
        time['month'] = t_input.month
        time['day'] = t_input.day
        time['hour'] = t_input.hour
        time['min'] = t_input.minute
        time['sec'] = t_input.second
    else:
        time = t_input

    if time['month'] == 1 or time['month'] == 2:
        Y = time['year'] - 1
        M = time['month'] + 12
    else:
        Y = time['year']
        M = time['month']

    ut_time = ((time['hour'] - time['UTC'])/24) + (time['min']/(60*24)) + (time['sec']/(60*60*24))   # time of day in UT time.
    D = time['day'] + ut_time   # Day of month in decimal time, ex. 2sd day of month at 12:30:30UT, D=2.521180556

    # In 1582, the gregorian calendar was adopted
    if time['year'] == 1582:
        if time['month'] == 10:
            if time['day'] <= 4:   # The Julian calendar ended on October 4, 1582
                B = (0)
            elif time['day'] >= 15:   # The Gregorian calendar started on October 15, 1582
                A = np.floor(Y/100)
                B = 2 - A + np.floor(A/4)
            else:
                print('This date never existed!. Date automatically set to October 4, 1582')
                time['month'] = 10
                time['day'] = 4
                B = 0
        elif time['month'] < 10:   # Julian calendar
            B = 0
        else: # Gregorian calendar
            A = np.floor(Y/100)
            B = 2 - A + np.floor(A/4)
    elif time['year'] < 1582:   # Julian calendar
        B = 0
    else:
        A = np.floor(Y/100)    # Gregorian calendar
        B = 2 - A + np.floor(A/4)

    julian = dict()
    julian['day'] = D + B + np.floor(365.25*(Y+4716)) + np.floor(30.6001*(M+1)) - 1524.5

    delta_t = 0   # 33.184;
    julian['ephemeris_day'] = (julian['day']) + (delta_t/86400)
    julian['century'] = (julian['day'] - 2451545) / 36525
    julian['ephemeris_century'] = (julian['ephemeris_day'] - 2451545) / 36525
    julian['ephemeris_millenium'] = (julian['ephemeris_century']) / 10

    return julian


def earth_heliocentric_position_calculation(julian):
    """
    % This function compute the earth position relative to the sun, using
    % tabulated values. 
    
    % Tabulated values for the longitude calculation
    % L terms  from the original code.
    """
    # Tabulated values for the longitude calculation
    # L terms  from the original code. 
    L0_terms = np.array([[175347046.0, 0, 0],
                        [3341656.0, 4.6692568, 6283.07585],
                        [34894.0, 4.6261, 12566.1517],
                        [3497.0, 2.7441, 5753.3849],
                        [3418.0, 2.8289, 3.5231],
                        [3136.0, 3.6277, 77713.7715],
                        [2676.0, 4.4181, 7860.4194],
                        [2343.0, 6.1352, 3930.2097],
                        [1324.0, 0.7425, 11506.7698],
                        [1273.0, 2.0371, 529.691],
                        [1199.0, 1.1096, 1577.3435],
                        [990, 5.233, 5884.927],
                        [902, 2.045, 26.298],
                        [857, 3.508, 398.149],
                        [780, 1.179, 5223.694],
                        [753, 2.533, 5507.553],
                        [505, 4.583, 18849.228],
                        [492, 4.205, 775.523],
                        [357, 2.92, 0.067],
                        [317, 5.849, 11790.629],
                        [284, 1.899, 796.298],
                        [271, 0.315, 10977.079],
                        [243, 0.345, 5486.778],
                        [206, 4.806, 2544.314],
                        [205, 1.869, 5573.143],
                        [202, 2.4458, 6069.777],
                        [156, 0.833, 213.299],
                        [132, 3.411, 2942.463],
                        [126, 1.083, 20.775],
                        [115, 0.645, 0.98],
                        [103, 0.636, 4694.003],
                        [102, 0.976, 15720.839],
                        [102, 4.267, 7.114],
                        [99, 6.21, 2146.17],
                        [98, 0.68, 155.42],
                        [86, 5.98, 161000.69],
                        [85, 1.3, 6275.96],
                        [85, 3.67, 71430.7],
                        [80, 1.81, 17260.15],
                        [79, 3.04, 12036.46],
                        [71, 1.76, 5088.63],
                        [74, 3.5, 3154.69],
                        [74, 4.68, 801.82],
                        [70, 0.83, 9437.76],
                        [62, 3.98, 8827.39],
                        [61, 1.82, 7084.9],
                        [57, 2.78, 6286.6],
                        [56, 4.39, 14143.5],
                        [56, 3.47, 6279.55],
                        [52, 0.19, 12139.55],
                        [52, 1.33, 1748.02],
                        [51, 0.28, 5856.48],
                        [49, 0.49, 1194.45],
                        [41, 5.37, 8429.24],
                        [41, 2.4, 19651.05],
                        [39, 6.17, 10447.39],
                        [37, 6.04, 10213.29],
                        [37, 2.57, 1059.38],
                        [36, 1.71, 2352.87],
                        [36, 1.78, 6812.77],
                        [33, 0.59, 17789.85],
                        [30, 0.44, 83996.85],
                        [30, 2.74, 1349.87],
                        [25, 3.16, 4690.48]])

    L1_terms = np.array([[628331966747.0, 0, 0],
                        [206059.0, 2.678235, 6283.07585],
                        [4303.0, 2.6351, 12566.1517],
                        [425.0, 1.59, 3.523],
                        [119.0, 5.796, 26.298],
                        [109.0, 2.966, 1577.344],
                        [93, 2.59, 18849.23],
                        [72, 1.14, 529.69],
                        [68, 1.87, 398.15],
                        [67, 4.41, 5507.55],
                        [59, 2.89, 5223.69],
                        [56, 2.17, 155.42],
                        [45, 0.4, 796.3],
                        [36, 0.47, 775.52],
                        [29, 2.65, 7.11],
                        [21, 5.34, 0.98],
                        [19, 1.85, 5486.78],
                        [19, 4.97, 213.3],
                        [17, 2.99, 6275.96],
                        [16, 0.03, 2544.31],
                        [16, 1.43, 2146.17],
                        [15, 1.21, 10977.08],
                        [12, 2.83, 1748.02],
                        [12, 3.26, 5088.63],
                        [12, 5.27, 1194.45],
                        [12, 2.08, 4694],
                        [11, 0.77, 553.57],
                        [10, 1.3, 3286.6],
                        [10, 4.24, 1349.87],
                        [9, 2.7, 242.73],
                        [9, 5.64, 951.72],
                        [8, 5.3, 2352.87],
                        [6, 2.65, 9437.76],
                        [6, 4.67, 4690.48]])

    L2_terms = np.array([[52919.0, 0, 0],
                        [8720.0, 1.0721, 6283.0758],
                        [309.0, 0.867, 12566.152],
                        [27, 0.05, 3.52],
                        [16, 5.19, 26.3],
                        [16, 3.68, 155.42],
                        [10, 0.76, 18849.23],
                        [9, 2.06, 77713.77],
                        [7, 0.83, 775.52],
                        [5, 4.66, 1577.34],
                        [4, 1.03, 7.11],
                        [4, 3.44, 5573.14],
                        [3, 5.14, 796.3],
                        [3, 6.05, 5507.55],
                        [3, 1.19, 242.73],
                        [3, 6.12, 529.69],
                        [3, 0.31, 398.15],
                        [3, 2.28, 553.57],
                        [2, 4.38, 5223.69],
                        [2, 3.75, 0.98]])

    L3_terms = np.array([[289.0, 5.844, 6283.076],
                        [35, 0, 0],
                        [17, 5.49, 12566.15],
                        [3, 5.2, 155.42],
                        [1, 4.72, 3.52],
                        [1, 5.3, 18849.23],
                        [1, 5.97, 242.73]])
    L4_terms = np.array([[114.0, 3.142, 0],
                        [8, 4.13, 6283.08],
                        [1, 3.84, 12566.15]])

    L5_terms = np.array([1, 3.14, 0])
    L5_terms = np.atleast_2d(L5_terms)    # since L5_terms is 1D, we have to convert it to 2D to avoid indexErrors

    A0 = L0_terms[:, 0]
    B0 = L0_terms[:, 1]
    C0 = L0_terms[:, 2]

    A1 = L1_terms[:, 0]
    B1 = L1_terms[:, 1]
    C1 = L1_terms[:, 2]

    A2 = L2_terms[:, 0]
    B2 = L2_terms[:, 1]
    C2 = L2_terms[:, 2]

    A3 = L3_terms[:, 0]
    B3 = L3_terms[:, 1]
    C3 = L3_terms[:, 2]

    A4 = L4_terms[:, 0]
    B4 = L4_terms[:, 1]
    C4 = L4_terms[:, 2]

    A5 = L5_terms[:, 0]
    B5 = L5_terms[:, 1]
    C5 = L5_terms[:, 2]

    JME = julian['ephemeris_millenium']

    # Compute the Earth Heliochentric longitude from the tabulated values.
    L0 = np.sum(A0 * np.cos(B0 + (C0 * JME)))
    L1 = np.sum(A1 * np.cos(B1 + (C1 * JME)))
    L2 = np.sum(A2 * np.cos(B2 + (C2 * JME)))
    L3 = np.sum(A3 * np.cos(B3 + (C3 * JME)))
    L4 = np.sum(A4 * np.cos(B4 + (C4 * JME)))
    L5 = A5 * np.cos(B5 + (C5 * JME))

    earth_heliocentric_position = dict()
    earth_heliocentric_position['longitude'] = (L0 + (L1 * JME) + (L2 * np.power(JME, 2)) +
                                                          (L3 * np.power(JME, 3)) +
                                                          (L4 * np.power(JME, 4)) +
                                                          (L5 * np.power(JME, 5))) / 1e8
    # Convert the longitude to degrees.
    earth_heliocentric_position['longitude'] = earth_heliocentric_position['longitude'] * 180/np.pi

    # Limit the range to [0,360]
    earth_heliocentric_position['longitude'] = set_to_range(earth_heliocentric_position['longitude'], 0, 360)

    # Tabulated values for the earth heliocentric latitude. 
    # B terms  from the original code. 
    B0_terms = np.array([[280.0, 3.199, 84334.662],
                        [102.0, 5.422, 5507.553],
                        [80, 3.88, 5223.69],
                        [44, 3.7, 2352.87],
                        [32, 4, 1577.34]])

    B1_terms = np.array([[9, 3.9, 5507.55],
                         [6, 1.73, 5223.69]])

    A0 = B0_terms[:, 0]
    B0 = B0_terms[:, 1]
    C0 = B0_terms[:, 2]
    
    A1 = B1_terms[:, 0]
    B1 = B1_terms[:, 1]
    C1 = B1_terms[:, 2]
    
    L0 = np.sum(A0 * np.cos(B0 + (C0 * JME)))
    L1 = np.sum(A1 * np.cos(B1 + (C1 * JME)))

    earth_heliocentric_position['latitude'] = (L0 + (L1 * JME)) / 1e8

    # Convert the latitude to degrees. 
    earth_heliocentric_position['latitude'] = earth_heliocentric_position['latitude'] * 180/np.pi

    # Limit the range to [0,360];
    earth_heliocentric_position['latitude'] = set_to_range(earth_heliocentric_position['latitude'], 0, 360)

    # Tabulated values for radius vector. 
    # R terms from the original code
    R0_terms = np.array([[100013989.0, 0, 0],
                        [1670700.0, 3.0984635, 6283.07585],
                        [13956.0, 3.05525, 12566.1517],
                        [3084.0, 5.1985, 77713.7715],
                        [1628.0, 1.1739, 5753.3849],
                        [1576.0, 2.8469, 7860.4194],
                        [925.0, 5.453, 11506.77],
                        [542.0, 4.564, 3930.21],
                        [472.0, 3.661, 5884.927],
                        [346.0, 0.964, 5507.553],
                        [329.0, 5.9, 5223.694],
                        [307.0, 0.299, 5573.143],
                        [243.0, 4.273, 11790.629],
                        [212.0, 5.847, 1577.344],
                        [186.0, 5.022, 10977.079],
                        [175.0, 3.012, 18849.228],
                        [110.0, 5.055, 5486.778],
                        [98, 0.89, 6069.78],
                        [86, 5.69, 15720.84],
                        [86, 1.27, 161000.69],
                        [85, 0.27, 17260.15],
                        [63, 0.92, 529.69],
                        [57, 2.01, 83996.85],
                        [56, 5.24, 71430.7],
                        [49, 3.25, 2544.31],
                        [47, 2.58, 775.52],
                        [45, 5.54, 9437.76],
                        [43, 6.01, 6275.96],
                        [39, 5.36, 4694],
                        [38, 2.39, 8827.39],
                        [37, 0.83, 19651.05],
                        [37, 4.9, 12139.55],
                        [36, 1.67, 12036.46],
                        [35, 1.84, 2942.46],
                        [33, 0.24, 7084.9],
                        [32, 0.18, 5088.63],
                        [32, 1.78, 398.15],
                        [28, 1.21, 6286.6],
                        [28, 1.9, 6279.55],
                        [26, 4.59, 10447.39]])

    R1_terms = np.array([[103019.0, 1.10749, 6283.07585],
                        [1721.0, 1.0644, 12566.1517],
                        [702.0, 3.142, 0],
                        [32, 1.02, 18849.23],
                        [31, 2.84, 5507.55],
                        [25, 1.32, 5223.69],
                        [18, 1.42, 1577.34],
                        [10, 5.91, 10977.08],
                        [9, 1.42, 6275.96],
                        [9, 0.27, 5486.78]])

    R2_terms = np.array([[4359.0, 5.7846, 6283.0758],
                        [124.0, 5.579, 12566.152],
                        [12, 3.14, 0],
                        [9, 3.63, 77713.77],
                        [6, 1.87, 5573.14],
                        [3, 5.47, 18849]])

    R3_terms = np.array([[145.0, 4.273, 6283.076],
                        [7, 3.92, 12566.15]])
    
    R4_terms = [4, 2.56, 6283.08]
    R4_terms = np.atleast_2d(R4_terms)    # since L5_terms is 1D, we have to convert it to 2D to avoid indexErrors

    A0 = R0_terms[:, 0]
    B0 = R0_terms[:, 1]
    C0 = R0_terms[:, 2]
    
    A1 = R1_terms[:, 0]
    B1 = R1_terms[:, 1]
    C1 = R1_terms[:, 2]
    
    A2 = R2_terms[:, 0]
    B2 = R2_terms[:, 1]
    C2 = R2_terms[:, 2]
    
    A3 = R3_terms[:, 0]
    B3 = R3_terms[:, 1]
    C3 = R3_terms[:, 2]
    
    A4 = R4_terms[:, 0]
    B4 = R4_terms[:, 1]
    C4 = R4_terms[:, 2]

    # Compute the Earth heliocentric radius vector
    L0 = np.sum(A0 * np.cos(B0 + (C0 * JME)))
    L1 = np.sum(A1 * np.cos(B1 + (C1 * JME)))
    L2 = np.sum(A2 * np.cos(B2 + (C2 * JME)))
    L3 = np.sum(A3 * np.cos(B3 + (C3 * JME)))
    L4 = A4 * np.cos(B4 + (C4 * JME))

    # Units are in AU
    earth_heliocentric_position['radius'] = (L0 + (L1 * JME) + (L2 * np.power(JME, 2)) +
                                             (L3 * np.power(JME, 3)) +
                                             (L4 * np.power(JME, 4))) / 1e8

    return earth_heliocentric_position


def sun_geocentric_position_calculation(earth_heliocentric_position):
    """
    % This function compute the sun position relative to the earth.
    """
    sun_geocentric_position = dict()
    sun_geocentric_position['longitude'] = earth_heliocentric_position['longitude'] + 180
    # Limit the range to [0,360];
    sun_geocentric_position['longitude'] = set_to_range(sun_geocentric_position['longitude'], 0, 360)

    sun_geocentric_position['latitude'] = -earth_heliocentric_position['latitude']
    # Limit the range to [0,360]
    sun_geocentric_position['latitude'] = set_to_range(sun_geocentric_position['latitude'], 0, 360)
    return sun_geocentric_position


def nutation_calculation(julian):
    """
    % This function compute the nutation in longtitude and in obliquity, in
    % degrees.
    :param julian:
    :return: nutation
    """

    # All Xi are in degrees.
    JCE = julian['ephemeris_century']

    # 1. Mean elongation of the moon from the sun
    p = np.atleast_2d([(1/189474), -0.0019142, 445267.11148, 297.85036])

    # X0 = polyval(p, JCE);
    X0 = p[0, 0] * np.power(JCE, 3) + p[0, 1] * np.power(JCE, 2) + p[0, 2] * JCE + p[0, 3]   # This is faster than polyval...

    # 2. Mean anomaly of the sun (earth)
    p = np.atleast_2d([-(1/300000), -0.0001603, 35999.05034, 357.52772])

    # X1 = polyval(p, JCE)
    X1 = p[0, 0] * np.power(JCE, 3) + p[0, 1] * np.power(JCE, 2) + p[0, 2] * JCE + p[0, 3]

    # 3. Mean anomaly of the moon
    p = np.atleast_2d([(1/56250), 0.0086972, 477198.867398, 134.96298])
    
    # X2 = polyval(p, JCE);
    X2 = p[0, 0] * np.power(JCE, 3) + p[0, 1] * np.power(JCE, 2) + p[0, 2] * JCE + p[0, 3]

    # 4. Moon argument of latitude
    p = np.atleast_2d([(1/327270), -0.0036825, 483202.017538, 93.27191])

    # X3 = polyval(p, JCE)
    X3 = p[0, 0] * np.power(JCE, 3) + p[0, 1] * np.power(JCE, 2) + p[0, 2] * JCE + p[0, 3]

    # 5. Longitude of the ascending node of the moon's mean orbit on the
    # ecliptic, measured from the mean equinox of the date
    p = np.atleast_2d([(1/450000), 0.0020708, -1934.136261, 125.04452])

    # X4 = polyval(p, JCE);
    X4 = p[0, 0] * np.power(JCE, 3) + p[0, 1] * np.power(JCE, 2) + p[0, 2] * JCE + p[0, 3]

    # Y tabulated terms from the original code
    Y_terms = np.array([[0, 0, 0, 0, 1],
                        [-2, 0, 0, 2, 2],
                        [0, 0, 0, 2, 2],
                        [0, 0, 0, 0, 2],
                        [0, 1, 0, 0, 0],
                        [0, 0, 1, 0, 0],
                        [-2, 1, 0, 2, 2],
                        [0, 0, 0, 2, 1],
                        [0, 0, 1, 2, 2],
                        [-2, -1, 0, 2, 2],
                        [-2, 0, 1, 0, 0],
                        [-2, 0, 0, 2, 1],
                        [0, 0, -1, 2, 2],
                        [2, 0, 0, 0, 0],
                        [0, 0, 1, 0, 1],
                        [2, 0, -1, 2, 2],
                        [0, 0, -1, 0, 1],
                        [0, 0, 1, 2, 1],
                        [-2, 0, 2, 0, 0],
                        [0, 0, -2, 2, 1],
                        [2, 0, 0, 2, 2],
                        [0, 0, 2, 2, 2],
                        [0, 0, 2, 0, 0],
                        [-2, 0, 1, 2, 2],
                        [0, 0, 0, 2, 0],
                        [-2, 0, 0, 2, 0],
                        [0, 0, -1, 2, 1],
                        [0, 2, 0, 0, 0],
                        [2, 0, -1, 0, 1],
                        [-2, 2, 0, 2, 2],
                        [0, 1, 0, 0, 1],
                        [-2, 0, 1, 0, 1],
                        [0, -1, 0, 0, 1],
                        [0, 0, 2, -2, 0],
                        [2, 0, -1, 2, 1],
                        [2, 0, 1, 2, 2],
                        [0, 1, 0, 2, 2],
                        [-2, 1, 1, 0, 0],
                        [0, -1, 0, 2, 2],
                        [2, 0, 0, 2, 1],
                        [2, 0, 1, 0, 0],
                        [-2, 0, 2, 2, 2],
                        [-2, 0, 1, 2, 1],
                        [2, 0, -2, 0, 1],
                        [2, 0, 0, 0, 1],
                        [0, -1, 1, 0, 0],
                        [-2, -1, 0, 2, 1],
                        [-2, 0, 0, 0, 1],
                        [0, 0, 2, 2, 1],
                        [-2, 0, 2, 0, 1],
                        [-2, 1, 0, 2, 1],
                        [0, 0, 1, -2, 0],
                        [-1, 0, 1, 0, 0],
                        [-2, 1, 0, 0, 0],
                        [1, 0, 0, 0, 0],
                        [0, 0, 1, 2, 0],
                        [0, 0, -2, 2, 2],
                        [-1, -1, 1, 0, 0],
                        [0, 1, 1, 0, 0],
                        [0, -1, 1, 2, 2],
                        [2, -1, -1, 2, 2],
                        [0, 0, 3, 2, 2],
                        [2, -1, 0, 2, 2]])

    nutation_terms = np.array([[-171996, -174.2, 92025, 8.9],
                                [-13187, -1.6, 5736, -3.1],
                                [-2274, -0.2, 977, -0.5],
                                [2062, 0.2, -895, 0.5],
                                [1426, -3.4, 54, -0.1],
                                [712, 0.1, -7, 0],
                                [-517, 1.2, 224, -0.6],
                                [-386, -0.4, 200, 0],
                                [-301, 0, 129, -0.1],
                                [217, -0.5, -95, 0.3],
                                [-158, 0, 0, 0],
                                [129, 0.1, -70, 0],
                                [123, 0, -53, 0],
                                [63, 0, 0, 0],
                                [63, 0.1, -33, 0],
                                [-59, 0, 26, 0],
                                [-58, -0.1, 32, 0],
                                [-51, 0, 27, 0],
                                [48, 0, 0, 0],
                                [46, 0, -24, 0],
                                [-38, 0, 16, 0],
                                [-31, 0, 13, 0],
                                [29, 0, 0, 0],
                                [29, 0, -12, 0],
                                [26, 0, 0, 0],
                                [-22, 0, 0, 0],
                                [21, 0, -10, 0],
                                [17, -0.1, 0, 0],
                                [16, 0, -8, 0],
                                [-16, 0.1, 7, 0],
                                [-15, 0, 9, 0],
                                [-13, 0, 7, 0],
                                [-12, 0, 6, 0],
                                [11, 0, 0, 0],
                                [-10, 0, 5, 0],
                                [-8, 0, 3, 0],
                                [7, 0, -3, 0],
                                [-7, 0, 0, 0],
                                [-7, 0, 3, 0],
                                [-7, 0, 3, 0],
                                [6, 0, 0, 0],
                                [6, 0, -3, 0],
                                [6, 0, -3, 0],
                                [-6, 0, 3, 0],
                                [-6, 0, 3, 0],
                                [5, 0, 0, 0],
                                [-5, 0, 3, 0],
                                [-5, 0, 3, 0],
                                [-5, 0, 3, 0],
                                [4, 0, 0, 0],
                                [4, 0, 0, 0],
                                [4, 0, 0, 0],
                                [-4, 0, 0, 0],
                                [-4, 0, 0, 0],
                                [-4, 0, 0, 0],
                                [3, 0, 0, 0],
                                [-3, 0, 0, 0],
                                [-3, 0, 0, 0],
                                [-3, 0, 0, 0],
                                [-3, 0, 0, 0],
                                [-3, 0, 0, 0],
                                [-3, 0, 0, 0],
                                [-3, 0, 0, 0]])

    # Using the tabulated values, compute the delta_longitude and
    # delta_obliquity.
    Xi = np.array([X0, X1, X2, X3, X4])    # a col mat in octave

    tabulated_argument = Y_terms.dot(np.transpose(Xi)) * (np.pi/180)

    delta_longitude = (nutation_terms[:, 0] + (nutation_terms[:, 1] * JCE)) * np.sin(tabulated_argument)
    delta_obliquity = (nutation_terms[:, 2] + (nutation_terms[:, 3] * JCE)) * np.cos(tabulated_argument)

    nutation = dict()    # init nutation dictionary
    # Nutation in longitude
    nutation['longitude'] = np.sum(delta_longitude) / 36000000

    # Nutation in obliquity
    nutation['obliquity'] = np.sum(delta_obliquity) / 36000000

    return nutation


def true_obliquity_calculation(julian, nutation):
    """
    This function compute the true obliquity of the ecliptic.

    :param julian:
    :param nutation:
    :return:
    """

    p = np.atleast_2d([2.45, 5.79, 27.87, 7.12, -39.05, -249.67, -51.38, 1999.25, -1.55, -4680.93, 84381.448])

    # mean_obliquity = polyval(p, julian.ephemeris_millenium/10);
    U = julian['ephemeris_millenium'] / 10
    mean_obliquity = p[0, 0] * np.power(U, 10) + p[0, 1] * np.power(U, 9) + \
                     p[0, 2] * np.power(U, 8) + p[0, 3] * np.power(U, 7) + \
                     p[0, 4] * np.power(U, 6) + p[0, 5] * np.power(U, 5) + \
                     p[0, 6] * np.power(U, 4) + p[0, 7] * np.power(U, 3) + \
                     p[0, 8] * np.power(U, 2) + p[0, 9] * U + p[0, 10]

    true_obliquity = (mean_obliquity/3600) + nutation['obliquity']

    return true_obliquity


def abberation_correction_calculation(earth_heliocentric_position):
    """
    This function compute the aberration_correction, as a function of the
    earth-sun distance.

    :param earth_heliocentric_position:
    :return:
    """
    aberration_correction = -20.4898/(3600*earth_heliocentric_position['radius'])
    return aberration_correction


def apparent_sun_longitude_calculation(sun_geocentric_position, nutation, aberration_correction):
    """
    This function compute the sun apparent longitude

    :param sun_geocentric_position:
    :param nutation:
    :param aberration_correction:
    :return:
    """
    apparent_sun_longitude = sun_geocentric_position['longitude'] + nutation['longitude'] + aberration_correction
    return apparent_sun_longitude


def apparent_stime_at_greenwich_calculation(julian, nutation, true_obliquity):
    """
    This function compute the apparent sideral time at Greenwich.

    :param julian:
    :param nutation:
    :param true_obliquity:
    :return:
    """

    JD = julian['day']
    JC = julian['century']

    # Mean sideral time, in degrees
    mean_stime = 280.46061837 + (360.98564736629*(JD-2451545)) + \
                 (0.000387933*np.power(JC, 2)) - \
                 (np.power(JC, 3)/38710000)

    # Limit the range to [0-360];
    mean_stime = set_to_range(mean_stime, 0, 360)

    apparent_stime_at_greenwich = mean_stime + (nutation['longitude'] * np.cos(true_obliquity * np.pi/180))
    return apparent_stime_at_greenwich


def sun_rigth_ascension_calculation(apparent_sun_longitude, true_obliquity, sun_geocentric_position):
    """
    This function compute the sun rigth ascension.
    :param apparent_sun_longitude:
    :param true_obliquity:
    :param sun_geocentric_position:
    :return:
    """

    argument_numerator = (np.sin(apparent_sun_longitude * np.pi/180) * np.cos(true_obliquity * np.pi/180)) - \
        (np.tan(sun_geocentric_position['latitude'] * np.pi/180) * np.sin(true_obliquity * np.pi/180))
    argument_denominator = np.cos(apparent_sun_longitude * np.pi/180);

    sun_rigth_ascension = np.arctan2(argument_numerator, argument_denominator) * 180/np.pi
    # Limit the range to [0,360];
    sun_rigth_ascension = set_to_range(sun_rigth_ascension, 0, 360)
    return sun_rigth_ascension


def sun_geocentric_declination_calculation(apparent_sun_longitude, true_obliquity, sun_geocentric_position):
    """

    :param apparent_sun_longitude:
    :param true_obliquity:
    :param sun_geocentric_position:
    :return:
    """

    argument = (np.sin(sun_geocentric_position['latitude'] * np.pi/180) * np.cos(true_obliquity * np.pi/180)) + \
        (np.cos(sun_geocentric_position['latitude'] * np.pi/180) * np.sin(true_obliquity * np.pi/180) *
         np.sin(apparent_sun_longitude * np.pi/180))

    sun_geocentric_declination = np.arcsin(argument) * 180/np.pi
    return sun_geocentric_declination


def observer_local_hour_calculation(apparent_stime_at_greenwich, location, sun_rigth_ascension):
    """
    This function computes observer local hour.

    :param apparent_stime_at_greenwich:
    :param location:
    :param sun_rigth_ascension:
    :return:
    """

    observer_local_hour = apparent_stime_at_greenwich + location['longitude'] - sun_rigth_ascension
    # Set the range to [0-360]
    observer_local_hour = set_to_range(observer_local_hour, 0, 360)
    return observer_local_hour


def topocentric_sun_position_calculate(earth_heliocentric_position, location,
                                       observer_local_hour, sun_rigth_ascension, sun_geocentric_declination):
    """
    This function compute the sun position (rigth ascension and declination)
    with respect to the observer local position at the Earth surface.

    :param earth_heliocentric_position:
    :param location:
    :param observer_local_hour:
    :param sun_rigth_ascension:
    :param sun_geocentric_declination:
    :return:
    """

    # Equatorial horizontal parallax of the sun in degrees
    eq_horizontal_parallax = 8.794 / (3600 * earth_heliocentric_position['radius'])

    # Term u, used in the following calculations (in radians)
    u = np.arctan(0.99664719 * np.tan(location['latitude'] * np.pi/180))

    # Term x, used in the following calculations
    x = np.cos(u) + ((location['altitude']/6378140) * np.cos(location['latitude'] * np.pi/180))

    # Term y, used in the following calculations
    y = (0.99664719 * np.sin(u)) + ((location['altitude']/6378140) * np.sin(location['latitude'] * np.pi/180))

    # Parallax in the sun rigth ascension (in radians)
    nominator = -x * np.sin(eq_horizontal_parallax * np.pi/180) * np.sin(observer_local_hour * np.pi/180)
    denominator = np.cos(sun_geocentric_declination * np.pi/180) - (x * np.sin(eq_horizontal_parallax * np.pi/180) *
                                                                    np.cos(observer_local_hour * np.pi/180))
    sun_rigth_ascension_parallax = np.arctan2(nominator, denominator)
    # Conversion to degrees.
    topocentric_sun_position = dict()
    topocentric_sun_position['rigth_ascension_parallax'] = sun_rigth_ascension_parallax * 180/np.pi

    # Topocentric sun rigth ascension (in degrees)
    topocentric_sun_position['rigth_ascension'] = sun_rigth_ascension + (sun_rigth_ascension_parallax * 180/np.pi)

    # Topocentric sun declination (in degrees)
    nominator = (np.sin(sun_geocentric_declination * np.pi/180) - (y*np.sin(eq_horizontal_parallax * np.pi/180))) * \
                np.cos(sun_rigth_ascension_parallax)
    denominator = np.cos(sun_geocentric_declination * np.pi/180) - (y*np.sin(eq_horizontal_parallax * np.pi/180)) * \
                                                                   np.cos(observer_local_hour * np.pi/180)
    topocentric_sun_position['declination'] = np.arctan2(nominator, denominator) * 180/np.pi
    return topocentric_sun_position


def topocentric_local_hour_calculate(observer_local_hour, topocentric_sun_position):
    """
    This function compute the topocentric local jour angle in degrees

    :param observer_local_hour:
    :param topocentric_sun_position:
    :return:
    """

    topocentric_local_hour = observer_local_hour - topocentric_sun_position['rigth_ascension_parallax']
    return topocentric_local_hour


def sun_topocentric_zenith_angle_calculate(location, topocentric_sun_position, topocentric_local_hour):
    """
    This function compute the sun zenith angle, taking into account the
    atmospheric refraction. A default temperature of 283K and a
    default pressure of 1010 mbar are used.

    :param location:
    :param topocentric_sun_position:
    :param topocentric_local_hour:
    :return:
    """

    # Topocentric elevation, without atmospheric refraction
    argument = (np.sin(location['latitude'] * np.pi/180) * np.sin(topocentric_sun_position['declination'] * np.pi/180)) + \
    (np.cos(location['latitude'] * np.pi/180) * np.cos(topocentric_sun_position['declination'] * np.pi/180) *
     np.cos(topocentric_local_hour * np.pi/180))
    true_elevation = np.arcsin(argument) * 180/np.pi

    # Atmospheric refraction correction (in degrees)
    argument = true_elevation + (10.3/(true_elevation + 5.11))
    refraction_corr = 1.02 / (60 * np.tan(argument * np.pi/180))

    # For exact pressure and temperature correction, use this,
    # with P the pressure in mbar amd T the temperature in Kelvins:
    # refraction_corr = (P/1010) * (283/T) * 1.02 / (60 * tan(argument * pi/180));

    # Apparent elevation
    apparent_elevation = true_elevation + refraction_corr

    sun = dict()
    sun['zenith'] = 90 - apparent_elevation

    # Topocentric azimuth angle. The +180 conversion is to pass from astronomer
    # notation (westward from south) to navigation notation (eastward from
    # north);
    nominator = np.sin(topocentric_local_hour * np.pi/180)
    denominator = (np.cos(topocentric_local_hour * np.pi/180) * np.sin(location['latitude'] * np.pi/180)) - \
    (np.tan(topocentric_sun_position['declination'] * np.pi/180) * np.cos(location['latitude'] * np.pi/180))
    sun['azimuth'] = (np.arctan2(nominator, denominator) * 180/np.pi) + 180

    # Set the range to [0-360]
    sun['azimuth'] = set_to_range(sun['azimuth'], 0, 360)
    return sun


def set_to_range(var, min_interval, max_interval):
    """
    Sets a variable in range min_interval and max_interval

    :param var:
    :param min_interval:
    :param max_interval:
    :return:
    """
    var = var - max_interval * np.floor(var/max_interval)

    if var < min_interval:
        var = var + max_interval
    return var
