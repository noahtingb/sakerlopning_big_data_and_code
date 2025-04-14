"""Adjusted code by noahtingb, orginally based on biglimp:s PETCalculatorWEB"""
import numpy as np
import datetime

def sun_position(time, location):
    def julian_calculation(t_input):
        time = t_input

        if time['month'] == 1 or time['month'] == 2:
            Y = time['year'] - 1
            M = time['month'] + 12
        else:
            Y = time['year']
            M = time['month']

        ut_time = ((time['hour'] - time['UTC'])/24) + (time['min']/(60*24)) + (time['sec']/(60*60*24))   # time of day in UT time.
        D = time['day'] + ut_time   # Day of month in decimal time, ex. 2sd day of month at 12:30:30UT, D=2.521180556

        #year after 1582
        julian = dict()
        julian['day'] = D + (2 - np.floor(Y/100) + np.floor(np.floor(Y/100)/4)) + np.floor(365.25*(Y+4716)) + np.floor(30.6001*(M+1)) - 1524.5
        julian['century'] = (julian['day'] - 2451545) / 36525
        julian['millenium'] = (julian['century']) / 10

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
        JME = julian['millenium']

        L_terms = [np.array([[175347046.0, 0, 0],
                            [3341656.0, 4.6692568, 6283.07585],
                            [34894.0, 4.6261, 12566.1517]]),np.array(
                            [[628331966747.0, 0, 0],
                            [206059.0, 2.678235, 6283.07585]]),
                            np.array([[52919.0, 0, 0]])]
        Lsum=0
        for i in range(3):
            A = L_terms[i][:, 0]
            B = L_terms[i][:, 1]
            C = L_terms[i][:, 2]
            Lsum+=np.sum(A * np.cos(B + (C * JME)))*JME**i
        # Compute the Earth Heliochentric longitude from the tabulated values.

        earth_heliocentric_position = dict()
        earth_heliocentric_position['longitude'] = Lsum*1e-8* 180/np.pi
        # Convert the longitude to degrees.

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
        
        L0 = np.sum(B0_terms[:, 0] * np.cos(B0_terms[:, 1] + (B0_terms[:, 2] * JME)))
        L1 = np.sum(B1_terms[:, 0] * np.cos(B1_terms[:, 1] + (B1_terms[:, 2] * JME)))

        earth_heliocentric_position['latitude'] = (L0+L1 * JME) / 1e8* 180/np.pi

        # Convert the latitude to degrees. 

        earth_heliocentric_position['latitude'] = set_to_range(earth_heliocentric_position['latitude'], 0, 360)
        # Tabulated values for radius vector. 
        # R terms from the original code

        R_terms = [np.array([[100013989.0, 0, 0],
                            [1670700.0, 3.0984635, 6283.07585],
                            [13956.0, 3.05525, 12566.1517]]),
                            np.array([[103019.0, 1.10749, 6283.07585],
                            [1721.0, 1.0644, 12566.1517]]),
                            np.array([[4359.0, 5.7846, 6283.0758]])]
        Rsum=0
        for i in range(3):
            A = R_terms[i][:, 0]
            B = R_terms[i][:, 1]
            C = R_terms[i][:, 2]
            Rsum+=np.sum(A * np.cos(B + (C * JME)))*JME**i
    
        # Units are in AU
        earth_heliocentric_position['radius'] = Rsum*1e8
        return earth_heliocentric_position


    def nutation_calculation(julian):
        """
        % This function compute the nutation in longtitude and in obliquity, in
        % degrees.
        :param julian:
        :return: nutation
        """
        # All Xi are in degrees.
        JCE = julian['century']
        #print("JCE",JCE)

        # Y tabulated terms from the original code
        Y_terms = np.array([[ 0,  0, 0, 0, 1],
                            [-2,  0, 0, 2, 2],
                            [ 0,  0, 0, 2, 2],
                            [ 0,  0, 0, 0, 2],
                            [ 0,  1, 0, 0, 0],
                            [ 0,  0, 1, 0, 0],
                            [-2,  1, 0, 2, 2],
                            [ 0,  0, 0, 2, 1],
                            [ 0,  0, 1, 2, 2],
                            [-2, -1, 0, 2, 2],
                            [-2,  0, 1, 0, 0],
                            [-2,  0, 0, 2, 1],
                            [ 0,  0,-1, 2, 2],
                            [ 2,  0, 0, 0, 0],
                            [ 0,  0, 1, 0, 1]])

        n0=np.array([-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217,-158, 129, 123, 63, 63])
        n1=np.array([-174.2, -1.6, -.2, .2, -3.4, .1, 1.2, .4, 0, -.5, 0, .1, 0, 0, .1])
        n2=np.array([92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, 0, -70, -53, 0, -33])
        n3=np.array([8.9, -3.1, -.5, .5, -.1, 0, -.6, 0, -.1, .3, 0, 0, 0, 0, 0])

        # Using the tabulated values, compute the delta_longitude and
        # delta_obliquity.

        # 1. Mean elongation of the moon from the sun
        X0 =  445267.11148 * JCE + 297.85036   # This is faster than polyval...

        # 2. Mean anomaly of the sun (earth)
        X1 =35999.05034 * JCE + 357.52772

        # 3. Mean anomaly of the moon    
        X2 = 477198.867398 * JCE + 134.96298

        # 4. Moon argument of latitude
        X3 = 483202.017538 * JCE + 93.27191

        # 5. Longitude of the ascending node of the moon's mean orbit on the
        # ecliptic, measured from the mean equinox of the date
        X4 = -1934.136261 * JCE + 125.04452

        Xi = np.array([X0, X1, X2, X3, X4])    # a col mat in octave

        tabulated_argument = Y_terms.dot(np.transpose(Xi)) * (np.pi/180)
        delta_longitude = np.sum((n0 + n1 * JCE) * np.sin(tabulated_argument))
        delta_obliquity = np.sum((n2 + n3* JCE) * np.cos(tabulated_argument))

        # Nutation in longitude
        # Nutation in obliquity
        nutation={'longitude': delta_longitude / 36000000,'obliquity':delta_obliquity / 36000000}
        #print(nutation['longitude'],nutation['obliquity'])
        return nutation

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

    def topocentric_sun_position_calculate(earth_heliocentric_position, location,
                                        observer_local_hour, sun_geocentric_declination):
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
    
        # Topocentric sun declination (in degrees)
        nominator = (np.sin(sun_geocentric_declination * np.pi/180) - (y*np.sin(eq_horizontal_parallax * np.pi/180))) * \
                    np.cos(sun_rigth_ascension_parallax)
        denominator = np.cos(sun_geocentric_declination * np.pi/180) - (y*np.sin(eq_horizontal_parallax * np.pi/180)) * \
                                                                    np.cos(observer_local_hour * np.pi/180)
        topocentric_sun_position['declination'] = np.arctan2(nominator, denominator) * 180/np.pi
        return topocentric_sun_position

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
    # 1. Calculate the Julian Day, and Century. Julian Ephemeris day, century
    # and millenium are calculated using a mean delta_t of 33.184 seconds.
    julian = julian_calculation(time)

    # 2. Calculate the Earth heliocentric longitude, latitude, and radius
    earth_heliocentric_position = earth_heliocentric_position_calculation(julian)

    # 3. Calculate the geocentric longitude and latitude
    sun_geocentric_position = {'longitude': set_to_range(earth_heliocentric_position['longitude'] + 180, 0, 360),\
                               'latitude':set_to_range(-earth_heliocentric_position['latitude'], 0, 360)}

    # 4. Calculate the nutation in longitude and obliquity (in degrees).
    nutation = nutation_calculation(julian)

    # 5. Calculate the true obliquity of the ecliptic (in degrees).
    U = float(julian['millenium'] / 10)
    true_obliquity = nutation['obliquity']+1/3600*(-249.67 * U**5 -51.38 * U**4 + 1999.25 * U**3 -1.55 * U**2 -4680.93 * U + 84381.448) 
    print(nutation['obliquity'],1/3600*(-249.67 * U**5 -51.38 * U**4 + 1999.25 * U**3 -1.55 * U**2 -4680.93 * U + 84381.448))
    # 6. Calculate the aberration correction (in degrees)
    aberration_correction = aberration_correction = -20.4898/(3600*earth_heliocentric_position['radius'])

    # 7. Calculate the apparent sun longitude in degrees)
    apparent_sun_longitude = sun_geocentric_position['longitude'] + nutation['longitude'] + aberration_correction
    # 8. Calculate the apparent sideral time at Greenwich (in degrees)
    mean_stime = 280.46061837 + (360.98564736629*(julian['day']-2451545)) + (0.000387933*np.power(julian['century'], 2)) - (np.power(julian['century'], 3)/38710000)
    apparent_stime_at_greenwich = set_to_range(mean_stime, 0, 360) + (nutation['longitude'] * np.cos(true_obliquity * np.pi/180))

    # 9. Calculate the sun rigth ascension (in degrees)
    sun_rigth_ascension = sun_rigth_ascension_calculation(apparent_sun_longitude, true_obliquity, sun_geocentric_position)

    # 10. Calculate the geocentric sun declination (in degrees). Positive or
    # negative if the sun is north or south of the celestial equator.
    sun_geocentric_declination = np.arcsin((np.sin(sun_geocentric_position['latitude'] * np.pi/180) * np.cos(true_obliquity * np.pi/180)) + 
                                           (np.cos(sun_geocentric_position['latitude'] * np.pi/180) * np.sin(true_obliquity * np.pi/180) *
                                            np.sin(apparent_sun_longitude * np.pi/180)))* 180/np.pi

    # 11. Calculate the observer local hour angle (in degrees, westward from south).
    observer_local_hour = set_to_range(apparent_stime_at_greenwich + location['longitude'] - sun_rigth_ascension, 0, 360)

    # 12. Calculate the topocentric sun position (rigth ascension, declination and
    # rigth ascension parallax in degrees)
    topocentric_sun_position = topocentric_sun_position_calculate(earth_heliocentric_position, location, observer_local_hour, sun_geocentric_declination)

    # 13. Calculate the topocentric local hour angle (in degrees)
    topocentric_local_hour = observer_local_hour - topocentric_sun_position['rigth_ascension_parallax']
        
    # 14. Calculate the topocentric zenith and azimuth angle (in degrees)
    sun = sun_topocentric_zenith_angle_calculate(location, topocentric_sun_position, topocentric_local_hour)

    return sun


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

    p = 1013.  # Pressure in millibars
    Itoa = 1370.0  # Effective solar constant
    b = 2.*np.pi*jday/365.
    D= np.sqrt((1.00011+np.dot(0.034221, np.cos(b))+np.dot(0.001280, np.sin(b))+np.dot(0.000719,np.cos((2.*b)))+np.dot(0.000077, np.sin((2.*b))))) # irradiance differences due to Sun-Earth distances
    m = 35. * np.cos(zen) * ((1224. * (np.cos(zen)**2) + 1) ** (-1/2.))     # optical air mass at p=1013
    Trpg = 1.021-0.084*(m*(0.000949*p+0.051))**0.5  # Transmission coefficient for Rayliegh scattering and permanent gases

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
    # b=I0==abs(zen)>np.pi/2
    # I0(b==1)=0
    # clear b;
    if not(np.isreal(I0)):
        I0 = 0

    corr=0.1473*np.log(90-(zen/np.pi*180))+0.3454  # 20070329

    CIuncorr = radG / I0
    CI = CIuncorr + (1-corr)
    I0et = Itoa*np.cos(zen)*D  # extra terrestial solar radiation
    Kt = radG / I0et
    if np.isnan(CI):
        CI = float('Inf')

    return I0, CI, Kt

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
    if altitude < 0:
        print("Airmass")
        print(AirMass)
        print(PerezBrightness)
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

def Lside_veg_v2020a(azimuth,altitude,Ta,Tw,SBC,Ldown,esky,t,F_sh,CI,Lup):
    
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
    Lrefl = (Ldown + Lup) * (viktrefl) * (1 - 0.9) * 0.5
    Lground = Lup * 0.5
    Lveg = SBC * 0.9 * ((Ta + 273.15) ** 4) * viktveg * 0.5

    if altitude > 0:  # daytime
        # betasun = np.arctan(0.5*np.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
        if (azimuth > (180-t))  and  (azimuth <= (360-t)):
            Lwallsun=SBC*0.9*((Ta+273.15+Tw*np.sin(aziE*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*0.5
    
    Lsum = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    if altitude>0: # daytime
        if (azimuth <= (90-t))  or  (azimuth > (270-t)):
            Lwallsun=SBC*0.9*((Ta+273.15+Tw*np.sin(aziS*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum+= Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        if (azimuth > (360-t))  or  (azimuth <= (180-t)):
            Lwallsun=SBC*0.9*((Ta+273.15+Tw*np.sin(aziW*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum += Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        if (azimuth > (90-t))  and  (azimuth <= (270-t)):
            Lwallsun=SBC*0.9*((Ta+273.15+Tw*np.sin(aziN*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*0.9*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum += Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        return Lsum
    return Lsum*4

def Kside_veg_v2019a(radI, radD, radG,altitude, F_sh, Kup, lv, diffsh):
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
    Keast = (0.2 * (svfviktbuveg * (radG * (1 - F_sh) + radD * F_sh)) + Kup) * 0.5
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

def daylen(doy, xlat):
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
    day_len = 12.0 + 24.0*np.arcsin(SOC)/np.pi
    sun_upp = 12.0 - day_len/2.0 
    sun_down = 12.0 + day_len/2.0
    return day_len, DEC, sun_down, sun_upp

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
    altmax = sunmaximum

    YMDHM = YMD + datetime.timedelta(hours=hour) + datetime.timedelta(minutes=minu) - datetime.timedelta(days=halftimestepdec)

    time = {"year":YMDHM.year,'month': YMDHM.month,'day':YMDHM.day,'hour': YMDHM.hour,'min':YMDHM.minute,'sec':0,"UTC":UTC}
    sun = sun_position(time, location)
    altitude = 90. - sun['zenith']
    azimuth = sun['azimuth']
    zen = sun['zenith'] * (3.141592653589793238/180.)
        # day of year and check for leap year
    if (time['year'] % 4) == 0 and ( ((time['year'] % 100)==0 and (time['year'] % 400) == 0) or ((time['year'] % 100)!=0)):
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    jday = sum(dayspermonth[0:time['month'] - 1]) + time['day']

    return None, altitude, azimuth, zen, jday, None, dectime, altmax

def Solweig1D_2020a_calc(Fside, Fup, Fcyl,
                        location, Ta, RH, 
                        year, month, day, hour, minu):
    RH*= 1e-2
    # Meteorological data, Should maybe be move out later on.
    #day of year
    if (year % 4) == 0 and ( ((year % 100)==0 and (year % 400) == 0) or ((year % 100)!=0)):
        dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy= sum(dayspermonth[0:month - 1]) + day

    _, altitude, azimuth, zen, jday, _, dectime, altmax = Solweig_2015a_metdata_noload(year,doy,hour,minu, location, 0)

    if altitude > 0.:
        I0, _, _= clearnessindex_2013b(zen, jday, Ta, RH, 40, location)
        radG = I0
    else:
        radG = 0.

    if (year % 4) == 0 and ( ((year % 100)==0 and (year % 400) == 0) or ((year % 100)!=0)):
        dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy=sum(dayspermonth[0:month - 1]) + day

    svfalfa = np.arcsin(np.exp((np.log((1.-0.6))/2.)))


    skyvaultalt = np.atleast_2d([])
    skyvaultaltint = [6, 18, 30, 42, 54, 66, 78]
    skyvaultaziint = [12, 12, 15, 15, 20, 30, 60]
    for j in range(7):
        for k in range(1, int(360/skyvaultaziint[j]) + 1):
            skyvaultalt = np.append(skyvaultalt, skyvaultaltint[j])

    skyvaultalt = np.append(skyvaultalt, 90)

    diffsh = np.zeros((145))
    svfalfadeg = svfalfa / (np.pi / 180.)
    for k in range(0, 145):
        if skyvaultalt[k] > svfalfadeg:
            diffsh[k] = 1
    # This is the core function of the SOLWEIG1D model, 2019-Jun-21
    # Fredrik Lindberg, fredrikl@gvc.gu.se, Goteborg Urban Climate Group, Gothenburg University, Sweden
    # Instrument offset in degrees
    t = 0.
    # Stefan Bolzmans Constant
    SBC = 5.67051e-8
    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(jday, location['latitude'])
    # Vapor pressure
    ea = 6.107 * 10 ** ((7.5 * Ta) / (237.3 + Ta)) * (RH)
    # Determination of clear - sky emissivity from Prata (1996)
    msteg = 46.5 * (ea / (Ta + 273.15))
    esky = (1 - (1 + msteg) * np.exp(-((1.2 + 3.0 * msteg) ** 0.5))) + 0  # -0.04 old error from Jonsson et al.2006

    if altitude > 0: # # # # # # DAYTIME # # # # # #
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        #  factor for low sun elevations after Lindberg et al.(2008)
        I0, CI, Kt = clearnessindex_2013b(zen, jday, Ta, RH, radG, location)
        if (CI > 1) or (CI == np.inf):
            CI = 1

        # Estimation of radD and radI if not measured after Reindl et al.(1990)
        if 1 == 1:
            I0, CI, Kt = clearnessindex_2013b(zen, jday, Ta, RH, radG, location)
            if (CI > 1) or (CI == np.inf):
                CI = 1

            radI, radD = diffusefraction(radG, altitude, Kt, Ta, RH)

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        if True:
            zenDeg = zen*(180/np.pi)
            lv = Perez_v3(zenDeg, azimuth, radD, radI, jday)   # Relative luminance

            aniLum = 0.
            for idx in range(0, 145):
                aniLum = aniLum + diffsh[idx] * lv[0][idx][2]     # Total relative luminance from sky into each cell

            dRad = aniLum * radD   # Total diffuse radiation from sky into each cell


        # # # Surface temperature parameterisation during daytime # # # #
        # new using max sun alt.instead of  dfm
        Tg = 0.37 * altmax * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (15 / 24 - SNUP / 24)) * np.pi / 2) -3.41 # 2015 a, based on max sun altitude
        Tgwall = 0.58/0.37*(Tg+3.41) -3.41 # 2015a, based on max sun altitude

        if Tgwall < 0:  # temporary for removing low Tg during morning 20130205
            # Tg = 0
            Tgwall = 0

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
        F_sh = float(cylindric_wedge(zen, svfalfa))  # Fraction shadow on building walls based on sun alt and svf
        #F_sh[np.isnan(F_sh)] = 0.5

        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        Kdown = radI * np.sin(altitude * (np.pi / 180)) + dRad  # *sin(altitude(i) * (pi / 180))
        
        Kup = 0.15 * (radI * np.sin(altitude * (np.pi / 180.))) + dRad 

        Kesnw, KsideI, KsideD = Kside_veg_v2019a(radI, radD, radG, altitude, F_sh, Kup, lv, diffsh)
    else:  # # # # # # # NIGHTTIME # # # # # # # #
        Tgwall = 0
        # Nocturnal K fluxes set to 0
        Knight = 0.
        Kdown = 0.
        Kup = 0.
        Kesnw = 0.#Keast, Kwest, Knorth, Ksouthö
        KsideI = 0.
        KsideD = 0.
        F_sh = 0.
        Tg = 0.
        CI = 1.
        # # # # Lup # # # #
        Lup = SBC * 0.95 * ((Knight + Ta + Tg + 273.15) ** 4)
        I0 = 0

    # # # # Ldown # # # # 
    Ldown = (0.6 + 1 - 1) * esky * SBC * ((Ta + 273.15) ** 4) + \
            (1 - 0.6) * 0.9 * SBC * ((Ta + 273.15 + Tgwall) ** 4) + \
            (2 - 0.6 - 1) * (1 - 0.9) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)

    if CI < 0.95:  # non - clear conditions
        c = 1 - CI
        Ldown = Ldown * (1 - c) + \
                c * (0.6 * SBC * ((Ta + 273.15) ** 4) + (1 - 0.6) * 0.9 * SBC * ((Ta + 273.15 + Tgwall) ** 4) +
                     (1 - 0.6) * (1 - 0.9) * esky * SBC * ((Ta + 273.15) ** 4))

    # # # # Lside # # # # FIXAD
    Lsum = Lside_veg_v2020a(azimuth, altitude, Ta, Tgwall, SBC, Ldown,esky, t, F_sh, CI, Lup)

    # # # # Calculation of radiant flux density and Tmrt # # # #
    # Human body considered as a cylinder with Perez et al. (1993)
    Sstr = 0.70 * ((KsideI + KsideD) * Fcyl + (Kdown + Kup) * Fup + (4*Kesnw) * Fside) + 0.95 * ((Ldown + Lup) * Fup + (Lsum) * Fside) #4 directions

    Tmrt = float(np.sqrt(np.sqrt((Sstr / (0.95 * SBC)))) - 273.15)
    return Tmrt
