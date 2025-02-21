import numpy as np
import pandas as pd
import os
#import matplotlib.pyplot as plt

def getlaskar2004(option=1, timeslice=(-np.inf, np.inf)):
    """
    tka, ecc, obl, lpe = getlaskar2004(option=1, timeslice=(-np.inf, np.inf))

    Open data files for the Laskar2004 solution (Laskar et al., 2004).
    Downloaded from http://vo.imcce.fr/insola/earth/online/earth/La2004/index.html

    Parameters
    ----------
    option : integer
        option = 1, 51 Ma to 0 Ma
        option = 2, 0 Ma to 21 Ma in the future
        option = 3, 101 Ma to 0 Ma
        option = 4, 249 Ma to 0 Ma
        option = 5, 51 Ma to 21 Ma in the future (concatenate options 1 & 2)
    timeslice : array-like, containing one or two values
        If one value, a single time interval. If two values, minimum and maximum time interval (in ka before 2000 CE)
        If not given, all time slices in the dataset will be returned.

    Returns
    -------
    tka, ecc, obl, lpe

    tka: ndarray
        time in ka before year 2000 CE (negative values = future from 2000 CE)
    ecc: ndarray
        eccentricity (dimensionless: https://en.wikipedia.org/wiki/Orbital_eccentricity)
    obl: ndarray
        obliquity (radians)
    lpe: ndarray
        longitude of perihelion from moving equinox (radians, heliocentric)

    Information
    -----------
    Script originally written in Matlab by B.C. Lougheed in 2020.
    Ported to python/pandas/numpy by B.C. Lougheed in Oct. 2024.

    Reference for the imported data:
    Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B., 2004.
    "A long-term numerical solution for the insolation quantities of the Earth. "
    A&A 428, 261-285. https://doi.org/10.1051/0004-6361:20041335
    """
    file_dir = os.path.dirname(os.path.abspath(__file__))
    dirloc = file_dir+'/laskar_et_al/'

    def scinot(val): # because Laskar et al output D (from fortran) instead of E for overflow
        return float(val.replace('D', 'E'))
 
    if option == 1:
        d = pd.read_csv(dirloc+'INSOLN.LA2004.BTL.ASC', sep=r'\s+', header=None, skiprows=3, converters={1: scinot, 2: scinot, 3: scinot}).values
    elif option == 2:
        d = pd.read_csv(dirloc+'INSOLP.LA2004.BTL.ASC', sep=r'\s+', header=None, skiprows=3, converters={1: scinot, 2: scinot, 3: scinot}).values
    elif option == 3:
        d = pd.read_csv(dirloc+'INSOLN.LA2004.BTL.100.ASC', sep=r'\s+', header=None, skiprows=3, converters={1: scinot, 2: scinot, 3: scinot}).values
    elif option == 4:
        d = pd.read_csv(dirloc+'INSOLN.LA2004.BTL.250.ASC', sep=r'\s+', header=None, skiprows=3, converters={1: scinot, 2: scinot, 3: scinot}).values
    elif option == 5:
        d1 = pd.read_csv(dirloc+'INSOLN.LA2004.BTL.ASC', sep=r'\s+', header=None, skiprows=3, converters={1: scinot, 2: scinot, 3: scinot}).values
        d2 = pd.read_csv(dirloc+'INSOLP.LA2004.BTL.ASC', sep=r'\s+', header=None, skiprows=3, converters={1: scinot, 2: scinot, 3: scinot}).values
        d = np.vstack((d1[1:], d2))
    
    # prep the data
    d[:, 0] = (d[:, 0] * -1) # geo style: make past positive, future negative
    d = d[d[:, 0].argsort()][::-1]  # sort rows by the first column in descending order

    # get the requested time slice
    timeslice = np.array([timeslice])
    d = d[(d[:, 0] >= np.min(timeslice)) & (d[:, 0] <= np.max(timeslice))]

    # return np arrays
    tka = np.array(d[:, 0])
    ecc = np.array(d[:, 1])
    obl = np.array(d[:, 2])
    lpe = np.array(d[:, 3])

    return tka, ecc, obl, lpe

def getlaskar2010(option=1, timeslice=(-np.inf, np.inf)):
    """
    tka, ecc = getlaskar2010(option)

    Open Laskar2010 eccentricity solution data files. Useful for looking at eccentricity > 30 Ma.
    Downloaded from http://vo.imcce.fr/insola/earth/online/earth/La2010/index.html

    Parameters
    ----------
    option: integer
        option = 1, La2010a (solution a)
        option = 2, La2010b (solution b)
        option = 3, La2010c (solution c)
        option = 4, La2010d (solution d)
    timeslice : array-like, containing one or two values
        If one value, a single time interval. If two values, minimum and maximum time interval (in ka before 2000 CE)
        If not given, all time slices in the dataset will be returned.

    Returns
    -------
    tka, ecc

    tka: ndarray
        time in ka before year 2000 CE (negative years = future from 2000 CE)
    ecc: ndarray
        eccentricity (dimensionless: https://en.wikipedia.org/wiki/Orbital_eccentricity)

    Information
    -----------
    Script originally written in Matlab 2019a by B.C. Lougheed in 2020.
    Ported to python/pandas/numpy by B.C. Lougheed in Oct. 2024.
    Python 3.12.4, pandas 2.2.2, numpy 1.26.4.

    Reference for the imported data:
    Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011.
    "La2010: a new orbital solution for the long-term motion of the Earth."
    A&A 532, A89. https://doi.org/10.1051/0004-6361/201116836
    """

    file_dir = os.path.dirname(os.path.abspath(__file__))
    dirloc = file_dir+'/laskar_et_al/'
    
    if option == 1:
        d = pd.read_csv(dirloc+'La2010a_ecc3L.dat', sep=r'\s+', header=None).values
    elif option == 2:
        d = pd.read_csv(dirloc+'La2010b_ecc3L.dat', sep=r'\s+', header=None).values
    elif option == 3:
        d = pd.read_csv(dirloc+'La2010c_ecc3L.dat', sep=r'\s+', header=None).values
    elif option == 4:
        d = pd.read_csv(dirloc+'La2010d_ecc3L.dat', sep=r'\s+', header=None).values

    d[:, 0] = (d[:, 0] * -1) - (50 / 1000)  # ka 1950
    d = d[d[:, 0].argsort()][::-1]  # sort rows by the first column in descending order

    # prep the data
    d[:, 0] = (d[:, 0] * -1) # geo style: make past positive, future negative
    d = d[d[:, 0].argsort()][::-1]  # sort rows by the first column in descending order

    # get the requested time slice
    timeslice = np.array([timeslice]) 
    d = d[(d[:, 0] >= np.min(timeslice)) & (d[:, 0] <= np.max(timeslice))]

    tka = np.array(d[:, 0])
    ecc = np.array(d[:, 1])

    return tka, ecc

def solvekeplerE(M, ecc, floatpp=64):
    """
    E, precision = solvekeplerE(M, ecc)

    Solves Kepler equation for E to within machine precision by
    using Sinnot (1985) binary search method.

    Suitable for eccentricity values between 0 and 0.98.

    Parameters
    ----------
    M : ndarray
        Mean anomaly (radians)
    ecc : ndarray
        Eccentricity of the ellipse (ratio)
    floatpp : integer (optional)
        Floating point precision you are using (default = 64)

    Returns
    -------
    E : ndarray
        Eccentric anomaly (radians)
    precision : ndarray
        The precision on the eccentric anomaly (difference between the final two loop iterations).

    Information
    -----------
    Roger Sinnott (1985) BASIC script, as suggested by Meeus (1998).
    Solves Kepler equation for E, using binary search, to within computer precision.
    Ported to Matlab by Tiho Kostadinov (Kostadinov and Gilb, 2014). 
    Matlab version vectorised to handle array input by Bryan Lougheed (Lougheed, 2022). 
    Subsequently ported to python/numpy in October 2024 by Bryan Lougheed.
    Python 3.12.4, numpy 1.26.4.

    References
    ----------
    R.W. Sinnott (1985), "A computer assault on Kepler's equation", Sky and Telescope, vol. 70, page 159.
    J. Meeus (1998). Chapter 30 in Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia.
    Kostadinov and Gilb, (2014): doi:10.5194/gmd-7-1051-2014.
    B.C. Lougheed (2022), doi:10.5334/oq.100.
    """
    M, ecc = np.broadcast_arrays(M, ecc) # not necessary in matlab, needed here in numpy

    F = np.sign(M)
    M = np.abs(M) / (2*np.pi)
    M = (M-np.floor(M)) * 2*np.pi * F
    M[M<0] += 2*np.pi  # put in same relative orbit
    F = np.ones_like(M)

    mask = M>np.pi
    F[mask] = -1
    M[mask] = 2*np.pi - M[mask]  # inbound
    
    Eo = np.full_like(ecc, np.pi/2)
    D = np.full_like(ecc, np.pi/4)
    # Converging loop
    # Sinnot says number of iterations is 3.30 * significant figures of system. 
    # Matlab double has 16 digit precision, so 16*3.30=53. Let's use 58
    # np int64 is the same
    iters = np.ceil((floatpp/4)*3.3+5)
    for i in range(int(iters)): 
        if i == int(iters-1):
            Eoo = np.copy(Eo)
        M1 = Eo - ecc * np.sin(Eo)
        Eo += D * np.sign(M - M1)
        D /= 2

    E = Eo * F
    precision = np.abs(Eo - Eoo)
        
    return E, precision

def sollon2time(sollon, ecc, lpe, tottime=365.24, obl=None):
    """
    time, eot = sollon2time(sollon, ecc, lpe, tottime=365.24, obl=None)

    Given a particular eccentricity and longitude of perihelion, get time of tropical year
    associated with a particular geocentric solar longitude, i.e. by accounting for 
    conservation of angular momentum during orbit (Kepler 2nd Law).

    Parameters
    ----------
    sollon : array-like
        Keplerian geocentric solar longitude in radians ('v' relative to NH spring equinox)
        Either 1 value (used as constant if other inputs are vector), or a vector of values.
    ecc : array-like
        Eccentricity (e.g., from Laskar et al.)
    lpe : arriay-like
        Heliocentric longitude of perihelion (a.k.a omega-bar) in radians (e.g., from Laskar et al.)
    tottime : float
        Total time in the year, single value, any time unit you want. Default value is 365.24.
    obl : array-like, optional
        Obliquity in radians (e.g., from Laskar et al.) for calculating the equation of time (eot)

    Returns
    -------
    time, eot

    time : ndarray
        Time interval of tropical year (where 0 is boreal spring equinox).
    eot : ndarray
        Equation of time (minutes). Returns empty if obl not supplied.

    Information
    -----------
    Bryan Lougheed, June 2020, Matlab 2019a
    Updated April 2023 to include eot.
    Converted to python/numpy October 2024 by Bryan Lougheed.
    Python 3.12.4, numpy 1.26.4.
    
    See following for background, as well as comments in the script:
    J. Meeus, (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30).
    https://dr-phill-edwards.eu/Science/EOT.html (for equation of time)
    """

    # Change lpe from heliocentric to geocentric
    omega = lpe + np.pi
    omega[omega >= 2*np.pi] -= 2*np.pi  # wrap to 360

    # Get day of anchor day (dz) relative to perihelion
    vz = 2*np.pi - omega  # v of spring equinox relative to perihelion
    vz[vz > 2*np.pi] -= 2*np.pi
    Ez = 2 * np.arctan(np.tan(vz / 2) * np.sqrt((1 - ecc) / (1 + ecc)))  # Meeus (1998) page 195, solve for E
    Mz = Ez - ecc * np.sin(Ez)  # Meeus page 195, solve for M (Kepler equation). M is the circular orbit equivalent of v
    Mz[Mz < 0] = np.pi + (np.pi - Mz[Mz < 0] * -1)  # inbound to perihelion
    dz = Mz / (2*np.pi) * tottime

    # Get day of target day (dx) relative to perihelion
    vx = vz + sollon
    vx[vx > 2*np.pi] -= 2 * np.pi
    Ex = 2 * np.arctan(np.tan(vx / 2) * np.sqrt((1 - ecc) / (1 + ecc)))  # Meeus (1998) page 195, solve for E
    Mx = Ex - ecc * np.sin(Ex)  # Solve for M (Kepler equation)
    Mx[Mx<0] = np.pi + (np.pi - Mx[Mx < 0] * -1)  # inbound to perihelion, (probably not necessary)
    dx = Mx / (2*np.pi) * tottime

    # Get day of target day (dx) relative to day of anchor day (dz)
    dx[dx<dz] += tottime  # for dz in next orbital period relative to perihelion, keep in same orbital period relative to NH spring equinox
    time = dx - dz

    # Eliminate rounding errors at zero
    sollon, time = np.broadcast_arrays(sollon, time)
    time[sollon == 0] = 0
    time[sollon == 2*np.pi] = 0

    # Calculate equation of time if obl is supplied
    # https://dr-phill-edwards.eu/Astrophysics/EOT.html (explains it very nicely) 
    if obl is not None:
        # eccentricity component
        dtecc = np.rad2deg(Mx-vx) * 4  # four minutes per degree (24 hrs * 60 mins / 360 degrees )
        # obliquity component
        alpha = np.arctan2(np.sin(sollon) * np.cos(obl), np.cos(sollon))
        alpha[alpha<0] += 2*np.pi
        dtobl = np.rad2deg(sollon-alpha) * 4 # same here
        # total EOT, time in minutes
        eot = dtecc + dtobl
    else:
        eot = np.array([])

    return time, eot

def time2sollon(time, ecc, lpe, tottime=365.24, obl=None, floatpp=64):
    """
    sollon, eot = time2sollon(time, ecc, lpe, tottime=365.24, obl=None, floatpp=64)

    Given a particular eccentricity and longitude of perihelion, get geocentric solar longitude 
    associated with a particular time of the tropical year i.e. by accounting for 
    conservation of angular momentum during orbit (Kepler 2nd Law).

    Parameters
    ----------
    time : ndarray
        Time interval of tropical year (where interval 0 is boreal spring equinox). Either on value, or vector of values.
    ecc : ndarray
        Eccentricity (e.g., from Laskar et al.)
    lpe : ndarray
        Heliocentric longitude of perihelion (a.k.a omega-bar) in radians (e.g., from Laskar et al.)
    tottime : float
        Total time in the year corresponding to 'time', single value, any time unit you want. Default value is 365.24.
    obl : ndarray, optional
        Obliquity in radians (e.g., from Laskar et al.) for calculating the equation of time (eot)
    floatpp : integer, optional
        Floating point precision you are using (default = 64 bit)

    Input can be vectorised in various ways, but double check output.

    Returns
    -------
    sollon : np.array
        Keplerian geocentric solar longitude in radians ('lambda', i.e. 'v' relative to boreal spring equinox) 
    eot : np.array
        Equation of time (minutes). Returns empty if obl not supplied.

    Info
    ----
    B.C. Lougheed, June 2020, Matlab 2019a
    Updated April 2023 to include eot.
    Ported to python/numpy in October 2024 by B.C. Lougheed.
    Python 3.12.4, numpy 1.26.4.

    See following for background, as well as comments in the script:
    R.W. Sinnott (1985), "A computer assault on Kepler's equation." Sky and Telescope, vol. 70, page 159.
    Meeus, J., (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30).
    Kostadinov and Gilb, (2014): doi: 10.5194/gmd-7-1051-2014
    B.C. Lougheed (2022), doi:10.5334/oq.100.
    https://dr-phill-edwards.eu/Science/EOT.html (for equation of time)
    """
    # convert input to numpy for speed and compatibility
    tottime = np.array([tottime])
    time = time.reshape(-1,1)

    # change lpe from heliocentric to geocentric
    omega = np.array(lpe + np.pi) # np.array needed for when doing only one timeslice
    omega[omega >= 2*np.pi] = omega[omega >= 2*np.pi] - 2*np.pi

    # NH spring equinox relative to perihelion
    veq = 2*np.pi - omega
    Eeq = 2 * np.arctan(np.tan(veq/2) * np.sqrt((1-ecc) / (1+ecc)))
    Meq = np.array(Eeq - ecc * np.sin(Eeq)) # as previous comment
    Meq[Meq<0] = np.pi + (np.pi - Meq[Meq<0] * -1)
    deq = Meq / (2*np.pi) * tottime

    # v of target (x) v relative to perihelion
    deq, time = np.broadcast_arrays(deq,time)
    dx = deq + time
    Mx = (dx / tottime) * 2*np.pi
    Ex, _ = solvekeplerE(Mx, ecc, floatpp=floatpp)  # Get Ex by solving Kepler equation
    vx = 2 * np.arctan(np.tan(Ex/2) * np.sqrt((1+ecc) / (1-ecc)))
    vx[vx<0] = np.pi + (np.pi - vx[vx<0] * -1)
    
    # target day's v relative to NH spring equinox v
    vx[vx<veq] += 2*np.pi
    sollon = vx - veq

    # eliminate rounding errors at 0
    sollon, time = np.broadcast_arrays(sollon, time)
    sollon[time == 0] = 0
    sollon[time == tottime] = 0
    sollon = np.array(sollon)

    if obl is not None:
        # eccentricity component
        dtecc = np.rad2deg(Mx - vx) * 4
        # obliquity component
        alpha = np.arctan2(np.sin(sollon) * np.cos(obl), np.cos(sollon))
        alpha[alpha<0] += 2*np.pi
        dtobl = np.rad2deg(sollon - alpha) * 4
        # total EOT, time in minutes
        eot = dtecc + dtobl
    else:
        eot = np.array([])

    return sollon, eot

def geographiclat(gclat, angles='rad'):
    """
    gplat = geographiclat(gclat, angles='rad')

    Convert geocentric latitude into geographic latitude
    assuming the WGS84 spheroid.

    Parameters
    ----------
    gclat : array-like
        Geocentric latitude.
    angles : string (optional)
        'rad' (default) or 'deg'. 
        Specify if gclat is in degrees or radians.

    Returns
    -------
    gplat : ndarray
        Geographic latitude in radians.

    Bryan Lougheed, February 2025.
    """
    gclat = np.array(gclat)

    if angles == 'rad':
        pass
    elif angles == 'deg':
        gclat = np.deg2rad(gclat)
    else:
        raise Exception("'angles' parameter should be set to either 'deg' or 'rad'")

    # calculate geographic latitude from geocentric latitude
    f = 1 / 298.257223563  # wgs84 flattening value
    re = 6378137.0  # wgs84 equatorial radius (metres)
    rp = re * (1 - f)  # calculate polar radius
    gplat = np.arctan((re / rp)**2 * np.tan(gclat))

    return gplat

def dailymeanwm2(lat, sollon, ecc, obl, lpe, con=1361, earthshape='sphere'):
    """
    irr, dayhrs, rx, tsi = dailymeanwm2(lat, sollon, ecc, obl, lpe, con=1361, earthshape='sphere')

    Calculate 24-hr mean irradiance (W/m²) at top of atmosphere and also length of daytime (in hours), 
    total solar irradiance (TSI; in W/m²) and distance from sun (in AU).

    Parameters
    ----------
    lat : array-like
        Geocentric latitude (plus for N, minus for S) on Earth, in radians.
    sollon : array-like
        Geocentric solar longitude (lambda), in radians.
    ecc : array-like
        Eccentricity. Numerical value(s). 1D array.
    obl : array-like
        Obliquity. Numerical value(s), radians. 1D array.
    lpe : array-like
        Longitude of perihelion from moving equinox (omega-bar). Numerical value(s), radians. 1D array.
    con : float or array-like, optional
        Solar constant in W/m². Single numerical value or 1D array. Default is 1361 W/m².
    earthshape : str, optional
        Shape of Earth, enter string 'sphere' or 'wgs84' (default is 'sphere').

    Returns
    -------
    irr, dayhrs, tsi, rx

    irr : ndarray
        Calculated mean daily (24 hr) irradiance (W/m²) at top of atmosphere. Array same size as ecc, obl and lpe.
    dayhrs : ndarray
        Hours of daylight. Array same size as ecc, obl and lpe.
    rx : ndarray
        Distance from Sun (AU). Insensitive to latitude, obliquity or earthshape.
    tsi : ndarray
        Calculated mean daily irradiance at top of atmosphere assuming 90 degree angle of incidence, W/m².
        Insensitive to latitude, obliquity or earthshape. Array same size as ecc, obl and lpe.

    Info
    ----
    B.C. Lougheed, May 2020, Matlab 2019a
    Updated to include Earth's oblateness Sep. 2020
    Ported to python/numpy Oct. 2024 by B.C. Lougheed
    Python 3.12.4, numpy 1.26.4.

    irr in Wm² based on equations in Berger (1978).
    Berger, AL. 1978. "Long-Term Variations of Daily Insolation and Quaternary Climatic Changes."
    J. Atmos. Sci., 35: 2362-2367.
 
    Part of script (specifically polar night and day) uses Ian Eisenman's
    Matlabification of Berger (1978) equations 8 and 9 (see comments in script)
    taken from here: http://eisenman.ucsd.edu/code/daily_insolation.m
    
    I added ability to take oblateness of Earth into account, validated against Van Hemelrijck (1983) solution. (see comments in script)
    
    I added daylight hours output following sunrise equation: 
    https://en.wikipedia.org/wiki/Sunrise_equation
    
    I added tsi by calculating distance from sun following Meeus (1998).
    Meeus, J., (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30)    
    """

    if earthshape == 'sphere':
        # geographic latitude = geocentric latitude
        pass
    elif earthshape == 'wgs84':
        lat = geographiclat(lat)
    else:
        raise ValueError('earthshape '+earthshape+' unrecognised')

    # Declination angle of the sun
    # https://en.wikipedia.org/wiki/Position_of_the_Sun
    dsun = np.arcsin(np.sin(obl) * np.sin(sollon))

    # Hour angle at sunrise/sunset
    # https://en.wikipedia.org/wiki/Sunrise_equation
    hangle = np.arccos(-np.tan(lat) * np.tan(dsun))

    # polar night / day. Berger (1978), eq. 8 and 9
    # following two lines come from Ian Eisenman Matlab script of equations 8 and 9
    # numpy prefers np.logical_and here
    hangle[np.logical_and(np.abs(lat) >= np.pi / 2 - np.abs(dsun), lat * dsun > 0)] = np.pi  # polar day
    hangle[np.logical_and(np.abs(lat) >= np.pi / 2 - np.abs(dsun), lat * dsun <= 0)] = 0    # polar night

    # Hours of daylight (https://en.wikipedia.org/wiki/Sunrise_equation)
    dayhrs = np.abs(hangle - hangle*-1) / (2*np.pi / 24)

    # Change lpe from heliocentric to geocentric (omega-bar to omega)
    omega = np.array(lpe + np.pi)  # add 180 degrees
    omega[omega >= 2*np.pi] -= 2*np.pi  # subtract 360 degrees from stuff >= 360 degrees

    # # Van Hemelrijck (1983) extended method for oblate Earth
    # if earthshape == 'wgs84'
    #     vangle = np.arctan((1 - f)**-2 * np.tan(lat)) - lat  # Van Hemelrijck (1983) eq. 9, f is wgs84 flattening
    #     # Van Hemelrijck (1983) eq. 11 second term
    #     hemelterm = (np.cos(vangle) * (hangle * np.sin(lat) * np.sin(dsun) + np.sin(hangle) * np.cos(lat) * np.cos(dsun)) + np.sin(vangle) * (-np.tan(lat) * (hangle * np.sin(lat) * np.sin(dsun) + np.sin(hangle) * np.cos(lat) * np.cos(dsun)) + hangle * np.sin(dsun) / np.cos(lat)))
    #     # Irradiation: Berger (1978) eq. 10, but replace final term with Van Hemelrijck (1983) eq. 11 second term
    #     irr = con / np.pi * (1 + ecc * np.cos(sollon - omega))**2 / (1 - ecc**2)**2 * hemelterm
    # # produces exact same output as simply inputting the corrected latitude (calculated at top of this function) into Berger equation below, which seems easier.. 

    # 24 hr mean irradiance: Berger (1978) eq (10)
    irr = con / np.pi * (1 + ecc * np.cos(sollon - omega))**2 / (1 - ecc**2)**2 * ( hangle * np.sin(lat) * np.sin(dsun) + np.cos(lat) * np.cos(dsun) * np.sin(hangle))

    # Calculate rx and tsi
    veq = 2*np.pi - omega  # v (true anomaly) of spring equinox relative to perihelion
    vx = veq + sollon  # v (true anomaly) of inputted sollon relative to perihelion
    vx[vx > 2*np.pi] -= 2*np.pi  # put back in 0-360 range
    rx = (1 - ecc**2) / (1 + ecc * np.cos(vx))  # Eq. 30.3 in Meeus (1998)
    tsi = con * (1 / rx)**2  

    return irr, dayhrs, rx, tsi

def intradaywm2(lat, ecc, obl, lpe, con=1361.0, daysinyear=365.240, dayres=0.001):
    """
    irr, days, dayher = intradaywm2(lat, ecc, obl, lpe, con=1361, daysinyear=365.240, dayres=0.001)

    Calculate a tropical year's worth of intraday irradiance (W/m²) for a particular
    latitude and orbital configuration. Calculates for a longitude where
    northern spring equinox occurs at day 0.0 (i.e., at exactly local midnight 
    on the first day of the tropical year). Script takes equation of time into account.

    Parameters
    ----------
    lat : ndarray
        geocentric latitude in radians (positive for north, negative for south)
    ecc : ndarray
        eccentricity of the ellipse, ratio (from, e.g., Laskar et al.)
    obl : ndarray
        obliquity in radians (from, e.g., Laskar et al.)
    lpe : ndarray
        heliocentric longitude of perihelion (from e.g., Laskar et al.)
        omega-bar (i.e., relative to NH autumn equinox) in radians.
    con : float
        solar constant in W/m², default is 1361 W/m²
    daysinyear : float
        number of mean solar days in the year, default is 365.240
    dayres : float
        mean solar day resolution for analysis, default is 0.001

    Returns
    -------
    irr, days, dayher

    irr : ndarray
        array of W/m² for every day interval calculated
    days : ndarray
        all the solar day intervals. 0 corresponds to the boreal spring equinox.
    dayhr : ndarray
        time of the the mean solar day (in hours), same dimensions as irr and days

    Bryan Lougheed, April 2023, Matlab 2019a
    Ported to python/numpy by Bryan Lougheed, Oct. 2024
    Python 3.12.4, numpy 1.26.4.
    """

    # convert some input to numpy for speed (might not be necessary)
    con = np.array([con])

    # day length in hours (placeholder for future development)
    # has not been properly implemented yet for day lengths other than 24
    dlen = 24 
    
    # Create time series of mean solar day fractions
    days = np.arange(0, daysinyear, dayres)
    dayhr = (days - np.floor(days)) * dlen

    # Calculate Earth's solar longitude (i.e., lambda) and equation of time for each day fraction
    sunlon, eot = time2sollon(days, ecc, lpe, daysinyear, obl)

    # Get local apparent solar hour (correct for eot)
    dayhr = dayhr.reshape(-1,1)
    eot, dayhr = np.broadcast_arrays(eot,dayhr)
    solhr = (eot/60) + dayhr # /60, mins -> hrs
    solhr[solhr<0] += dlen
    solhr[solhr>dlen] -= dlen

    # Declination of the sun
    dsun = np.arcsin(np.sin(obl) * np.sin(sunlon))

    # Local hour angle (-180 to +180 deg, midday = 0 deg)
    hangles = (2*np.pi / dlen) * (solhr - dlen/2)

    # Solar elevation
    elev = np.arcsin(np.sin(dsun) * np.sin(lat) + np.cos(dsun) * np.cos(lat) * np.cos(hangles))

    # Calculate distance from Sun in AU
    omega = lpe + np.pi
    veq = 2 * np.pi - omega  # v (true anomaly) of NH spring equinox relative to perihelion
    vx = veq + sunlon  # v (true anomaly) of inputted sunlon relative to perihelion
    vx[vx > 2*np.pi] -= 2 * np.pi  # put back in 0-360 range
    rx = (1 - ecc**2) / (1 + ecc * np.cos(vx))  # Eq. 30.3 in Meeus (1998)

    # Calculate tsi as function of con relative to 1 AU
    tsi = con * (1/rx)**2

    # Calculate W/m2, vertical component of tsi
    irr = tsi * np.sin(elev)
    irr[irr < 0] = 0  # sun under horizon, night time

    return irr, days, dayhr

def thresholdjm2(thresh, lat, ecc, obl, lpe, con=1361, timeres=0.01, tottime=365.24, earthshape='sphere'):
    """
    intirr, ndays = thresholdjm2(thresh, lat, ecc, obl, lpe, con=1361, timeres=0.01, tottime=365.24, earthshape='sphere')

    Calculate integrated irradiation (J/m²) at top of atmosphere for all day intervals 
    exceeding a certain threshold in mean daily irradiance (W/m²).
    Can be used to emulate analysis by Huybers (2006; 10.1126/science.1125249)
    
    Parameters
    ----------
    thresh : float or array-like
        Threshold value (W/m2). Single value, or vector of values.
    lat : float
        Geocentric latitude (in deg. N, negative for S) on Earth. Single value.
    con : float or array-like, optional
        Solar constant. Single numerical value or 1D array, W/m2. Default is 1361.
    dayres : float, optional
        Day resolution for the integration. Default is 0.01.
    ecc : array-like
        Eccentricity. Numerical value(s). 1D array.
    obl : array-like
        Obliquity. Numerical value(s), radians. 1D array.
    lpe : array-like
        Longitude of perihelion from moving equinox. (omega-bar) Numerical value(s), radians. 1D array.
    earthshape : str (optional)
        Shape of Earth, 'sphere' (default) or 'wgs84'.

    Returns
    -------
    intirr : ndarray
        Integrated irradiation at top of atmosphere for days exceeding thresh. J/m2. Array same dimensions as ecc, obl, and lpe.
    ndays : ndarray
        Time (in days) exceeding thresh. Same dimensions as intirr.
    """

    timerange = np.arange(0, tottime, timeres)
    intirr = np.full_like(ecc, np.nan)
    ndays = np.full_like(ecc, np.nan)

    for i in range(len(ecc)):
        sollons, _ = time2sollon(timerange, ecc[i], lpe[i], tottime)
        irrs, _, _, _ = dailymeanwm2(lat, sollons, con, ecc[i], obl[i], lpe[i], earthshape)
        ndays[i] = np.sum(irrs >= thresh) * timeres
        intirr[i] = np.mean(irrs[irrs >= thresh]) * (ndays[i] * 24 * 60 * 60)  # W/m2 to J/m2

    return intirr, ndays

def areaquad(lat1, lat2, lon1, lon2, shape='sphere', angles='rad'):
    """
    aq = areaquad(lat1, lat2, lon1, lon2, shape='sphere', angles='rad')
    
    Calculate the surface area of a lat/lon bounding box on Earth.

    Parameters
    ----------
    lat1 : float
        A bounding geocentric latitude.
    lat2 : float
        The other bounding geocentric latitude.
    lon1 : float
        A bounding geocentric longitude.
    lon2 : float
        The other bounding gecentric longitude.
    shape : string (optional)
        'sphere' (default) or 'wgs84'
        'sphere' will assume a sphere with a radius of 6371008.7714 metres.
        'wgs84' will assume an oblate Earth with a semi-major axis 
        of 6378137.0 metres and a first eccentricity of 0.0818191908426215.
    angles : string (optional)
        'rad' (default) or 'deg'. 
        Specify if lat1, lat2, lon1 and lon2 are in degrees or radians.

    Returns
    -------
    aq : float
        The area of the bounding box, given in square metres.

    Bryan Lougheed, February 2025

    This a python/numpy simplified port of the Octave function areaquad.m from the 
    Octave "mapping" package (v.1.4.2) (https://gnu-octave.github.io/packages/mapping/),
    which included the following license:

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Full license text available at: http://www.gnu.org/licenses/
    """
    if angles == 'rad':
        pass
    elif angles == 'deg':
        lat1 = np.deg2rad(lat1)
        lat2 = np.deg2rad(lat2)
        lon1 = np.deg2rad(lon1)
        lon2 = np.deg2rad(lon2)
    else:
        raise Exception("'angles' parameter should be set to either 'deg' or 'rad'")

    if shape == 'sphere':
        a = 6371008.7714
        e = 0
    elif shape == 'wgs84':
        a = 6378137.0
        e = 0.0818191908426215
    else:
        raise Exception("'shape' parameter should be set to either 'sphere' or 'wgs84'")

    s1 = np.sin(lat1)
    s2 = np.sin(lat2)
    lonwidth = lon1 - lon2
    
    if e < np.finfo(float).eps:
        aq = abs((lonwidth * a**2) * (s2 - s1))
    else:
        e2 = e**2
        f = 1 / (2 * e)
        e2m1 = 1 - e2

        s21 = s1**2
        s22 = s2**2
        se1 = 1 - e2 * s21
        se2 = 1 - e2 * s22

        c = (lonwidth * a**2 * e2m1) / 2
        t1 = 1 + e * s1
        t2 = 1 + e * s2
        b1 = 1 - e * s1
        b2 = 1 - e * s2

        g = f * (np.log(t2 / b2) - np.log(t1 / b1))

        aq = np.abs( c * ((s2 / se2) - (s1 / se1) + g) )
  
    return aq

##### under construction

# def plotorbit(ecc, lpe, savename):
#     """
#     Creates a scale plot of the Earth's orbit and saves to hard drive.

#     Parameters
#     ----------
#     ecc : float
#         Eccentricity of the ellipse (ratio).
#     lpe : float
#         Longitude of perihelion as given by, e.g., Laskar: omega-bar 
#         (i.e. relative to NH autumn equinox) in radians.
#     savename : str
#         The name of the output .png file.
    
#     Output
#     ------
#     A scale plot of the orbit in current plotting window, or new plotting
#     window if there is none open.
#     """

#     # true anomaly v of seasons
#     omega = lpe + np.pi  # perihelion (geocentric)
#     nse = 2 * np.pi - omega  # nh spring equinox
#     nss = nse + np.pi/2  # nh summer solstice
#     nae = nss + np.pi/2  # nh autumn equinox
#     nws = nae + np.pi/2  # nh winter solstice

#     # convert true anomaly v to ellipse angle E (Meeus 1998 Equation 30.1 solve for E)
#     nse = 2 * np.arctan(  np.tan(nse / 2) * np.sqrt((1 - ecc)/(1 + ecc))  )
#     nss = 2 * np.arctan(  np.tan(nss / 2) * np.sqrt((1 - ecc)/(1 + ecc))  )
#     nae = 2 * np.arctan(  np.tan(nae / 2) * np.sqrt((1 - ecc)/(1 + ecc))  )
#     nws = 2 * np.arctan(  np.tan(nws / 2) * np.sqrt((1 - ecc)/(1 + ecc))  )

#     # semi-minor and semi-major axis of ellipse
#     a = 149.5978707  # au in 10^6 km
#     b = a*(1-ecc**2)**0.5
#     # perihelion and aphelion distances
#     rper = a * (1 - ecc)
#     raph = a * (1 + ecc)
#     fd = a - rper  # focal distance offset from centre of ellipse

#     # ellipse circumference points x and y coords
#     E = np.linspace(0, 2*np.pi, 1000)  # central angle E of ellipse, 0 to 360
#     xell = a*np.cos(E)-fd  # -fd to place sun at centre of image
#     yell = b*np.sin(E)

#     # plot
#     plt.figure()
#     plt.clf()

#     # plot the ellipse
#     plt.plot(xell, yell, '-', color=[0.6, 0.6, 0.6], linewidth=1)

#     # plot the sun
#     plt.plot(0, 0, 'yo', markersize=15, markerfacecolor=[255/255, 221/255, 66/255], markeredgecolor=[255/255, 143/255, 66/255])

#     # plot the astronomical seasons
#     ofs = 1.13  # text offset
    
#     def plot_season(v, label, color):
#         plt.plot(a * np.cos(v) - fd, b * np.sin(v), 'ko', markersize=7, markerfacecolor=color, markeredgecolor=[0, 100 / 255, 0])
#         plt.text(a * np.cos(v) * ofs - fd, b * np.sin(v) * ofs, label, horizontalalignment='center', rotation=np.rad2deg(v))

#     plot_season(nse, 'NH spring equinox', [0, 82 / 255, 162 / 255])
#     plot_season(nss, 'NH summer solstice', [0, 82 / 255, 162 / 255])
#     plot_season(nae, 'NH autumn equinox', [0, 82 / 255, 162 / 255])
#     plot_season(nws, 'NH winter solstice', [0, 82 / 255, 162 / 255])

#     # plot the per and aph distance
#     plt.arrow(xell[0], 0, -rper, 0, head_width=2, head_length=3, fc='r', ec='r')
#     plt.arrow(-fd, 0, raph, 0, head_width=2, head_length=3, fc='r', ec='r')
#     plt.text(np.mean([xell[0], 0]), 13, f'Perihelion\n{rper:.1f} million km', horizontalalignment='center')
#     plt.text(np.mean([a * np.cos(np.pi) - fd, 0]), 13, f'Aphelion\n{raph:.1f} million km', horizontalalignment='center')

#     plt.xlim([-170, 170])
#     plt.ylim([-170, 170])

#     plt.axis('off')
#     plt.gcf().set_size_inches(7.5, 7.5)
#     plt.savefig(savename, dpi=150, bbox_inches='tight')
#     plt.close()



