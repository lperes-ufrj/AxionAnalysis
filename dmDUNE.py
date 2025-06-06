from astropy import coordinates as coords
from astropy import time
import astropy.units as u
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors 
import dunestyle.matplotlib as dunestyle
import ROOT


####
#Here are DUNE FD values:
#
#FD1 (NE): 1417.62 m, 3986.75 mwe; 103.7546341 W, 44.3459914 N
#
#FD2 (SE): 1392.17 m, 3866.89 mwe; 103.7548042 W, 44.3450261 N
#
#FD3 (NW): 1405.13 m, 3942.87 mwe; 103.7556732 W, 44.3460857 N
#
#FD4 (SW): 1404.21 m, 3885.87 mwe; 103.7558433 W, 44.3451204 N
#
#The highest point of the detector caverns is 126.62 m (415.42 ft) above sea level.
#
####


FluxRotValues = np.array([[+0.9877, -0.1564, +0.0000],  # new x axis in old coordinates: be nice and fix things (ends up the same)
                        [+0.0000, +0.0000, +1.0000],  # new y axis in old coordinates: vertical
                        [-0.1564, -0.9877, +0.0000]])   # new z axis in old coordinates: away from Batavia, IL

def Calc_GC_AngDist():
    # Define observer's location
    observer_location = coords.EarthLocation(lat=44.34 * u.deg, lon=-103.75 * u.deg, height= 100 * u.m)

    # Define the Galactic coordinates
    galactic_coords = coords.SkyCoord(l=0 * u.deg, b=0* u.deg, frame='galactic')

    # Define the time range
    current_time = time.Time.now() + 10 * u.year # 10 years from now
    end_time = current_time + 1 * u.year  # 20 years from now
    #test_oneday =current_time + 0.5 * u.day
    times = time.Time(np.linspace(current_time.jd, end_time.jd, 2000), format='jd')  # 200,000 time points


    # Prepare to store results
    data_honda = []
    data_det = []
    data_check = []

    # Loop over all the times to calculate directions
    for t in times:
        # Convert Galactic coordinates to ICRS
        icrs_coords = galactic_coords.transform_to(coords.ICRS)

        # Convert ICRS coordinates to AltAz (topocentric horizontal) at each time
        altaz_coords = icrs_coords.transform_to(coords.AltAz(obstime=t, location=observer_location))

        # Convert AltAz coordinates to Cartesian (custom direction system)
        alt = altaz_coords.alt.to(u.rad).value  # Altitude in radians
        az = altaz_coords.az.to(u.rad).value    # Azimuth in radians


        #z is orthogonal to earth surface pointing to the zenith, x points toward south and y toward east 
        # Calculate Cartesian coordinates where:
        # x-axis points towards geographic south
        # y-axis points towards geographic east
        # z-axis points towards zenith
        x = np.sin(az) * np.cos(alt)
        y = np.cos(az) * np.cos(alt)
        z = np.sin(alt)

        arr = np.array([x,y,z])

        arr_det = FluxRotValues.dot(arr)

        # Append the results
        data_check.append([t.iso,alt, az])

        data_honda.append([t.iso, x, y, z, math.atan2(y, x), math.acos(z)])

        data_det.append([t.iso, arr_det[0], arr_det[1], arr_det[2], math.atan2(arr_det[1], arr_det[0]), math.acos(arr_det[2])])

    # Save the data to a CSV file
    df_Honda = pd.DataFrame(data_honda, columns=['Time', 'X (South)', 'Y (East)', 'Z (Zenith)', 'Phi (Honda)', 'Theta (Honda)'])
    df_Honda.to_csv('HondaCoordSys_directions_over_10_years.csv', index=False)

    # Save the data to a CSV file
    df_Detector = pd.DataFrame(data_det, columns=['Time', 'X (RightHanded)', 'Y (Zenith)', 'Z (Beam Direction)', 'Phi (Detector)', 'Theta (Detector)'])
    df_Detector.to_csv('DetectorCoordSys_directions_over_10_years.csv', index=False)

    # Save the data to a CSV file
    df_Check = pd.DataFrame(data_check, columns=['Time', 'Alt (Altitude)', 'Az (Azimuth)'])
    df_Check.to_csv('altaz_coords_directions_over_10_years.csv', index=False)
    
    