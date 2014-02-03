'''
Script to run drifters from river inputs forward for 4 months
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')
        
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc)

    overallstartdate = datetime(2007, 5, 1, 12, 1)
    overallstopdate = datetime(2007, 9, 1, 12, 1)
    # overallstartdate = datetime(2008, 5, 1, 12, 1)
    # overallstopdate = datetime(2008, 9, 1, 12, 1)

    date = overallstartdate

    # Start from the beginning and add days on for loop
    # keep running until we hit the next month
    while date < overallstopdate:

        name = date.isoformat()[0:13]

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc') and \
            not os.path.exists('tracks/' + name + 'gc.nc'):

            # Read in simulation initialization
            nstep, N, ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, do3d, doturb, \
                    grid, dostream, T0, U, V = init.init(date, loc, grid=grid)
            # pdb.set_trace()
            # Run tracpy
            # Save directly to grid coordinates
            lonp, latp, zp, t, grid, T0, U, V \
                = tracpy.run.run(loc, nstep, ndays, ff, date, tseas, ah, av, \
                                    lon0, lat0, z0, zpar, do3d, doturb, name, N=N,  \
                                    grid=grid, dostream=dostream, T0=T0, U=U, V=V, savell=False)

        # # If basic figures don't exist, make them
        # if not os.path.exists('figures/' + name + '*.png'):

            # # Read in and plot tracks
            # d = netCDF.Dataset('tracks/' + name + '.nc')
            # lonp = d.variables['lonp'][:]
            # latp = d.variables['latp'][:]
            # # tracpy.plotting.tracks(lonp, latp, name, grid=grid)
            # # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
            # d.close()
            # # # Do transport plot
            # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=date.isoformat()[0:13], 
            #     extraname=date.isoformat()[0:13], 
            #     Title='Transport on Shelf, for a week from ' + date.isoformat()[0:13], dmax=1.0)

        # Increment by 24 hours for next loop, to move through more quickly
        # nh = nh + 24
        date = date + timedelta(hours=24)

    # # Do transport plot
    # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=startdate.isoformat()[0:7] + '*', 
    #     extraname=startdate.isoformat()[0:7], Title='Transport on Shelf', dmax=1.0)


if __name__ == "__main__":
    run()    
