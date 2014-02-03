'''
Functions to initialize various numerical experiments.

Make a new init_* for your application.

loc     Path to directory of grid and output files
nsteps  Number of steps to do between model outputs (iter in tracmass)
ndays   number of days to track the particles from start date
ff      ff=1 to go forward in time and ff=-1 for backward in time
date    Start date in datetime object
tseas   Time between outputs in seconds
ah      Horizontal diffusion in m^2/s. 
        See project values of 350, 100, 0, 2000. For -turb,-diffusion
av      Vertical diffusion in m^2/s.
do3d    for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
doturb  turbulence/diffusion flag. 
        doturb=0 means no turb/diffusion,
        doturb=1 means adding parameterized turbulence
        doturb=2 means adding diffusion on a circle
        doturb=3 means adding diffusion on an ellipse (anisodiffusion)
lon0    Drifter starting locations in x/zonal direction.
lat0    Drifter starting locations in y/meridional direction.
z0/zpar Then z0 should be an array of initial drifter depths. 
        The array should be the same size as lon0 and be negative
        for under water. Currently drifter depths need to be above 
        the seabed for every x,y particle location for the script to run.
        To do 3D but start at surface, use z0=zeros(ia.shape) and have
         either zpar='fromMSL'
        choose fromMSL to have z0 starting depths be for that depth below the base 
        time-independent sea level (or mean sea level).
        choose 'fromZeta' to have z0 starting depths be for that depth below the
        time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
        Then: 
        set z0 to 's' for 2D along a terrain-following slice
         and zpar to be the index of s level you want to use (0 to km-1)
        set z0 to 'rho' for 2D along a density surface
         and zpar to be the density value you want to use
         Can do the same thing with salinity ('salt') or temperature ('temp')
         The model output doesn't currently have density though.
        set z0 to 'z' for 2D along a depth slice
         and zpar to be the constant (negative) depth value you want to use
        To simulate drifters at the surface, set z0 to 's' 
         and zpar = grid['km']-1 to put them in the upper s level
         z0='s' is currently not working correctly!!!
         In the meantime, do surface using the 3d set up option but with 2d flag set
xp      x-locations in x,y coordinates for drifters
yp      y-locations in x,y coordinates for drifters
zp      z-locations (depths from mean sea level) for drifters
t       time for drifter tracks
name    Name of simulation to be used for netcdf file containing final tracks

'''

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import glob
from datetime import datetime, timedelta
from matplotlib.mlab import *
import tracpy
import time

units = 'seconds since 1970-01-01'

def init(date, loc, grid=None):
    '''
    Initialization for seeding drifters at all shelf model grid points to be run
    forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        loc     Location of model output
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    tic = time.time()

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 25 

    # Number of steps to divide model output for outputting drifter location
    N = 5

    # Number of days
    ndays = 30*4

    # This is a forward-moving simulation
    ff = 1 

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 5.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    ## Initial lon/lat locations and info for drifters ##

    # Read locs in from TXLA_river_4dyes_2011.nc
    r = netCDF.Dataset('TXLA_river_4dyes_2011.nc')

    # rho grid index for river starting locations, Miss and Atch, xi and eta directions
    xi = r.variables['river_Xposition'][:51] 
    eta = r.variables['river_Eposition'][:51] 

    # river direction: 0 is east-west and 1 is north-south
    rdir = r.variables['river_direction'][:51]

    # river discharge rate and time
    Q = r.variables['river_transport'][:,:51]
    rt = r.variables['river_time'][:] # daily
    runits = 'days since 1970-01-01'

    # Starting positions will be at the center of the grid cell wall where the river is input
    # which is determined by the starting grid cell (xi, eta) then the direction (rdir), then 
    # whether it is positive or negative (positive is east/north)
    # position should be rho - 0.5 to get onto the cell wall
    xstart0 = xi.copy(); ystart0 = eta.copy() # at rho points
    # if east-west and negative, subtract 1 in the east-west direction to put at left cell edge on u grid
    # subtract 0.5 to put north-south in the middle of the cell on vgrid
    ind = (rdir==0)*(np.sign(Q[0,:])==-1) # Q stays the same sign for each locations
    xstart0[ind] = xstart0[ind] - 1
    ystart0[ind] = ystart0[ind] - 0.5
    # if east-west and positive. rho grid index is the same as the right side u grid index
    ind = (rdir==0)*(np.sign(Q[0,:])==1)
    ystart0[ind] = ystart0[ind] - 0.5
    # north-south and negative
    ind = (rdir==1)*(np.sign(Q[0,:])==-1)
    xstart0[ind] = xstart0[ind] - 0.5
    ystart0[ind] = ystart0[ind] - 1
    # north-south and positive
    ind = (rdir==1)*(np.sign(Q[0,:])==1)
    xstart0[ind] = xstart0[ind] - 0.5
    # Change to lat/lon
    lon0temp, lat0temp, _ = tracpy.tools.interpolate2d(xstart0, ystart0, grid, 'm_ij2ll')

    ## Determine how many drifters to start based on the size of the transport at the start time ##

    # Find index in river time
    rind = find(netCDF.num2date(rt, runits)==date) # startdate should exactly equal a river time
    # find ndrifters based on the transports at that time for all rivers
    ndrifters = abs(np.round(Q[rind,:]/5)[0])

    # make ndrifters number for each starting location
    lon0 = []; lat0 = []; T0 = []
    for iloc in xrange(ndrifters.size):
        lon0.extend(np.ones(ndrifters[iloc])*lon0temp[iloc])
        lat0.extend(np.ones(ndrifters[iloc])*lat0temp[iloc])
        T0.extend(np.ones(ndrifters[iloc])*(Q[rind,iloc]/ndrifters[iloc]))

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 2

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 1

    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')



    return nsteps, N, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, grid, dostream, np.asarray(T0), np.asarray(U), np.asarray(V)
