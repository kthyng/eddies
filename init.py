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
    ndays = 30*3

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
    # go up to 50 only because 51 is on land
    idrift = 50
    xi = r.variables['river_Xposition'][:idrift].astype(int)
    eta = r.variables['river_Eposition'][:idrift].astype(int)

    # river direction: 0 is east-west and 1 is north-south
    rdir = r.variables['river_direction'][:idrift]

    # river discharge rate and time
    Q = r.variables['river_transport'][:,:idrift]
    Vshape = r.variables['river_Vshape'][-1,:idrift] # only take surface transport (-1)
    Q = Q*Vshape # river transport from surface layer only
    rt = r.variables['river_time'][:] # daily
    runits = 'days since 1970-01-01'
    rdates = netCDF.num2date(rt, runits)
    r.close()

    # Starting positions are at grid cell walls
    # which is determined by the starting grid cell (xi, eta) then the direction (rdir), then 
    # whether it is positive or negative (positive is east/north)
    # position should be rho - 0.5 to get onto the cell wall
    xstart0 = np.zeros(xi.shape)
    ystart0 = np.zeros(xi.shape)

    iu = find(rdir==0) # then the xi/eta points are u grid-referencing
    iv = find(rdir==1) # then the xi/eta points are v grid-referencing

    # The east-west forcings are referenced to u grid, and in the x-direction should be shifted
    # by one to get indexing right
    # pdb.set_trace()
    xstart0[iu] = grid['xu'][xi[iu]-1,eta[iu]]
    ystart0[iu] = grid['yu'][xi[iu]-1,eta[iu]]

    # The north-south forcings are referenced to v grid and in the y direction should be shifted
    # by one
    xstart0[iv] = grid['xv'][xi[iv],eta[iv]-1]
    ystart0[iv] = grid['yv'][xi[iv],eta[iv]-1]

    # Do some manual shifting slightly more into cell to avoid getting masked out
    # west side of miss delta
    xstart0[:11] = xstart0[:11]-10 # shift due west 10 meters
    # south side of delta
    ystart0[11:24] = ystart0[11:24]-10 # shift due south 10 meters
    # east side of delta
    xstart0[24:43] = xstart0[24:43]+10 # shift due east 10 meters
    # atchafalaya
    ystart0[43:] = ystart0[43:]-10 # shift due south 10 meters

    # Change to lat/lon
    lon0temp, lat0temp = grid['basemap'](xstart0, ystart0, inverse=True)


    ## Determine how many drifters to start based on the size of the transport at the start time ##

    # For this part, linearly interpolate the discharge rate data between outputs to the model
    # output timing, every 4 hours, so it can be selected by index easily
    im = find(rdates>=date)[0]-1 # index to the left
    ip = im + 1
    # interpolation weighting (from right)
    rintp = (netCDF.date2num(date, units)-netCDF.date2num(rdates[im], units))/ \
            (netCDF.date2num(rdates[ip], units)-netCDF.date2num(rdates[im], units))
    rintm = 1-rintp
    Qinterp = Q[im,:]*rintm + Q[ip,:]*rintp

    # Qinterp = np.interp(netCDF.date2num(date, units), netCDF.date2num(rdates, units), Q)
    # # Find index in river time
    # rind = find(netCDF.num2date(rt, runits)==date) # startdate should exactly equal a river time
    # find ndrifters based on the transports at that time for all rivers
    # ndrifters = abs(np.round(Q[rind,:]/5)[0])
    T0max = abs(Q).min()
    ndrifters = abs(np.round(Qinterp/T0max)) # about 3 m^3/s per drifter

    # make ndrifters number for each starting location
    lon0 = []; lat0 = []; T0 = []
    for iloc in xrange(ndrifters.size):
        lon0.extend(np.ones(ndrifters[iloc])*lon0temp[iloc])
        lat0.extend(np.ones(ndrifters[iloc])*lat0temp[iloc])
        T0.extend(np.ones(ndrifters[iloc])*(Qinterp[iloc]/ndrifters[iloc]))

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
