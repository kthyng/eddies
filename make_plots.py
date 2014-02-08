'''
Run plots for drifters from these simulations.
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
import cm_pong
from matplotlib.mlab import find

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 26})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, llcrnrlat=27.01, 
        urcrnrlat=30.5, llcrnrlon=-97.8, urcrnrlon=-87.7)
xr = np.asanyarray(grid['xr'].T, order='C')
yr = np.asanyarray(grid['yr'].T, order='C')

## Model output ##
m = netCDF.Dataset(loc)

# Model time period to use
units = m.variables['ocean_time'].units
year = 2007
starttime = netCDF.date2num(datetime(year, 5, 1, 12, 0, 0), units)
endtime = netCDF.date2num(datetime(year, 9, 1, 12, 0, 0), units)
dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
dates = netCDF.num2date(m.variables['ocean_time'][:], units)
months = ['January', 'February', 'March', 'April', 'May', 'June',
    'July', 'August', 'September', 'October', 'November', 'December']

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity
cmap = cm_pong.salinity('YlGnBu', levels)
ilevels = [0,1,2,3,4,5,8] # which levels to label
ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
##

    # # Change axis and label color
    # ax.spines['bottom'].set_color('0.2')
    # ax.spines['top'].set_color('0.2')
    # ax.spines['left'].set_color('0.2')
    # ax.spines['right'].set_color('0.2')
    # ax.xaxis.label.set_color('0.2')
    # ax.yaxis.label.set_color('0.2')
    # ax.tick_params(axis='x', colors='0.2')
    # ax.tick_params(axis='y', colors='0.2')


# Loop through times that simulations were started
for t in ts[:29]:

    # Set up before plotting

    x14 = []; y14 = []; x7 = []; y7 = []; x = []; y = []

    itmodel = find(t==m.variables['ocean_time'][:])[0] # index for model output at this time

    ## Drifter file set up ##
    # Read in all drifter sets starting before or up to the current datetime
    lastsimdate = netCDF.num2date(t, units) # start date we're examining
    # Read in drifters from simulations that would reach this date; they run for 90 days
    firstsimdate = (netCDF.num2date(t, units)-timedelta(days=90))
    if firstsimdate < datetime(year, 5, 1, 12, 0, 0):
        firstsimdate = datetime(year, 5, 1, 12, 0, 0) # earliest simulation
    nsims = 1 + (lastsimdate-firstsimdate).total_seconds()/(3600*4.) # number of sims between start and end
    simdates = [firstsimdate + i*timedelta(hours=4) for i in xrange(int(nsims))]
    Files = []
    for i,simdate in enumerate(simdates):
        Files.append('tracks/' + simdates[i].isoformat()[0:13] + 'gc.nc')
    ##

    figname = 'figures/' + str(year) + '/' + lastsimdate.isoformat()[0:13] + '.png'

    # Don't redo plot
    if os.path.exists(figname):
        continue

    # Set up plot
    fig = plt.figure(figsize=(17,9))
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax, outline=False)

    # Date
    date = str(dates[itmodel].year) + ' ' + months[dates[itmodel].month-1] \
         + ' ' + str(dates[itmodel].day) + ' ' + str(dates[itmodel].hour).zfill(2) + ':00'
    ax.text(0.75, 0.02, date, fontsize=24, color='0.2', transform=ax.transAxes)

    # Plot surface salinity
    salt = np.squeeze(m.variables['salt'][itmodel,-1,:,:])
    # ax.contour(xr, yr, salt, [33], colors='k')
    mappable = ax.pcolormesh(xr, yr, salt, cmap=cmap, vmin=0, vmax=35)

    # Colorbar in upper left corner
    cax = fig.add_axes([0.15, 0.75, 0.3, 0.03]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=20)
    cb.ax.tick_params(labelsize=18) 
    cb.set_ticks(ticks)

    # I think I need to loop to use the right amount of time from each set of drifters
    for j,File in enumerate(Files):

        d = netCDF.Dataset(File)
        xg = d.variables['xg'][:]
        yg = d.variables['yg'][:]
        tg = d.variables['tp'][:]
        d.close()

        itdrifter = find(tg==t)[0] # find index in drifter time for drifters

        # Days back from this time
        days = (tg-tg[0])/(3600.*24)

        # Change to projected drifter locations now
        nanind = np.isnan(xg)*(xg==-1) # indices where nans are location in xg, yg; for reinstitution of nans
        xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
        xp[nanind] = np.nan; yp[nanind] = np.nan
        del(xg,yg) # don't need grid info anymore

        # Plot drifter tails for 7 days one color, 14 days another lighter color (find indices)
        # these are 0 if enough time hasn't passed, so works either way
        i7 = find(days>days[itdrifter]-7)[0]
        i14 = find(days>days[itdrifter]-14)[0]

        # Save drifter info and plot all at once so things don't overlap
        x14.append(xp[:,i14:itdrifter+1].T)
        y14.append(yp[:,i14:itdrifter+1].T)
        x7.append(xp[:,i7:itdrifter+1].T)
        y7.append(yp[:,i7:itdrifter+1].T)
        x.append(xp[:,itdrifter])
        y.append(yp[:,itdrifter])

    # Plot drifters
    for i in xrange(len(x14)):
        ax.plot(x14[i], y14[i], '-', color='0.7', zorder=7, linewidth=1, alpha=0.5)
    for i in xrange(len(x14)):
        ax.plot(x7[i], y7[i], '-', color='0.5', zorder=7, linewidth=2, alpha=0.7)
    for i in xrange(len(x14)):
        ax.plot(x[i], y[i], 'o', color='0.3', zorder=7)

    plt.savefig(figname, bbox_inches='tight', dpi=150)
    plt.close(fig)
