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

## Drifters ##
# Tracks to plot
Files = glob.glob('tracks/2007' + '*gc.nc')

# # Read in info
# d = netCDF.MFDataset(Files, aggdim='ntrac')
# xg = d.variables['xg'][:]
# yg = d.variables['yg'][:]
# tg = d.variables['tp'][:]
# d.close()

# # Change to projected drifter locations now
# nanind = np.isnan(xg) # indices where nans are location in xg, yg; for reinstitution of nans
# xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
# xp[nanind] = np.nan; yp[nanind] = np.nan
# del(xg,yg) # don't need grid info anymore
##

## Model output ##
m = netCDF.Dataset(loc)

# Model time period to use
starttime = netCDF.date2num(datetime(2007, 5, 1, 12, 0, 0), m.variables['ocean_time'].units)
endtime = netCDF.date2num(datetime(2007, 9, 1, 12, 0, 0), m.variables['ocean_time'].units)
dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
dates = netCDF.num2date(m.variables['ocean_time'][:], m.variables['ocean_time'].units)
months = ['January', 'February', 'March', 'April', 'May', 'June',
    'July', 'August', 'September', 'October', 'November', 'December']

# Colormap for model output
levels = (37-np.exp(linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity
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
for t in ts[5:6]:

    # Set up plot
    fig = plt.figure(figsize=(17,9))
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax, outline=False)

    itmodel = find(t==m.variables['ocean_time'][:])[0]

    # # Plot surface salinity
    # salt = np.squeeze(m.variables['salt'][itmodel,-1,:,:])
    # ax.contour(xr, yr, salt, [33], colors='k')
    # mappable = ax.pcolormesh(xr, yr, salt, cmap=cmap, vmin=0, vmax=35)

    # # Date
    # date = str(dates[itmodel].year) + ' ' + months[dates[itmodel].month-1] \
    #      + ' ' + str(dates[itmodel].day) + ' ' + str(dates[itmodel].hour).fill(2) + ':00'
    # text(0.75, 0.02, date, 
    #         fontsize=24, color='0.2', transform=ax.transAxes)

    # # Colorbar in upper left corner
    # cax = fig.add_axes([0.15, 0.75, 0.3, 0.03]) #colorbar axes
    # cb = colorbar(mappable, cax=cax, orientation='horizontal')
    # cb.set_label('Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=20)
    # cb.ax.tick_params(labelsize=18) 
    # cb.set_ticks(ticks)

    # Read in all drifter sets starting before or up to the current datetime
    for j,File in enumerate(Files[0:1]):

        d = netCDF.Dataset(File)
        xg = d.variables['xg'][:]
        yg = d.variables['yg'][:]
        tg = d.variables['tp'][:]
        d.close()

        if tg[0]>t: # this drifter simulation starts after this time
            continue
        else:
            itdrifter = find(tg==t)[0] # find index in time for drifters

        # Change to projected drifter locations now
        nanind = np.isnan(xg) # indices where nans are location in xg, yg; for reinstitution of nans
        xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
        xp[nanind] = np.nan; yp[nanind] = np.nan
        # del(xg,yg) # don't need grid info anymore

        # Plot drifters
        ax.plot(xp[:,:itdrifter].T, yp[:,:itdrifter].T, '-', color='0.5', zorder=7)
        ax.plot(xp[:,itdrifter], yp[:,itdrifter], 'o', color='0.5', zorder=7)

    # ax.plot(xp[::5,tind], yp[::5,tind], 'o', color='darkcyan')
    # plt.savefig('figures/2010-07-01T00/' + str(tind) + '.png', dpi=150)
    # plt.close(fig)
