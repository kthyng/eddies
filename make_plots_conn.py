'''
Script to do analysis for where drifters from the river inputs go.

Run from eddies directory.
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
from glob import glob
from matplotlib.mlab import find, Path


## ---------- Miss vs. Atch ---------- ##



## ---------- Miss vs. Atch ---------- ##




## ---------- Year ---------- ##

mpl.rcParams.update({'font.size': 22})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# Find files to run through
years = np.arange(2007,2009)
Files = []
for year in years:
    Files.append(glob('tracks/' + str(year) + '*gc.nc'))

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, llcrnrlat=27.01, 
        urcrnrlat=30.5, llcrnrlon=-97.8, urcrnrlon=-87.7)

# Time vector for first simulation, then 
d = netCDF.Dataset(Files[0][0])
tp = d.variables['tp'][:]
d.close()
days = (tp-tp[0])/(3600.*24)

# Different advection times
taus = np.array([1,10,20,30,60,90]) # days
# Use histogram for sample space
# H = np.zeros((taus.size,100,100))
for j, year in enumerate(years):

    d = netCDF.MFDataset(Files[j], aggdim='ntrac')

    for i, tau in enumerate(taus):
        tind = find(days==tau)
        xg = d.variables['xg'][:,tind]
        yg = d.variables['yg'][:,tind]

        # ntrac = xg.size
        # weights = (1./ntrac)*np.ones(ntrac)

        nanind = np.isnan(xg) + (xg==-1)
        nnanind = ~nanind
        xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
        xp = xp[nnanind]; yp = yp[nnanind]
        # xp[nanind] = np.nan; yp[nanind] = np.nan
        del(xg,yg) # don't need grid info anymore

        # # add up locations of drifters
        # Htemp, xe, ye = np.histogram2d(xp[:,0], yp[:,0], bins=(100,100), weights=weights,
        #                         range=[[grid['xpsi'].min(), grid['xpsi'].max()],
        #                                 [grid['ypsi'].min(), grid['ypsi'].max()]])

        # H[i,:,:] = Htemp

        # Just directly plot with hexbin
        fig = plt.figure(figsize=(17,9))
        ax = fig.add_subplot(111)
        tracpy.plotting.hist(xp, yp, isll=False, fname='conn/' + str(year) + 'day' + str(tau), tind='other', 
                            which='hexbin', bins=(40,40), grid=grid, vmax=8,
                            fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
        plt.close()
        # fig.show()


    d.close()