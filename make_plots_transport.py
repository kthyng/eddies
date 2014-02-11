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
import op

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
# actually using psi grid here despite the name
xr = grid['xpsi']
yr = grid['ypsi']


## Read in transport info ##
Files = glob.glob('tracks/2007-05-*')
U = 0; V = 0
for File in Files:
	d = netCDF.Dataset(File)
	U += d.variables['U'][:]
	V += d.variables['V'][:]
	d.close


## Do plot ##
fig = plt.figure(figsize=(17,9))
ax = fig.add_subplot(111)
tracpy.plotting.background(grid=grid, ax=ax)

# S is at cell centers, minus ghost points
S = np.sqrt(op.resize(U[:,1:-1],0)**2 + op.resize(V[1:-1,:],1)**2)
mappable = ax.pcolormesh(xr, yr, S/S.max(), cmap='Blues', vmax=0.1)

# Colorbar in upper left corner
cax = fig.add_axes([0.15, 0.75, 0.3, 0.03]) #colorbar axes
cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
cb.set_label('Surface transport', fontsize=20)
cb.ax.tick_params(labelsize=18) 
# cb.set_ticks(ticks)

ax.set_title('2007-05')
plt.show()

