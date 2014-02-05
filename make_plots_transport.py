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
grid = tracpy.inout.readgrid(loc, llcrnrlat=27, urcrnrlat=30.5, llcrnrlon=-98)
xr = np.asanyarray(grid['xr'].T, order='C')
yr = np.asanyarray(grid['yr'].T, order='C')

Files = glob.glob('tracks/2008-08-*')
U = 0; V = 0
for File in Files:
	d = netCDF.Dataset(File)
	U += d.variables['U'][:]
	V += d.variables['V'][:]
	d.close



fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)
tracpy.plotting.background(grid=grid, ax=ax)

ax.pcolormesh(grid['xr'], grid['yr'], np.sqrt(op.resize(U,1)**2 + op.resize(V,0)**2), 
			cmap='Blues', vmax=30000)
ax.set_title('2008-08')
plt.show()

# ## Drifters ##
# # Tracks to plot
# Files = glob.glob('tracks/2007' + '*gc.nc')

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
# ##

# ## Model output ##
# d = netCDF.Dataset(loc)

# # Time period
# starttime = netCDF.date2num(datetime(2007, 5, 1, 12, 0, 0), d.variables['ocean_time'].units)
# endtime = netCDF.date2num(datetime(2007, 9, 1, 12, 0, 0), d.variables['ocean_time'].units)
# dt = d.variables['ocean_time'][1] - d.variables['ocean_time'][0] # 4 hours in seconds
# ts = np.arange(starttime, endtime, dt)
# itshift = find(starttime==d.variables['ocean_time'][:]) # shift to get to the right place in model output

# # Loop through time
# for it,t in enumerate(ts):

# 	# Set up plot
# 	fig = plt.figure(figsize=(20,10))
# 	ax = fig.add_subplot(111)
# 	tracpy.plotting.background(grid=grid, ax=ax)

# 	# Plot surface salinity
# 	pcolormesh(xr, yr, np.squeeze(d.variables['salt'][it+itshift,-1,:,:]), 
# 				cmap='Blues', vmin=5, vmax=35)

# 	# Colorbar in upper left corner
# 	cax = fig.add_axes([0.15, 0.75, 0.3, 0.03]) #colorbar axes
# 	cb = colorbar(cax=cax, orientation='horizontal')
# 	cb.set_label('Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=20)
# 	cb.ax.tick_params(labelsize=20) 

# 	# Plot drifters
# 	plot(xp[:,tind], yp[:,tind], 'o', color='darkcyan')

# 	# ax.plot(xp[::5,tind], yp[::5,tind], 'o', color='darkcyan')
# 	# plt.savefig('figures/2010-07-01T00/' + str(tind) + '.png', dpi=150)
# 	# plt.close(fig)
