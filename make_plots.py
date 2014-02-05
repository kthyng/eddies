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
d = netCDF.Dataset(loc)

# Model time period to use
starttime = netCDF.date2num(datetime(2007, 5, 1, 12, 0, 0), d.variables['ocean_time'].units)
endtime = netCDF.date2num(datetime(2007, 9, 1, 12, 0, 0), d.variables['ocean_time'].units)
dt = d.variables['ocean_time'][1] - d.variables['ocean_time'][0] # 4 hours in seconds
ts = np.arange(starttime, endtime, dt)
itshift = find(starttime==d.variables['ocean_time'][:]) # shift to get to the right place in model output
##

## Colormap ##
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.4, 1.0, 0.7),
                 (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),
                   (0.4, 1.0, 0.0),
                   (1.0, 1.0, 1.0)),
         'blue': ((0.0, 0.0, 0.0),
                  (0.4, 1.0, 0.0),
                  (1.0, 0.5, 1.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
pcolor(rand(10,10),cmap=my_cmap)
colorbar()

np.array([0, 13, 21, 26, 30, 32, 34, 35, 35.5, 36])

cdict = {'red': ((0/36., 0.0, 0.1),
                 (13/36., 0.1, 0.2),
                 (21/36., 0.2, 0.3),
                 (26/36., 0.3, 0.4),
                 (30/36., 0.4, 0.5),
                 (32/36., 0.5, 0.6),
                 (34/36., 0.6, 0.7),
                 (35/36., 0.7, 0.8),
                 (35.5/36., 0.8, 0.9),
                 (36/36., 0.9, 1.0)),
		'green': ((0/36., 0.0, 0.1),
                 (13/36., 0.1, 0.2),
                 (21/36., 0.2, 0.3),
                 (26/36., 0.3, 0.4),
                 (30/36., 0.4, 0.5),
                 (32/36., 0.5, 0.6),
                 (34/36., 0.6, 0.7),
                 (35/36., 0.7, 0.8),
                 (35.5/36., 0.8, 0.9),
                 (36/36., 0.9, 1.0)),
		'blue': ((0/36., 1, 1),
                 (13/36., 1, 1),
                 (21/36., 1, 1),
                 (26/36., 1, 1),
                 (30/36., 1, 1),
                 (32/36., 1, 1),
                 (34/36., 1, 1),
                 (35/36., 1, 1),
                 (35.5/36., 1, 1),
                 (36/36., 1, 1.0))}


# Loop through time
for it,t in enumerate(ts):

	# Set up plot
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)
tracpy.plotting.background(grid=grid, ax=ax)

# Plot surface salinity
salt = np.squeeze(d.variables['salt'][it+itshift,-1,:,:])
plt.contourf(xr, yr, salt, levels=np.array([ 1.        ,  12.82424751,  20.76480532,  26.09727644,
        29.67829038,  32.08311557,  33.69807275,  34.78259511,
        35.51090468,  36.]),
			cmap='Blues')#, vmin=5, vmax=35)

# Colorbar in upper left corner
cax = fig.add_axes([0.15, 0.75, 0.3, 0.03]) #colorbar axes
cb = colorbar(cax=cax, orientation='horizontal')
cb.set_label('Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=20)
cb.ax.tick_params(labelsize=20) 

	# Read in all drifter sets starting before or up to the current datetime
	for File for Files:
		d = netCDF.Dataset(File)
		xg = d.variables['xg'][:]
		yg = d.variables['yg'][:]
		tg = d.variables['tp'][:]
		d.close()

		# Change to projected drifter locations now
		nanind = np.isnan(xg) # indices where nans are location in xg, yg; for reinstitution of nans
		xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
		xp[nanind] = np.nan; yp[nanind] = np.nan
		del(xg,yg) # don't need grid info anymore

	# Plot drifters
	plot(xp[:,tind], yp[:,tind], 'o', color='darkcyan')

	# ax.plot(xp[::5,tind], yp[::5,tind], 'o', color='darkcyan')
	# plt.savefig('figures/2010-07-01T00/' + str(tind) + '.png', dpi=150)
	# plt.close(fig)
