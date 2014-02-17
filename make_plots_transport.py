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
from matplotlib.colors import LogNorm
import prettyplotlib as ppl

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 20})
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
grid = tracpy.inout.readgrid(loc, llcrnrlat=27.5, urcrnrlat=30.5, llcrnrlon=-93.9, urcrnrlon=-88.6)
# actually using psi grid here despite the name
xr = grid['xpsi']
yr = grid['ypsi']

# modifiers for pulling in files in order
fmods = ['2007-05', '2007-06', '2007-07', '2007-08', 
		'2008-05', '2008-06', '2008-07', '2008-08']

S = []
Smax = 0

for i,fmod in enumerate(fmods):

	## Read in transport info ##
	Files = glob.glob('tracks/' + fmod + '-*')

	U = 0; V = 0
	for File in Files:
		d = netCDF.Dataset(File)
		U += d.variables['U'][:]
		V += d.variables['V'][:]
		d.close

	# S is at cell centers, minus ghost points
	Stemp = np.sqrt(op.resize(U[:,1:-1],0)**2 + op.resize(V[1:-1,:],1)**2)
	S.append(Stemp)

	Smax = max((Smax,Stemp.max()))



## Plot ##
fig = plt.figure(figsize=(22,7))
fig.suptitle('Surface transport', fontsize=24)
for i in xrange(len(S)):

	ax = fig.add_subplot(2,4,i+1)

	if i==0:
		ax.set_title('May')
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0])

	elif i==1:
		ax.set_title('June')
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0], parslabels=[0,0,0,0])

	elif i==2:
		ax.set_title('July')
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0], parslabels=[0,0,0,0])

	elif i==3:
		ax.set_title('August')
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0], parslabels=[0,0,0,0])
		ax.yaxis.set_label_position("right")
		ax.set_ylabel('2007')

	elif i==4:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))

	elif i==5:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), parslabels=[0,0,0,0])

	elif i==6:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), parslabels=[0,0,0,0])

	elif i==7:
		tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
		ax.yaxis.set_label_position("right")
		ax.set_ylabel('2008')

	C =  np.log(S[i]/Smax)
	C = np.ma.masked_where(np.isinf(C), C)
	C = np.ma.masked_where(np.isnan(C), C)
	mappable = ax.pcolormesh(xr, yr, C, cmap='Blues', vmax=-1.5, vmin=-7)
	# mappable = ax.pcolormesh(xr, yr, S[i]/Smax, cmap='Blues', vmax=0.09)

# Colorbar in upper left corner
cax = fig.add_axes([0.375, 0.5, 0.3, 0.015]) #colorbar axes
cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
# cb.set_label('Surface transport', fontsize=20)
cb.ax.tick_params(labelsize=16) 
clim = cb.get_clim()
# cb.set_clim(exp(clim[0]), exp(clim[1]))
cb.set_ticks(np.linspace(clim[0],clim[1],6))
cb.set_ticklabels(exp(np.linspace(clim[0],clim[1],6)).round(3))
# cb.set_ticks(exp(np.linspace(clim[0],clim[1],8)))

fig.savefig('figures/transport/all_log.png', bbox_inches='tight', dpi=100)
