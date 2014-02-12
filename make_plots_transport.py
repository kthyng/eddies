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
mpl.rcParams.update({'font.size': 24})
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

# modifiers for pulling in files in order
fmods = ['2007-05', '2008-05', '2007-06', '2008-06',
		'2007-07', '2008-07', '2007-08', '2008-08']

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
fig = plt.figure(figsize=(17,15))

for i in xrange(len(S)):
	ax = fig.add_subplot(4,2,i+1)

	if i==0:
		ax.set_title('2007')
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0])

	elif i==1:
		ax.set_title('2008')
		ax.yaxis.set_label_position("right")
		ax.set_ylabel('May')
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0], parslabels=[0,0,0,0])

	if i==2:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0])

	elif i==3:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0], parslabels=[0,0,0,0])
		ax.yaxis.set_label_position("right")
		ax.set_ylabel('June')

	elif i==4:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0])

	elif i==5:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), parslabels=[0,0,0,0], merslabels=[0,0,0,0])
		ax.yaxis.set_label_position("right")
		ax.set_ylabel('July')

	elif i==6:
		tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))

	elif i==7:
		tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
		ax.yaxis.set_label_position("right")
		ax.set_ylabel('August')

	mappable = ax.pcolormesh(xr, yr, S[i]/Smax, cmap='Blues', vmax=0.1)

# Colorbar in upper left corner
cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
cb.set_label('Surface transport', fontsize=24)
cb.ax.tick_params(labelsize=20) 

fig.savefig('figures/transport/all.png', bbox_inches='tight', dpi=100)
