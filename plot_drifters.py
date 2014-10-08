'''
Plot a set of drifters in map view, in time.
'''

import matplotlib.pyplot as plt
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from matplotlib.mlab import find
import pdb
import numpy as np
import matplotlib as mpl
import os
import glob

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

# read in grid
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# loc = '/home/kthyng/shelf/grid.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True, llcrnrlat=27)

# Read in drifter tracks
dd = 1 #500 # drifter decimation
# Recommended options: 
# '2007-05-30T12', '2007-05-15T12', '2007-05-01T12'
# '2008-05-30T12', '2008-05-15T12', '2008-05-01T12'
# Read in groups of start dates of simulations for better representation
startdates = '2008-05-??T00'
files = glob.glob('tracks/' + startdates + 'gc.nc')
d = netCDF.MFDataset(files)
xg = d.variables['xg'][::dd,:]
yg = d.variables['yg'][::dd,:]
ind = (xg == -1)
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
xp[ind] = np.nan; yp[ind] = np.nan
tp = d.variables['tp'][:]
# d.close()

color = (0.03137255,  0.30191466,  0.58840448) # the almost dark blue from Blues

# Plot drifters, starting 5 days into simulation
# 2 days for one tail, 3 days for other tail
# t = tp-tp[0]
days = (tp-tp[0])/(3600.*24)
dates = netCDF.num2date(tp, d.variables['tp'].units)
# Find indices relative to present time
i5daysago = 0 # keeps track of index 5 days ago
i1dayago = find(days>=4)[0] # index for 1 day ago, which starts as 4 days in
i2daysago = find(days>=3)[0] # index for 2 days ago, which starts as 3 days in
i5days = find(days>=5)[0] # index for 5 days in
nt = tp.size # total number of time indices
# for i in np.arange(0,nt+1,5):
for i in np.arange(i5days,nt+1,5):

    fname = 'figures/drifters/' + startdate + '/' + dates[i].isoformat()[:-6] + '.png'

    if os.path.exists(fname):
        # Still need to update indices
        i5daysago += 5
        i2daysago += 5
        i1dayago += 5
        continue

    # pdb.set_trace()
    
    # Plot background
    fig = plt.figure(figsize=(18,6.6))
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax)

    # Plot 5 days ago to 2 days ago
    ax.plot(xp[:,i5daysago:i2daysago].T, yp[:,i5daysago:i2daysago].T, color=color, lw=2, alpha=0.4)

    # Plot 1-2 day tail
    ax.plot(xp[:,i2daysago:i1dayago].T, yp[:,i2daysago:i1dayago].T, color=color, lw=2, alpha=0.7)

    # Plot 0-1 day tail
    ax.plot(xp[:,i1dayago:i].T, yp[:,i1dayago:i].T, color=color, lw=3, alpha=0.9)

    # Plot drifter locations
    ax.plot(xp[:,i].T, yp[:,i].T, 'o', color=color, ms=10)

    # Time
    ax.text(0.075, 0.85, dates[i].isoformat()[:-6], transform=ax.transAxes, fontsize=20)

    # Drifter legend
    ax.plot(0.0895, 0.8, 'o', color=color, ms=10, transform=ax.transAxes) # drifter head
    ax.plot([0.075, 0.1], [0.76, 0.76], color=color, lw=3, alpha=0.9, transform=ax.transAxes) # drifter tail #1
    ax.plot([0.075, 0.1], [0.72, 0.72], color=color, lw=2, alpha=0.7, transform=ax.transAxes) # drifter tail #2
    ax.plot([0.075, 0.1], [0.68, 0.68], color=color, lw=2, alpha=0.4, transform=ax.transAxes) # drifter tail #3
    ax.text(0.125, 0.79, 'Drifter location', color=color, transform=ax.transAxes, fontsize=16)
    ax.text(0.125, 0.75, '1 day prior', color=color, alpha=0.9, transform=ax.transAxes, fontsize=16)
    ax.text(0.125, 0.71, '2 days prior', color=color, alpha=0.7, transform=ax.transAxes, fontsize=16)
    ax.text(0.125, 0.67, '5 days prior', color=color, alpha=0.4, transform=ax.transAxes, fontsize=16)
    # ax.plot(0.0895, 0.9, 'or', ms=10, transform=ax.transAxes) # drifter head
    # ax.plot([0.075, 0.1], [0.875, 0.875], 'darkcyan', lw=3, alpha=0.9, transform=ax.transAxes) # drifter tail #1
    # ax.plot([0.075, 0.1], [0.85, 0.85], 'darkcyan', lw=2, alpha=0.7, transform=ax.transAxes) # drifter tail #2
    # ax.plot([0.075, 0.1], [0.825, 0.825], 'darkcyan', lw=2, alpha=0.4, transform=ax.transAxes) # drifter tail #3
    # ax.text(0.125, 0.89, 'Drifter location', color='r', transform=ax.transAxes, fontsize=16)
    # ax.text(0.125, 0.866, '1 day prior', color='darkcyan', alpha=0.9, transform=ax.transAxes, fontsize=16)
    # ax.text(0.125, 0.842, '2 days prior', color='darkcyan', alpha=0.7, transform=ax.transAxes, fontsize=16)
    # ax.text(0.125, 0.818, '5 days prior', color='darkcyan', alpha=0.4, transform=ax.transAxes, fontsize=16)

    # Update indices
    i5daysago += 5
    i2daysago += 5
    i1dayago += 5

    # plt.show()
    # pdb.set_trace()

    fig.savefig(fname, bbox_inches='tight')

    plt.close()
