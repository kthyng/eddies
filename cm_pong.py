'''
PONG colormaps.
'''

import numpy as np
import matplotlib
from pylab import *

def salinity(cmap, levels):
    '''
    Colormap for salinity, with bigger chunks of salinity per color
    section at lower salinity than higher.
    From http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps

    Inputs:
        cmap        Colormap name to use
        levels      edges of colors, as in contourf

    Outputs:
        my_cmap     colormap instance
    '''

    N = levels.size

    # Colors on either side of the edges
    rgb0 = cm.get_cmap(cmap)(linspace(0.0, 0.9, N))[:,0:3]
    rgb1 = cm.get_cmap(cmap)(linspace(0.1, 1.0, N))[:,0:3]

    red = np.vstack((levels/levels.max(), 
                    rgb0[:,0], 
                    rgb1[:,0])).T
    red = tuple(map(tuple, red))

    green = np.vstack((levels/levels.max(), 
                    rgb0[:,1], 
                    rgb1[:,1])).T
    green = tuple(map(tuple, green))

    blue = np.vstack((levels/levels.max(), 
                    rgb0[:,2], 
                    rgb1[:,2])).T
    blue = tuple(map(tuple, blue))

    cdict = {'red':red, 'green':green, 'blue':blue}

    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

    return my_cmap

def test(cmap):
    '''
    Test colormap.
    '''

    figure()
    pcolor(rand(10,10), cmap=cmap)
    cb = colorbar()

    # how to set the labels to reflect the stretching
