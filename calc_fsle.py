'''
Calculate FSLE of drifters starting from the same locations.
'''

import numpy as np
import pdb
from matplotlib.mlab import find
import netCDF4 as netCDF
from scipy import ndimage
import time
from glob import glob
import tracpy
import os

def get_dist(lon1, lons, lat1, lats): 
    '''
    Function to compute great circle distance between point lat1 and lon1 
    and arrays of points given by lons, lats or both same length arrays.
    Uses Haversine formula. Distance is in km.
    '''

    lon1 = lon1*np.pi/180.
    lons = lons*np.pi/180.
    lat1 = lat1*np.pi/180.
    lats = lats*np.pi/180.

    earth_radius = 6373.
    distance = earth_radius*2.0*np.arcsin(np.sqrt(np.sin(0.50*(lat1-lats))**2 \
                                       + np.cos(lat1)*np.cos(lats) \
                                       * np.sin(0.50*(lon1-lons))**2))
    return distance

def calc_fsle(lonpc, latpc, lonp, latp, tp, alpha=np.sqrt(2)):
    '''
    Calculate the relative dispersion of tracks lonp, latp as directly compared with
    the tracks described by lonpc, latpc. The two sets of tracks must start in the same 
    locations since this is assumed for making "pairs" of drifters for comparison (and 
    therefore pairs do not need to be found). The relative dispersion in this case is a
    measure of the difference between the two simulations, and is aimed at being used
    for examining differences in tracks due to changes in the numerical simulation.
    The tracks should also be coincident in time, but the script will find a way to match
    them up for the overlap periods.


    Inputs:
        lonpc, latpc    Longitude/latitude of the control drifter tracks [ndrifter,ntime]
        lonp, latp      Longitude/latitude of the drifter tracks [ndrifter,ntime]
        squared         Whether to present the results as separation distance squared or 
                        not squared. Squared by default.

    Outputs:
        D2              Relative dispersion (squared or not) averaged over drifter 
                        pairs [ntime]. In km (squared or not).
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of relative dispersion, first multiply by nnans, then
    combine with other relative dispersion calculations, then divide by the total number
    of nnans.

    Example call:
    dc = netCDF.Dataset('tracks/tseas_use300_nsteps1.nc') (5 min, control output)
    d = netCDF.Dataset('tracks/tseas_use1200_nsteps1.nc') (20 min, comparison output)
    tracpy.calcs.rel_dispersion_comp(dc.variables['lonp'][:], dc.variables['latp'][:], dc.variables['tp'][:],
                                     d.variables['lonp'][:], d.variables['latp'][:], d.variables['tp'][:],
                                     squared=True)
    '''

    # We know that drifters from the two sets have a one to one correspondence
    dist = get_dist(lonpc, lonp, latpc, latp)

    Rs = np.asarray([0.1*alpha**i for i in np.arange(28)])

    # Find first time dist>delta and dist>delta*alpha for each delta to
    # then linearly interpolate to find the corresponding time
    # FOR ONE DRIFTER TO START AND ONE DELTA
    tau = np.zeros(Rs.size)
    nnans = np.zeros(Rs.size) # not nans
    for idrifter in xrange(dist.shape[0]):
    # idrifter = 0
        # delta = Rs[10]
        for i, R in enumerate(Rs):

            if R<=np.nanmax(dist[idrifter,:]) \
                and Rs[i+1]<=np.nanmax(dist[idrifter,:]) \
                and R>=dist[idrifter,0]:

                # for delta
                ind = find(dist[idrifter,:]>=R)[0]
                time1 = np.interp(R, dist[idrifter, ind-1:ind+1], tp[ind-1:ind+1])
                # print R, dist[idrifter,ind-1:ind+1]

                # for delta*alpha
                ind = find(dist[idrifter,:]>=Rs[i+1])[0]
                time2 = np.interp(Rs[i+1], dist[idrifter, ind-1:ind+1], tp[ind-1:ind+1])
                # print Rs[i+1], dist[idrifter,ind-1:ind+1]

                dt = time2-time1

            else:
                dt = np.nan

            if not np.isnan(dt):
                # print R, dt
                tau[i] += dt
                nnans[i] += 1 # counting not-nan entries for averaging later

    return tau, nnans, Rs


def run():
    '''
    Run FSLE calculation for shelf eddy drifter simulations.
    '''

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc)

    Files = glob('tracks/2008-*gc.nc')

    for File in Files:

        fname = 'calcs/' + File[:-5].split('/')[-1] + 'fsle.npz'

        if os.path.exists(fname): # don't redo if already done
            continue

        d = netCDF.Dataset(File)
        xg = d.variables['xg'][:]
        yg = d.variables['yg'][:]
        tp = d.variables['tp'][:]
        d.close()

        nanind = np.isnan(xg) + (xg==-1) + (np.ceil(xg)<=5) + (np.ceil(xg)>=grid['xr'].shape[0]-5) + (np.ceil(yg)<=5)
        lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll') 
        lonp[nanind] = np.nan; latp[nanind] = np.nan
        del(xg,yg) # don't need grid info anymore

        # Loop through the unique starting locations for the drifters, 
        # to calculate dispersion based on drifters starting in same place
        # Need to sub-loop through every other drifter starting at that location
        # to get all of the drifter pairs

        # Find unique initial locations
        _, iunique = np.unique(lonp[:,0], return_index=True)
        iunique = np.sort(iunique)

        # Save altogether for a single simulation
        fsle = np.zeros((iunique.size, 28))
        nnans = np.zeros((iunique.size, 28))
        ntrac = lonp.shape[0] # num drifters

        for i, istartloc in enumerate(iunique): # loop through the unique starting location indices

            if i==iunique.size-1: # if last unique starting position index
                iendloc = lonp.shape[0]
            else: 
                iendloc = iunique[i+1]

            for j in xrange(ntrac-1): # loop over drifters
                # pdb.set_trace()
                fsletemp, nnanstemp, Rs = calc_fsle(lonp[j,:], latp[j,:], 
                                            lonp[j+1:iendloc,:], latp[j+1:iendloc,:], tp)
                # pdb.set_trace()
                fsle[i,:] += fsletemp
                nnans[i,:] += nnanstemp

                # NOT Now average all pairs starting at this unique location
                # fsle = fsle/nnans
        pdb.set_trace()
        # save: fsle in time, averaged over all combinations of drifters starting at
        # a unique river input point for a unique starting time
        np.savez(fname, fsle=fsle, nnans=nnans, Rs=Rs)



if __name__ == "__main__":
    run()    
