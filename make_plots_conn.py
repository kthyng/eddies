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


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# 'rivers', 'rivers1year' 'months+years', 'years'
whichplot = 'months+years' 


## ---------- Miss vs. Atch ---------- ##

if whichplot=='rivers':

    imiss = np.arange(0,44)
    iatch = np.arange(43,50)

    # Find files to run through
    # years = np.arange(2007,2009)
    # Files = []
    # for year in years:
    Files = glob('tracks/*gc.nc')

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, llcrnrlat=27.01, 
            urcrnrlat=30.5, llcrnrlon=-97.8, urcrnrlon=-87.7)

    # Time vector for first simulation, then 
    d = netCDF.Dataset(Files[0])
    tp = d.variables['tp'][:]
    d.close()
    days = (tp-tp[0])/(3600.*24)

    # Different advection times
    taus = np.array([1,10,20,30,60,90]) # days

    for i, tau in enumerate(taus):

        xpsavemiss = []; ypsavemiss = []
        xpsaveatch = []; ypsaveatch = []

        # Run through years just because we can't open all files at once
        for File in Files: # in enumerate(years):

            d = netCDF.Dataset(File, aggdim='ntrac')

            tind = find(days==tau)
            xg = d.variables['xg'][:,tind]
            yg = d.variables['yg'][:,tind]

            # initial locations, for finding starting locations
            xg0 = d.variables['xg'][:,0]

            nanind = np.isnan(xg) + (xg==-1)
            nnanind = ~nanind
            xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
            xp = xp[nnanind]; yp = yp[nnanind]
            del(xg,yg) # don't need grid info anymore

            # Find unique initial locations for different inputs
            _, iunique = np.unique(xg0, return_index=True)
            iunique = np.sort(iunique)

            # save up drifter locations from all files for the rivers separately and taus separately

            ## Mississippi ##
            xpsavemiss.extend(xp[iunique[imiss[0]]:iunique[imiss[-1]]])
            ypsavemiss.extend(yp[iunique[imiss[0]]:iunique[imiss[-1]]])

            ## Atchafalaya ##
            xpsaveatch.extend(xp[iunique[iatch[0]]:])
            ypsaveatch.extend(yp[iunique[iatch[0]]:])

            d.close()


        # Just directly plot with hexbin
        ## Miss ##
        fig = plt.figure(figsize=(17,9))
        ax = fig.add_subplot(111)
        ax.set_frame_on(False) # kind of like it without the box
        tracpy.plotting.hist(np.asarray(xpsavemiss), np.asarray(ypsavemiss), isll=False, fname='conn/miss' + 'day' + str(tau), tind='other', 
                            which='hexbin', bins=(40,40), grid=grid, binscale='log', vmax=0.3, 
                            fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
        plt.close()

        ## Atch ##
        fig = plt.figure(figsize=(17,9))
        ax = fig.add_subplot(111)
        ax.set_frame_on(False) # kind of like it without the box
        tracpy.plotting.hist(np.asarray(xpsaveatch), np.asarray(ypsaveatch), isll=False, fname='conn/atch' + 'day' + str(tau), tind='other', 
                            which='hexbin', bins=(40,40), grid=grid, binscale='log', vmax=0.3, 
                            fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
        plt.close()



## ---------- Miss vs. Atch 1 year ---------- ##

elif whichplot=='rivers1year':

    imiss = np.arange(0,44)
    iatch = np.arange(43,50)

    year = 2008

    # Find files to run through
    # years = np.arange(2007,2009)
    # Files = []
    # for year in years:
    Files = glob('tracks/' + str(year) + '*gc.nc')

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, llcrnrlat=27.01, 
            urcrnrlat=30.5, llcrnrlon=-97.8, urcrnrlon=-87.7)

    # Time vector for first simulation, then 
    d = netCDF.Dataset(Files[0])
    tp = d.variables['tp'][:]
    d.close()
    days = (tp-tp[0])/(3600.*24)

    # Different advection times
    taus = np.array([1,10,20,30,60,90]) # days

    for i, tau in enumerate(taus):

        xpsavemiss = []; ypsavemiss = []
        xpsaveatch = []; ypsaveatch = []

        # Run through years just because we can't open all files at once
        for File in Files: # in enumerate(years):

            d = netCDF.Dataset(File, aggdim='ntrac')

            tind = find(days==tau)
            xg = d.variables['xg'][:,tind]
            yg = d.variables['yg'][:,tind]

            # initial locations, for finding starting locations
            xg0 = d.variables['xg'][:,0]

            nanind = np.isnan(xg) + (xg==-1)
            nnanind = ~nanind
            xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
            xp = xp[nnanind]; yp = yp[nnanind]
            del(xg,yg) # don't need grid info anymore

            # Find unique initial locations for different inputs
            _, iunique = np.unique(xg0, return_index=True)
            iunique = np.sort(iunique)

            # save up drifter locations from all files for the rivers separately and taus separately

            ## Mississippi ##
            xpsavemiss.extend(xp[iunique[imiss[0]]:iunique[imiss[-1]]])
            ypsavemiss.extend(yp[iunique[imiss[0]]:iunique[imiss[-1]]])

            ## Atchafalaya ##
            xpsaveatch.extend(xp[iunique[iatch[0]]:])
            ypsaveatch.extend(yp[iunique[iatch[0]]:])

            d.close()


        # Just directly plot with hexbin
        ## Miss ##
        fig = plt.figure(figsize=(17,9))
        ax = fig.add_subplot(111)
        ax.set_frame_on(False) # kind of like it without the box
        tracpy.plotting.hist(np.asarray(xpsavemiss), np.asarray(ypsavemiss), isll=False, 
                            fname='conn/miss' + str(year) + 'day' + str(tau), tind='other', 
                            which='hexbin', bins=(40,40), grid=grid, binscale='log', vmax=0.3, 
                            fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
        plt.close()

        ## Atch ##
        fig = plt.figure(figsize=(17,9))
        ax = fig.add_subplot(111)
        ax.set_frame_on(False) # kind of like it without the box
        tracpy.plotting.hist(np.asarray(xpsaveatch), np.asarray(ypsaveatch), isll=False, 
                            fname='conn/atch' + str(year) + 'day' + str(tau), tind='other', 
                            which='hexbin', bins=(40,40), grid=grid, binscale='log', vmax=0.3, 
                            fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
        plt.close()



## ---------- Months and Years ---------- ##

elif whichplot=='months+years':

    # Find files to run through
    years = np.arange(2007,2009)
    months = np.arange(5,9)

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, llcrnrlat=27.01, 
            urcrnrlat=30.5, llcrnrlon=-97.8, urcrnrlon=-87.7)

    # Time vector for first simulation, then 
    d = netCDF.Dataset('tracks/2007-05-01T12gc.nc')
    tp = d.variables['tp'][:]
    d.close()
    days = (tp-tp[0])/(3600.*24)

    # Different advection times
    taus = np.array([1,10,20,30,60,90]) # days
    # Use histogram for sample space
    # H = np.zeros((taus.size,100,100))
    for j, year in enumerate(years):

        for month in months:

            # aggregate files from the same month/year
            Files = glob('tracks/' + str(year) + '-' + str(month).zfill(2) + '*gc.nc')

            d = netCDF.MFDataset(Files, aggdim='ntrac')

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
                fig = plt.figure(figsize=(12,8))
                ax = fig.add_subplot(111)
                ax.set_frame_on(False) # kind of like it without the box
                tracpy.plotting.hist(xp, yp, isll=False, fname='conn/' + str(year) + str(month).zfill(2) + 'day' + str(tau), tind='other', 
                                    which='hexbin', bins=(40,40), grid=grid, binscale='log', vmax=0.3, 
                                    fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
                plt.close()
                # fig.show()


            d.close()



## ---------- Year ---------- ##

elif whichplot=='years':

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
            ax.set_frame_on(False) # kind of like it without the box
            tracpy.plotting.hist(xp, yp, isll=False, fname='conn/' + str(year) + 'day' + str(tau), tind='other', 
                                which='hexbin', bins=(40,40), grid=grid, binscale='log', vmax=0.3, 
                                fig=fig, ax=ax, Label='Drifter location (%%), day=%i' % tau)
            plt.close()
            # fig.show()


        d.close()
