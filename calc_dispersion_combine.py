'''
Combine already-calculated dispersions in various ways.

Things to investigate:

* Compare Mississippi and Atchafalaya
* Compare months
* Compare years

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


## -------- River Comparison ---------- ##

Files = glob('calcs/200*.npz') # ALL the files

# Which indices are from Miss and which from Atch
imiss = np.arange(0,43)
iatch = np.arange(43,50)

# Have to go through every file and combine the D2's for Atch and Miss separately
D2miss = np.zeros(2701); D2atch = np.zeros(2701)
nnansmiss = np.zeros(2701); nnansatch = np.zeros(2701)
for File in Files:

    d = np.load(File)
    D2 = d['D2']; nnans = d['nnans']; t = d['t']
    d.close()

    # Add up after un-averaging
    nanind = np.isnan(D2)
    D2[nanind] = 0
    D2miss += np.nansum(D2[imiss,:]*nnans[imiss,:], axis=0)
    nnansmiss += nnans[imiss,:].sum(axis=0)
    D2atch += np.nansum(D2[iatch,:]*nnans[iatch,:], axis=0)
    nnansatch += nnans[iatch,:].sum(axis=0)

# re-average now, all together
D2miss /= nnansmiss
D2atch /= nnansatch

# Save results
np.savez('calcs/riversD2.npz', D2miss=D2miss, D2atch=D2atch, nnansmiss=nnansmiss, nnansatch=nnansatch, t=t)

## --------------------------------- ##



## -------- Month Comparison ---------- ##

months = np.arange(5,9)

D2 = np.zeros((months.size,2701)); 
nnans = np.zeros((months.size,2701)); 

for i,month in enumerate(months):

    Files = glob('calcs/*-' + str(month).zfill(2) + '-*.npz') # files by month

    for File in Files:

        d = np.load(File)
        D2temp = d['D2']; nnanstemp = d['nnans']; t = d['t']
        d.close()

        # Add up after un-averaging
        nanind = np.isnan(D2temp)
        D2temp[nanind] = 0
        D2[i,:] += np.nansum(D2temp*nnanstemp, axis=0)
        nnans[i,:] += nnanstemp.sum(axis=0)


# re-average now, all together
D2 /= nnans

# Save results
np.savez('calcs/monthsD2.npz', D2=D2, nnans=nnans, t=t)


## --------------------------------- ##



## -------- Year Comparison ---------- ##

years = np.arange(2007,2009)

D2 = np.zeros((years.size,2701)); 
nnans = np.zeros((years.size,2701)); 

for i,year in enumerate(years):

    Files = glob('calcs/' + str(year) + '-*.npz') # files by year

    for File in Files:

        d = np.load(File)
        D2temp = d['D2']; nnanstemp = d['nnans']; t = d['t']
        d.close()

        # Add up after un-averaging
        nanind = np.isnan(D2temp)
        D2temp[nanind] = 0
        D2[i,:] += np.nansum(D2temp*nnanstemp, axis=0)
        nnans[i,:] += nnanstemp.sum(axis=0)


# re-average now, all together
D2 /= nnans

# Save results
np.savez('calcs/yearsD2.npz', D2=D2, nnans=nnans, t=t)


## --------------------------------- ##



## ---------- Everything separately ----------- ##

months = np.arange(5,9)
years = np.arange(2007,2009)

# Which indices are from Miss and which from Atch
imiss = np.arange(0,43)
iatch = np.arange(43,50)

D2 = np.zeros((months.size, years.size, 2, 2701)); 
nnans = np.zeros((months.size, years.size, 2, 2701)); 

for j, month in enumerate(months):
    for i, year in enumerate(years):

        Files = glob('calcs/' + str(year) + '-' + str(month).zfill(2) + '-*.npz') # files by year and month

        for File in Files:

            d = np.load(File)
            D2temp = d['D2']; nnanstemp = d['nnans']; t = d['t']
            d.close()

            # Add up after un-averaging
            nanind = np.isnan(D2temp)
            D2temp[nanind] = 0
            # miss
            D2[j,i,0,:] += np.nansum(D2temp[imiss,:]*nnanstemp[imiss,:], axis=0)
            nnans[j,i,0,:] += nnanstemp[imiss,:].sum(axis=0)
            # atch
            D2[j,i,1,:] += np.nansum(D2temp[iatch,:]*nnanstemp[iatch,:], axis=0)
            nnans[j,i,1,:] += nnanstemp[iatch,:].sum(axis=0)

# re-average now, all together
D2 /= nnans

# Save results
np.savez('calcs/allD2.npz', D2=D2, nnans=nnans, t=t)
