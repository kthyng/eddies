#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py
Created on 2008-04-24

Copyright (c) 2008 Robert D. Hetland
"""
__docformat__ = "restructuredtext en"

import os
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import netcdftime

import octant

time_units = 'days since 1970-01-01 00:00:00'
utime = netcdftime.utime(time_units)

river_data = np.load('tarbert_landing2010.npy')
#river_atchafalaya = np.load('Simmesport.npy')
river_brazos = np.load('brazos2010.npy')
river_nueces = np.load('nueces2010.npy')
river_sanantonio = np.load('sanantonio2010.npy')
river_lavaca = np.load('lavaca2010.npy')
river_trinity = np.load('trinity2010.npy')
river_sabine = np.load('sabine2010.npy')
river_calcasieu = np.load('calcasieu2010.npy')

river_jd = utime.date2num(river_data['date'])

river = np.arange(1, 59)

river_Xposition = \
    np.array([ 574., 574., 575., 579., 581., 582.,  
               582., 587., 583., 581., 580.,
               583., 584., 585., 586., 587., 589.,590., 
               591., 592., 593., 594., 595., 596., 
               604., 604., 605., 603., 
               602., 599., 600., 599., 
               601., 601., 601, 600., 601., 601., 
               601., 600., 600., 597., 596., 
               442., 443., 444.,
               454., 455., 456., 457., 458.,
               338., 316., 284., 242., 202., 190., 161])

river_Eposition = \
    np.array([ 41.,  42., 43., 44., 45., 46., 
               47.,  50., 51., 52., 53., 
               46., 46., 46., 45., 45., 39., 37., 
               40., 41., 41., 41., 42., 43.,   
               45.,  46., 47., 48., 49., 50., 51., 52.,  
               56., 57., 58., 59., 60., 61., 
               62., 63.,  64., 66., 67., 
               175., 175., 175., 
               175., 175., 175., 175., 175.,
               190., 187., 188., 150, 182., 180., 185])

river_direction = \
    np.array([ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  
               1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
               0., 0., 0., 0., 0., 0., 0., 0., 
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  
               1., 1., 1., 
               1., 1., 1., 1., 1.,
               1., 1., 1., 1., 1., 1., 1.])

# multiplication factor for each point source.  Multiply by Mississippi River
# discharge to get the right point source flux.
river_factor = \
    np.array([ -0.11666688, -0.11666688,-0.11666688,
            -0.01624998, -0.01624998, -0.01624998, -0.01624998, -0.01624998,
            -0.01624998,-0.01624998,-0.01624998,-0.01624998,-0.01624998,
            -0.01624998,-0.01624998,-0.01624998,-0.01624998,-0.01624998,
            -0.01624998,-0.01624998,-0.01624998,-0.01624998,-0.01624998,
            -0.01624998, 0.01624998, 0.01624998, 0.01624998, 0.01624998,
             0.01624998, 0.01624998, 0.01624998, 0.01624998, 0.01624998,
             0.01624998, 0.01624998, 0.01624998, 0.01624998, 0.01624998,
             0.01624998, 0.01624998, 0.01624998, 0.01624998, 0.01624998,   
            -0.100, -0.100, -0.100, 
            -0.140, -0.140, -0.140, -0.140, -0.140,
            -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000])

river_flag = 5

# Set river temp equal to local air temp
nc = netCDF4.Dataset('ss_air_temp.cdf')
air_temp = nc.variables['sat']
air_temp.set_auto_maskandscale(True)
air_temp = air_temp[:,118,270]
nc.close()

sat = []
for yr in range(1986, 2011):
    for mo in range(1, 13):
        jd = utime.date2num(datetime(yr, mo, 15))
        sat.append((jd, air_temp[mo-1]))

sat = np.asarray(sat, dtype=[('jd', 'd'), ('air_temp', 'd')])
river_temp = np.interp(river_jd, sat['jd'], sat['air_temp'])
river_temp = river_temp[:,np.newaxis, np.newaxis] * np.ones((30, len(river)))

river_salt = 0.0

river_dye_01 = octant.ocean.o2_sat(river_temp, river_salt)

river_dye_02 = \
    np.array([ 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1.,
            0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0.])
river_dye_02 = river_dye_02[:,np.newaxis,np.newaxis]*np.ones((len(river_jd),30))
river_dye_02 = river_dye_02.swapaxes(0,1)
river_dye_02 = river_dye_02.swapaxes(1,2)

river_dye_03 = \
    np.array([ 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            1., 1., 1.,
            1., 1., 1., 1., 1.,
            0., 0., 0., 0., 0., 0., 0.])
river_dye_03 = river_dye_03[:,np.newaxis,np.newaxis]*np.ones((len(river_jd),30))
river_dye_03 = river_dye_03.swapaxes(0,1)
river_dye_03 = river_dye_03.swapaxes(1,2)

river_dye_04 = \
    np.array([ 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0.,
            0., 0., 0., 0., 0.,
            0., 0., 0., 1., 0., 0., 0.])
river_dye_04 = river_dye_04[:,np.newaxis,np.newaxis]*np.ones((len(river_jd),30))
river_dye_04 = river_dye_04.swapaxes(0,1)
river_dye_04 = river_dye_04.swapaxes(1,2)

river_Miss = np.zeros_like(river_temp)
river_Miss[:,:,river_Xposition > 400] = 1.0
river_Atch = np.zeros_like(river_temp)
river_Atch[:,:,river_Xposition < 400] = 1.0

vshape_default = \
    np.array([ 0.0022 ,  0.0043,  0.0065,  0.0086,  0.0108,
               0.0129,   0.0151,  0.0172,  0.0194,  0.0215,
               0.0237,   0.0258,  0.0280,  0.0301,  0.0323,
               0.0344,   0.0366,  0.0387,  0.0409,  0.0430,
               0.0452,   0.0473,  0.0495,  0.0516,  0.0538,
               0.0559,   0.0581,  0.0602,  0.0624,  0.0645,
             ])[:,np.newaxis]

vshape_swpass = \
    np.array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.       , 0.         ,
            0.        ,  0.        ,  0.        ,  0.0128    , 0.0256     ,
            0.0385    , 0.0513     , 0.0641     , 0.0769    ,  0.0897     ,
            0.1026    , 0.1154     , 0.1282     , 0.1410    , 0.1538,
            ])[:,np.newaxis]


river_Vshape = vshape_default * np.ones( (1, len(river)), 'd')
river_Vshape[:, 0:3] = vshape_swpass

river_transport = np.ones((len(river_jd),len(river)))
river_transport[:,0:43] = river_data['q'][:, np.newaxis] * river_factor[np.newaxis, 0:43]
river_transport[:,43:51] = 0.5*river_data['q'][:, np.newaxis] * river_factor[np.newaxis, 43:51]
river_transport[:,51:52] = river_calcasieu['q'][:,np.newaxis] * river_factor[np.newaxis, 51:52]
river_transport[:,52:53] = river_sabine['q'][:,np.newaxis] * river_factor[np.newaxis, 52:53]
river_transport[:,53:54] = river_trinity['q'][:,np.newaxis] * river_factor[np.newaxis, 53:54]
river_transport[:,54:55] = river_brazos['q'][:,np.newaxis] * river_factor[np.newaxis, 54:55]
river_transport[:,55:56] = river_lavaca['q'][:,np.newaxis] * river_factor[np.newaxis, 55:56]
river_transport[:,56:57] = river_sanantonio['q'][:,np.newaxis] * river_factor[np.newaxis, 56:57]
river_transport[:,57:58] = river_nueces['q'][:,np.newaxis] * river_factor[np.newaxis, 57:58]
#river_transport[:,51:52] = river_brazos['q'][:,np.newaxis] * river_factor[np.newaxis, 51:52]
#river_transport[:,52:53] = river_nueces['q'][:,np.newaxis] * river_factor[np.newaxis, 52:53]
#river_transport[:,53:54] = river_sanantonio['q'][:,np.newaxis] * river_factor[np.newaxis, 53:54]
#river_transport[:,54:55] = river_lavaca['q'][:,np.newaxis] * river_factor[np.newaxis, 54:55]
#river_transport[:,55:56] = river_trinity['q'][:,np.newaxis] * river_factor[np.newaxis, 55:56]
#river_transport[:,56:57] = river_sabine['q'][:,np.newaxis] * river_factor[np.newaxis, 56:57]
#river_transport[:,57:58] = river_calcasieu['q'][:,np.newaxis] * river_factor[np.newaxis, 57:58]

###################################
# Crate netCDF file

# dimensions:
#   xi_rho = 575 ;
#   eta_rho = 191 ;
#   s_rho = 30 ;
#   river = 10 ;
#   river_time = 3652 ;
# variables:
#   double river(river) ;
#   double river_Xposition(river) ;
#   double river_Eposition(river) ;
#   double river_direction(river) ;
#   double river_Vshape(s_rho, river) ;
#   double river_time(river_time) ;
#   double river_transport(river_time, river) ;
#   double river_temp(river_time, s_rho, river) ;
#   double river_salt(river_time, s_rho, river) ;
#   double river_flag(river) ;
#   double river_O2(river_time, s_rho, river) ;
#   double river_Miss(river_time, s_rho, river) ;
#   double river_Atch(river_time, s_rho, river) ;
# 
# // global attributes:
#       :title = "TGLO river forcing file" ;
#       :type = "FORCING file" ;
#       :author = "Robert Hetland" ;
#       :date = "20-Nov-2003 18:35:37" ;
# }

nc = netCDF4.Dataset('TXLA_river_4dyes_2010.nc', 'w', format='NETCDF3_CLASSIC')

nc.createDimension('x_rho', 671)
nc.createDimension('eta_rho', 191)
nc.createDimension('s_rho', 30)
nc.createDimension('river', len(river))
nc.createDimension('river_time', len(river_jd))

def write_nc_var(name, dimensions, var, units=None):
    nc.createVariable(name, 'f8', dimensions)
    if units is not None:
        nc.variables[name].units = units
    nc.variables[name][:] = var

write_nc_var('river', ('river', ), river, 'River index')
write_nc_var('river_Xposition', ('river', ), river_Xposition, 'River xi index')
write_nc_var('river_Eposition', ('river', ), river_Eposition, 'River eta index')
write_nc_var('river_direction', ('river', ), river_direction, 'River direction')
write_nc_var('river_flag', ('river', ), river_flag, 'River flag')

write_nc_var('river_Vshape', ('s_rho', 'river'), river_Vshape, 'River vertical shape')
write_nc_var('river_time', ('river_time', ), river_jd, time_units)
write_nc_var('river_transport', ('river_time', 'river'), river_transport, 
             'River transport [m3/s]')
write_nc_var('river_temp', ('river_time', 's_rho', 'river'), 
              river_temp,
             'River temperature [deg C]')
write_nc_var('river_salt', ('river_time', 's_rho', 'river'), river_salt,
             'River salinity [psu]')

# Tracers
write_nc_var('river_dye_01', ('river_time', 's_rho', 'river'), river_dye_01,
            'River oxygen [mmol/l]')
write_nc_var('river_dye_02', ('river_time', 's_rho', 'river'), river_dye_02,
            'Mississippi river tag [none]')
write_nc_var('river_dye_03', ('river_time', 's_rho', 'river'), river_dye_03,
            'Atchafalaya river tag [none]')
write_nc_var('river_dye_04', ('river_time', 's_rho', 'river'), river_dye_04,
            'Brazos river tag [none]')

#write_nc_var('river_Miss', ('river_time', 's_rho', 'river'), river_Miss,
#            'Mississippi river tag [none]')
#write_nc_var('river_Atch', ('river_time', 's_rho', 'river'), river_Atch,
#            'Atchafalaya river tag [none]')

nc.close()
