#!/usr/bin/env python
# coding: utf-8

# Update 08/14/2021
# All MET forcing needs to be hourly values
# For upscaling, take average of all MET fields

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator
import matplotlib.dates as mdates
plt.rcParams['font.size']=12



# Desired fields for CLM forcing files
units = {'rad_s':'Short Wave Radiation [W/m^2]',
         'rad_l':'Long Wave Radiation [W/m^2]',
         'prcp':'Precipitation [mm/s]',
         'temp':'Air Temperature [K]',
         'wnd_u':'East-to-West Wind [m/s]',
         'wnd_v':'South-to-North Wind [m/s]',
         'press':'Atmospheric Pressure [pa]',
         'vap':'Water-vapor Specific Humidity [kg/kg]'}



# Function to write MET forcing for CLM
def write_met(df, fname):
    # Write a MET file for CLM from a dataframe
    fmt = '%.3f', '%.2f', '%.2E', '%.2f', '%.2f', '%.2f', '%.1f', '%.5f'
    np.savetxt(fname, df.to_numpy(), fmt=fmt, delimiter='\t')



#-------------------------------------------------------------
# MET Hourly Forcing Organized in Water Years
# Note! this was just for timestep testing purposes
#-------------------------------------------------------------
# combine into single array
wy = [15]
met_ = [np.loadtxt('NLDAS_met_hourly/wy{}_forcing.txt'.format(i)) for i in wy]
met = np.vstack(met_)


# Move to pandas
tstart = pd.to_datetime('20{}-10-01 00'.format(min(wy)-1), format='%Y-%m-%d %H')
tend = pd.to_datetime('20{}-09-30 23'.format(max(wy)), format='%Y-%m-%d %H')
hours = pd.Series(pd.date_range(tstart, tend, freq='H'))


# Dataframe based on hourly NLDAS-2
met_hr = pd.DataFrame(data=met, index=hours, columns=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])
print (met_hr.shape)


# Use built-in pandas utilities to resmaple timeseries
met_3h = met_hr.resample('3H').agg(dict(rad_s='mean',rad_l='mean',prcp='mean',temp='mean',wnd_u='mean',wnd_v='mean',press='mean',vap='mean'))
print (met_3h.shape)


# 3hr
write_met(met_hr,  'wy15_forcing.txt.1h')
write_met(met_3h, 'wy15_forcing.txt.3h')







#---------------------------------------------------
# Hourly NLDAS forcing 
# Converted to 3 hr timesteps for ParFlow
#---------------------------------------------------
# Combine NLDAS forcing arrays -- these are hourly right now
# wy 1980 to 2021
wy   = np.arange(0,22)
met1 = [np.loadtxt('./NLDAS_met_hourly/wy{:02d}_forcing.txt'.format(i)) for i in wy]
met1 = np.vstack(met1)
met2 = [np.loadtxt('./NLDAS_met_hourly/wy{}_forcing.txt'.format(i)) for i in ['80_84','85_89','90_94','95_99']]
met2 = np.vstack(met2)
met  = np.concatenate((met2,met1))


# Move to pandas
#tstart = pd.to_datetime('1999-10-01 00', format='%Y-%m-%d %H')
tstart = pd.to_datetime('1979-10-01 00', format='%Y-%m-%d %H')
tend = pd.to_datetime('2021-08-30 12', format='%Y-%m-%d %H') # Water year 21 is not over yet
hours = pd.Series(pd.date_range(tstart, tend, freq='H'))


# Dataframe based on hourly NLDAS-2
met_hr  = pd.DataFrame(data=met, index=hours, columns=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])
met_3hr = met_hr.resample('3H').mean()
# Drop Leap years again
met_3hr = met_3hr[~((met_3hr.index.month == 2) & (met_3hr.index.day == 29))]


#----
# Get averages for spinup
#---
# This is 1 year of forcing data at 3hr timesteps
# has the average across all NLDAS
met_3hr_avg_ = met_3hr.groupby([met_3hr.index.month, met_3hr.index.day,met_3hr.index.hour]).mean()
met_3hr_avg = met_3hr_avg_.reindex([10,11,12,1,2,3,4,5,6,7,8,9], level=0)
met_3hr_avg_10yr = pd.concat([met_3hr_avg]*10)
met_3hr_avg_50yr = pd.concat([met_3hr_avg]*50)
write_met(met_3hr_avg_10yr, 'met.avg.10yr.3hr.txt')
write_met(met_3hr_avg_50yr, 'met.avg.50yr.3hr.txt')



#------
# Clip to smaller time regions -- clips from hour 00 on first day to hour 23 on last day?
#------
# wy1980 to wy2016
#met_3h_8016 = met_3hr['1979-10-01' : '2016-09-30'] 
#write_met(met_3h_8016, 'met.1980-2016.3hr.txt')


# wy2000 to wy2012
met_3h_0016 = met_3hr['1999-10-01' : '2012-09-30'] 
write_met(met_3h_0016, 'met.2000-2012.3hr.txt')


# wy2012 to wy2021
met_3hr_1721 = met_3hr['2012-10-01' :]  
write_met(met_3hr_1721, 'met.2013-2021.3hr.txt')





