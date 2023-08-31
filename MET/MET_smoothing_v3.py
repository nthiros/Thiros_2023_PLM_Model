# This script generates the ParFlow-CLM forcing conditions from the NLDAS-2 dataset
# Hourly transience from 1980 through 2021
# Raw NLDAS-2 data in NLDAS_met_hourly directory

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



#
# Hourly NLDAS forcing 
#
# Combine NLDAS forcing arrays -- these are hourly right now
wy   = np.arange(0,22)
met1 = [np.loadtxt('./NLDAS_met_hourly/wy{:02d}_forcing.txt'.format(i)) for i in wy]
met1 = np.vstack(met1)
met2 = [np.loadtxt('./NLDAS_met_hourly/wy{}_forcing.txt'.format(i)) for i in ['80_84','85_89','90_94','95_99']]
met2 = np.vstack(met2)
met  = np.concatenate((met2,met1))


# Move to pandas
tstart = pd.to_datetime('1979-10-01 00', format='%Y-%m-%d %H')
tend   = pd.to_datetime('2021-09-30 23', format='%Y-%m-%d %H') 
hours  = pd.Series(pd.date_range(tstart, tend, freq='H'))

# Dataframe based on hourly NLDAS-2
met_hr  = pd.DataFrame(data=met, index=hours, columns=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])


# Get averages for spinup
# This is 1 year of forcing data at 1hr timesteps
# has the average across all NLDAS - e.g. average precip on Jan.1 over the 30 years of NLDAS
met_1hr_avg_ = met_hr.groupby([met_hr.index.month, met_hr.index.day, met_hr.index.hour]).mean()
met_1hr_avg = met_1hr_avg_.reindex([10,11,12,1,2,3,4,5,6,7,8,9], level=0)
met_1hr_avg_10yr = pd.concat([met_1hr_avg]*10) # 10 years repeated
met_1hr_avg_50yr = pd.concat([met_1hr_avg]*50) # 50 years repeated
write_met(met_1hr_avg_10yr, 'met.avg.10yr.1hr.txt')
write_met(met_1hr_avg_50yr, 'met.avg.50yr.1hr.txt')



# Hourly Forcing Files 
met_1h_8021 = met_hr['1979-10-01':] 
write_met(met_1h_8021, 'met.1980-2021.1hr.txt')

met_1h_8016 = met_hr['1979-10-01' : '2016-09-30'] 
write_met(met_1h_8016, 'met.1980-2016.1hr.txt')


# wy2000 to wy2016
met_1h_0016 = met_hr['1999-10-01' : '2016-09-30'] 
write_met(met_1h_0016, 'met.2000-2016.1hr.txt')


# wy2016 to present
met_1h_1721 = met_hr['2016-10-01' :]  
write_met(met_1h_1721, 'met.2017-2021.1hr.txt')


# wy2012 to present -- for water level MC simulations
met_1h_1221 = met_hr['2012-10-01' :]  
write_met(met_1h_1221, 'met.2012-2021.1hr.txt')



