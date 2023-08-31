# Parse a directory full of nldas grb files into single timeseries.
# Need to first run ./parse_nldas_grb.sc
# This assumes that the grb have only single point in them, not a full raster of points.

import pandas as pd
import numpy as np
import os 
import subprocess


#-----------------------------------------
# Read in the metadata file
#-----------------------------------------
with open('metadata_full.txt', 'r') as f:
    ll = f.readlines()

unit_code_l = [] # Full description of the fields
unit_code   = [] # Short codes for the fields

for i in range(len(ll)):
    for j in range(len(ll[i].split())):
        if (ll[i].split()[j][:3] == 'N/S'):
            uglyline = ll[i].split()[j:]
            unit_code_l.append(' '.join([uglyline[0][5:]] + uglyline[1:]))
    unit_code.append(ll[i].split()[0].split(':')[3])
    
    
    
#-----------------------------------------
# Read in the metadata file
#-----------------------------------------
# Generate a list of all the nldas txt files in currently directory
def find_grb(dir_loc):
    '''Gathers all *.grb.txt files in directory with path dir_loc into a list.
       Returns a list of lists'''
    ff = []
    for file in os.listdir(dir_loc):
        if file.endswith(".grb.txt"):
            ff.append(os.path.join(dir_loc, file))
    return ff
grb_files = find_grb('./')

dates_ = [] # hold the dates for the data
data_  = [] # the actual data

# Loop through and pull the data
for i in grb_files:
    # First pull the dates from file name
    d = i.split('.')[-7][1:]
    h = i.split('.')[-6]
    dates_.append(d+'-'+h)
    # Pull data
    data_.append(np.loadtxt(i))
    
    
#-----------------------------------------
# Make a Dataframe
#-----------------------------------------
df = pd.DataFrame(data=np.array(data_), columns=unit_code, index=pd.to_datetime(dates_, format='%Y%m%d-%H%M'))
df.sort_index(inplace=True)

# Put into order that CLM wants
df_ = df[['DSWRF','DLWRF','APCP','TMP','UGRD','VGRD','PRES','SPFH']]

# Convert precip from kg/m2 to mm/s 
# NLDAS pixels are 1/8 degree x 1/8 degree
# Use https://www.engr.scu.edu/~emaurer/tools/calc_cell_area_cgi.pl
# 1° = 98 km at Crested Butte Latitude (vs 1° = 111 km at equator)
A = 98*98/(8**2) * 1000**2 #m^2
rho = 1000 #kg/m3
##sec_per_month = df_.index.daysinmonth * 60*60.
df_.loc[:,'APCP'] = df_.loc[:,'APCP']/rho*1000/3600


# Write a file
fmt = '%.3f', '%.2f', '%.2E', '%.2f', '%.2f', '%.2f', '%.1f', '%.5f'
np.savetxt('wy20_forcing.txt', df_.to_numpy(), fmt=fmt, delimiter='\t')
#df_.to_csv('wy20_forcing-.csv', index=True, header=True)
