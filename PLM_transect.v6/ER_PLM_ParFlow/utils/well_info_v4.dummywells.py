# script to generate some info about well locations.
# doing it here because tcl is hard

# Updated to 102 m deep model

import numpy as np
import pandas as pd
#from parflowio.pyParflowio import PFData


# Import DEM info -- Lat, Long, Elvations
#cx = pd.read_excel('./wksht_slopes_v2.xlsx', sheet_name='Sheet1')
#cx = cx[['X','lat','long','elevation']]

cx = pd.read_csv('DEM_transect_regrid10_extended.txt', delim_whitespace=True, index_col=0, names=['X','lat','long','elev'])
cx.index = cx.index.astype(int)

#-----------------------------------------
# Obsevation Well Coordinates
# lat, lon, easting, northing
#-----------------------------------------
# read in well locations
well_mint = pd.read_excel('well_lat_lon.xlsx', 'Sheet1', index_col='well')
well = well_mint.copy().loc[['PLM1','PLM7','PLM6'],:]

# find closest cx position to wells
# this only makes sense for wells on pumpouse transect
for i in range(len(well)):
    dist = np.sqrt(((cx[['lat','long']].to_numpy() - well[['utm_east','utm_north']].to_numpy()[i])**2).sum(axis=1))
    w = well.index[i]
    well.loc[w,'X']     = cx.loc[dist.argmin(), 'X']
    well.loc[w,'Cell_X']  = cx.index[dist.argmin()]
    well.loc[w,'land_surf_cx'] =  cx.loc[dist.argmin(), 'elev']

# Negative numbers mean the simulated wl will be biased too high
dem_cx_elev_diff = well['land_surf_dem_m'] - well['land_surf_cx']


well.drop(labels=['land_surf_lit_m'], axis=1, inplace=True)




#----------------------------
# ParFlow Grid
#----------------------------
Xmin = 0.0
Ymin = 0.0
Zmin = 0.0

NX = 559
NY = 1
NZ = 32

DX = 1.5125
DY = 1.5125
DZ = 10.0

#
# Variable Z scaling
# Note - index 0 is the bottom on the domain
dz_scale = np.array([1.00, 1.00, 1.00, 1.00, 1.00,       # 52.0 - 102.0
                     0.80, 0.80,                         # 36.0 - 52.0
                     0.60, 0.60,                         # 24.0 - 36.0
                     0.40, 0.40,                         # 16.0 - 24.0
                     0.20, 0.20, 0.20,                   # 10.0 - 16.0
                     0.10, 0.10, 0.10,                   # 7.0  - 10.0
                     0.05, 0.05, 0.05, 0.05, 0.05,       # 4.5  - 7.0  --> 0.0 to 7.0 (N=15 in 0.5m steps)
                     0.05, 0.05, 0.05, 0.05,             # 2.5  - 4.5
                     0.05, 0.05, 0.05, 0.05,             # 0.5  - 2.5  
                     0.025,0.025]) 


# Layer thickness, counting from bottom up
dz_scaled = np.array(dz_scale) * DZ
# Cumulative thickness, counting from bottom up
dz_cumsum = dz_scaled.cumsum()

# Depth below land surface, cumulative thickness
# Counts top down, because makes more sense to me
dbs = np.flip(dz_scaled).cumsum() 

# cell centered z values, counting from bottom of domain to the top of the domain 
dbs_c = dz_cumsum.max() - (dz_cumsum) + dz_scaled/2



#-------------------
# Update 11/3/21
# Add dummy wells upslope for comparison purposes
#-------------------
# Xcell index number
#Xpos = [404, 494, 528]
Xpos =[404, 494, 508]

#top_screen = 0.25
#bot_screen = 2.75
#smp_depth  = 1.75

top_screen = 0.0
bot_screen = 2.0
#bot_screen = 1.5
smp_depth  = 0.75


for i in range(len(Xpos)):
    well.loc['X{}'.format(Xpos[i]), 'Cell_X']       = Xpos[i]
    well.loc['X{}'.format(Xpos[i]), 'X']            = Xpos[i]*1.5125
    
    well.loc['X{}'.format(Xpos[i]), 'Cell_Z']       = abs(dbs_c-smp_depth).argmin()
    well.loc['X{}'.format(Xpos[i]), 'smp_depth_m']  = smp_depth
    well.loc['X{}'.format(Xpos[i]), 'screen_top_m'] = top_screen
    well.loc['X{}'.format(Xpos[i]), 'screen_bot_m'] = bot_screen
    
    well.loc['X{}'.format(Xpos[i]), 'land_surf_cx'] = cx.loc[Xpos[i], 'elev']
    well.loc['X{}'.format(Xpos[i]), 'land_surf_dem_m'] = cx.loc[Xpos[i], 'elev']





#----------------------------------------------
# Put well locations into ParFlow Model space
#----------------------------------------------
# Find what z layer each well is in
# using the sample depth 
for i in range(len(well)):
    wn  = well.index[i]
    wz  = well['smp_depth_m'][i]
    pfz = abs(dbs_c - wz).argmin()
    well.loc[wn, 'Cell_Z'] = pfz



#----
# ParFlow grid cell number of wells, makes Paraview plotting easier
well_pf_ind = well['Cell_X'] + NX*well['Cell_Z']
well['Cell_ind'] = well_pf_ind 


#---
# cell centered domain info
cc = np.column_stack((np.arange(len(dbs_c)), dbs_c, np.ones(len(dbs_c))*dbs.max()))


#-----
# Writes a file for later use
well.to_csv('wells_2_pf_v4.dummy.csv', index=True)

np.savetxt('plm_grid_info_v4.dummy.csv', cc, fmt='%.4f', delimiter=',', header='Z_layer_num,Depth_bls,Total_Depth',comments='')




"""
#---------------
# Now for shorter hillslope

# Trim it
# make sure matches 'transect_slope_setup.py'
cxt = cx.copy().loc[165:, :]
X_t = np.arange(0, 1.5125*len(cxt), 1.5125)
cxt.loc[:,'X'] = X_t
cxt.reset_index(drop=True, inplace=True)


#-----------------------------------------
# Obsevation Well Coordinates
# lat, lon, easting, northing
#-----------------------------------------
# read in well locations
#well = pd.read_excel('well_lat_lon.xlsx', 'Sheet1', index_col='well')well = well_mint.copy().loc[['PLM1','PLM7','PLM6'],:]

# find closest cx position to wells
# this only makes sense for wells on pumpouse transect
for i in range(len(well)):
    dist = np.sqrt(((cx[['lat','long']].to_numpy() - well[['utm_east','utm_north']].to_numpy()[i])**2).sum(axis=1))
    w = well.index[i]
    well.loc[w,'X']     = cx.loc[dist.argmin(), 'X']
    well.loc[w,'Xind']  = cx.index[dist.argmin()]
    well.loc[w,'land_surf_cx'] =  cx.loc[dist.argmin(), 'elev']






#----------------------------
# ParFlow Grid
#----------------------------
Xmin = 0.0
Ymin = 0.0
Zmin = 0.0

NX = len(cxt)
NY = 1
NZ = 26

DX = 1.5125
DY = 1.5125
DZ = 10.0

#
# Variable Z scaling
# Note - index 0 is the bottom on the domain
dz_scale = np.array([0.4, 0.4, 0.4, 0.4, 0.4,       
                     0.2, 0.2, 0.2, 0.2, 0.2,      
                     0.1, 0.1, 0.1, 0.1, 0.1,       
                     0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
                     0.025,0.025])    


# Layer thickness, counting from bottom up
dz_scaled = np.array(dz_scale) * DZ
# Cumulative thickness, counting from bottom up
dz_cumsum = dz_scaled.cumsum()

# Depth below land surface, cumulative thickness
# Counts top down, because makes more sense to me
dbs = np.flip(dz_scaled).cumsum() 


#----------------------------------------------
# Put well locations into ParFlow Model space
#----------------------------------------------
# Cell-centered depths, starts at land surface
dbs_c  = dbs - np.flip(dz_scaled/2)
dbs_c_ = np.flip(dbs_c) # make index 0 bottom of domain

# Find what z layer each well is in
for i in range(len(well)):
    wn  = well.index[i]
    wz  = well['smp_depth_m'][i]
    pfz = abs(dbs_c_ - wz).argmin()
    well.loc[wn, 'pf_Z_num'] = pfz



#----
# ParFlow grid cell number of wells, makes Paraview plotting easier
well_pf_ind = well['X'] + NX*well['pf_Z_num']
well['cell_ind'] = well_pf_ind 



#-----
# Write a file for later use
well.to_csv('plm_grid_info_t.csv', index=True)
""" 
    
    
