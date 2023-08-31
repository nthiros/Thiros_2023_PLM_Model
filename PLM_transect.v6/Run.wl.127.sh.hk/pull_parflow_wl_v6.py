from parflowio.pyParflowio import PFData

import numpy as np
import pandas as pd
import os
import pdb
import glob

class pull_pfb_pressure():
    def __init__(self, pf_wells_df, pf_dz_cumsum_c, pf_zmax, dir_loc, run_name):
        self.pf_wells_df = pf_wells_df     # Where in ParFlow grid the wells are -- run ER_PLM_ParFlow/utils/well_info_v4.dummywells.py first
        self.dz_cumsum   = pf_dz_cumsum_c  # ParFlow variable z depths -- cell centered value from bottom to top of domain
        self.zmax        = pf_zmax         # ParFlow total depth
        self.dir_loc     = dir_loc         # Location of parflow .pfb files
        self.run_name    = run_name        

        self.wt  = None
        self.wt_ = None
    
    def find_pfb(self):
        '''Gathers all press.pfb files in directory with path dir_loc into a list.
           Returns a lists'''
        #pdb.set_trace()
        ff = []
        for file in os.listdir(self.dir_loc):  
            if file.endswith(".pfb"):
                #if 'press' in file.split('.'):
                if all(x in file.split('.') for x in ['press',self.run_name]):
                    ff.append(os.path.join(self.dir_loc, file))
        ff.sort()
        return ff 

    def read_pfb(self, pfb_fname):
        '''parflowio magic to read the pfb files into numpy array.
           Returns 2D array with shape NZ,NX'''
        pfdata = PFData(pfb_fname)
        pfdata.loadHeader()
        pfdata.loadData()
        dd = pfdata.copyDataArray()
        return dd[:,0,:] # because our y-dim is zero, ie. a 2D model
        
    def pf_pressure(self, pfb_fname):
        '''Pull simulated water table below land surface at well locations.
           This is referenced to the cell-centered grid.
           Returns array with distance below land surface in m'''
        #pdb.set_trace()
        # Mapping indices numbers to well names to be sure indexing done correctly
        ind_map = {}
        wells = self.pf_wells_df.index.to_list()
        for i in range(len(wells)):
            ind_map[wells[i]] = i
        
        # ParFlow cell indices at X and sampling depth Z
        wellinds = self.pf_wells_df[['Cell_Z', 'Cell_X']].to_numpy().astype(int)
        
        # Pressure at all grid points 
        dd = self.read_pfb(pfb_fname)
        
        # Calculate wl elev using first cell with postive pressure head -- moving from land surface down
        # Pressure head profile along column located at well location X
        dd_col  = dd[:, wellinds[:,1]]
        # Z-cell that is closest to pressure head >= 0 (water table)
        dd_zero = (dd_col>0.0).argmin(axis=0) - 1
        # Pressure head at this cell
        pr_wt = dd_col[dd_zero, np.arange(len(dd_zero))]
        # Depth below land surface to first fully saturated cell
        h_bls= self.dz_cumsum[dd_zero]
        
        # Now pull pressure at sampling location
        # Because PLM1 and PLM6 are peizometers, want screen depth for reference height
        # cell indices for sampling loc
        #z_samps = [abs(z_info['Depth_bls'] - pf_wells['screen_top_m'][i]).idxmin() for i in list(ind_map.keys())] # at top of screen elevation
        z_samps   = self.pf_wells_df['Cell_Z'].astype(int) # at sampling depth
        # Pressure head at sampling depth z
        pr_pz = dd_col[z_samps.to_numpy().T, list(ind_map.values())] 
        # Depth below land surface to the piezometer cell
        h_bls_pz  = self.dz_cumsum[z_samps]
        
        # combine for PLM1, PLM7, PLM6
        h  = h_bls_pz.copy() # depth to cell
        pr = pr_pz.copy() # pressure at cell
        
        # Account for pressure head at this cell
        wt_bls = h - pr
        return wt_bls
        
    
    def get_water_table_ts(self):
        '''Loop through all times and pull water table depth (below top of model).
           Also calculates the water table elevation (meters above sea-level).
           Returns arrays of water table below land surface and water table elevation at well X locations'''
        wt_df  = pd.DataFrame(index=self.pf_wells_df.index)
        wt_df_ = pd.DataFrame(index=self.pf_wells_df.index)
        
        pfb_list = self.find_pfb()
        for i in pfb_list:
            wt_bls = self.pf_pressure(i)
            #wt_ele = self.pf_wells_df['land_surf_dem_m'].to_numpy() - wt_bls
            t  = float(i.split('.')[-2])
            
            wt_df.loc[:, t] = wt_bls
            #wt_df_.loc[:, t] = wt_ele
            
        self.wt  = wt_df.T
        #self.wt_ = wt_df_.T
        #return wt_df.T, wt_df_.T
        return wt_df.T

    def mangle_dates(self, startdate):
        '''Put timeseries into acutal dates.
           Assumes 24 hour timesteps right now'''
        startd = pd.to_datetime(startdate)
        endd   = startd + pd.Timedelta(len(self.wt)-1, 'D')
        date_list = pd.date_range(start=startd, end=endd, freq='D')
        
        self.wt.index  = date_list
        #self.wt_.index = date_list
        
        return self.wt #, self.wt_
        

# Read in well info dataframe
# Produced by ER_PLM_ParFlow/utils/well_info_v4.dummywells.py
#pf_wells = pd.read_csv('../ER_PLM_ParFlow/utils/wells_2_pf_v4.dummy.csv', index_col='well')
#z_info   = pd.read_csv('../ER_PLM_ParFlow/utils/plm_grid_info_v4.dummy.csv') # these are cell centered values
#pf_wells = pd.read_csv('../utils/wells_2_pf_v4.dummy.csv', index_col='well')
#z_info   = pd.read_csv('../utils/wells_2_pf_v4.dummy.csv') # these are cell centered values


# ParFlow Grid Info
DZ = 10.0
# Variable Z scaling
# Note - index 0 is the bottom on the domain
#dz_scale = np.array([1.00, 1.00, 1.00, 1.00, 1.00,       # 52.0 - 102.0
#                     0.80, 0.80,                         # 36.0 - 52.0
#                     0.60, 0.60,                         # 24.0 - 36.0
#                     0.40, 0.40,                         # 16.0 - 24.0
#                     0.20, 0.20, 0.20,                   # 10.0 - 16.0  -- 2m layers down to 16 m
#                     0.10, 0.10,                         # 8.0  - 10.0
#                     0.10, 0.05, 0.05,                   # 6.0  - 8.0   -- 0.5m res possible down to 7.0 m
#                     0.05, 0.05, 0.05, 0.05,             # 4.0  - 6.0
#                     0.05, 0.05, 0.05, 0.05,             # 2.0  - 4.0
#                     0.05, 0.05, 0.05, 0.025, 0.025])    # 0.0  - 2.0  
# Read dz_scale in directly from file
with open(glob.glob('*.mat*')[0], 'r') as f:
    l = f.readlines()
cell_id  = []
dz = []
for ll in l:
    try:
        ll.split()[1].split('.')[1]
    except IndexError:
        pass
    else:
        if ll.split()[1].split('.')[0] == 'Cell':
            cell_id.append(ll.split()[1].split('.')[1])
            dz.append(float(ll.split()[2]))

#Layer thickness, counting from bottom up
dz_scaled = np.array(dz) * DZ
# Cumulative thickness, counting from bottom up
dz_cumsum = dz_scaled.cumsum()
# Cell centered z values, counting from bottom of domain to the top of the domain 
dbs_c = dz_cumsum.max() - (dz_cumsum) + dz_scaled/2


# Manually building the well position array for now...
wells      = ['PLM1','PLM7','PLM6','X404','X494','X508']
smp_depth  = [7.0, 20.0, 9.0, 0.025, 0.025, 0.025] # where to calculate pressure head from
xpos       = [404, 424, 494, 404, 494, 508] # index for parflow cell in x direction
# find parflow cell closest to sample depth
pf_wells = pd.DataFrame()
for i in range(len(wells)):
    pfz = abs(dbs_c - smp_depth[i]).argmin()
    pf_wells.loc[wells[i], 'Cell_Z'] = pfz
pf_wells['Cell_X'] = xpos







# Pull Simulated water levels
print (' ')
for ff in ['wy_spinup','wy_2000_2016','wy_2017_2021']:
    if ff in os.listdir():
        print ('working on {} water levels'.format(ff))
        pf_wt = pull_pfb_pressure(pf_wells, dbs_c, dz_cumsum.max(), ff, ff)
        #pf_wt_bls, pf_wt_elev = pf_wt.get_water_table_ts()
        pf_wt_bls = pf_wt.get_water_table_ts()
        pf_wt_bls.index.name = 'Timestep'
        #pf_wt_elev.index.name = 'Timestep'
        #pf_wt_bls, pf_wt_elev = pf_wt.mangle_dates('2014-10-01',)

        # save it
        if not os.path.exists('parflow_out'):
            os.makedirs('parflow_out')
        pf_wt_bls.to_csv('./parflow_out/{}_wt_bls.csv'.format(ff))
        #pf_wt_elev.to_csv('./parflow_out/{}_wt_elev.csv'.format(ff))
    else:
        pass
print (' ')










