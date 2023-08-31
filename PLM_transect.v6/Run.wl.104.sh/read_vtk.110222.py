# Update 04/19/2022
# Need to verify that particles that are pulled are in the fractured bedrock units
# Update 11/02/2022
# Pulls mean age dynamics in soil, saprolite, and bedrock lithologies
# Defines boreholes and pulls full RTDs across full model depth at select times


import pandas as pd
import numpy as np
import os
import pickle
import glob

from parflowio.pyParflowio import PFData
import pyvista as pv
import parflow.tools as pftools

import matplotlib.pyplot as plt

import pdb



#------------------------------------------------------
#
# External info needed to parse EcoSlim vtk
#
#------------------------------------------------------   
# Read in well info -- created by well_info_v2.py
well_df = pd.read_csv('../ER_PLM_ParFlow/utils/wells_2_pf_v4.dummy.csv', index_col=('well'))


# Update some of this on the fly  to ccount for deeper and shallower models
# ParFlow Grid Info
DZ = 10.0
NX = 559
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

NZ = len(dz)
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
for i in range(len(wells)):
    pfz = abs(dbs_c - smp_depth[i]).argmin()
    well_df.loc[wells[i], 'Cell_Z'] = pfz
well_df['Cell_X'] = xpos
well_df['Cell_ind'] = well_df['Cell_X'] + NX*well_df['Cell_Z']




#----------------------------------------------------
#
# Class to pull from EcoSlim_cgrid.vtk files
#
#----------------------------------------------------
def read_pfb(fname):
    '''Read in a .pfb file'''
    pfdata = PFData(fname)
    pfdata.loadHeader()
    pfdata.loadData()
    return pfdata.copyDataArray()


class ecoslim_grid_vtk():
    def __init__(self, wells_df):
        
        # List of vtk files to process
        self.vtk_files = None
        
        # Observation wells locations
        self.wells = wells_df   # output from well_info_v2.py
        
        # vtk info
        self.cell_xyz = None   # cell centers -- from: read_vtk_grid
        self.var = None        # variable names from the vtk file -- from: read_vtk_grid
        
        # Age df
        self.age_df = None     # final dataframe with ecoslim output at wells
               
    def find_cgrid_vtk(self, dir_loc):
        '''Gathers all *.vtk files in directory with path dir_loc into a list.
           Returns a list of lists'''
        ff = []
        for file in os.listdir(dir_loc):
            if file.endswith(".vtk"):
                if file.split('.')[0].split('_')[-1] == 'cgrid':
                    ff.append(os.path.join(dir_loc, file))
        ff.sort()
        self.vtk_files = ff
        return ff
        
    def read_vtk_grid(self, vtkfile):
        '''This generates center coordinates for all cells.
           Only needs to be run once.'''
        #pdb.set_trace()
        mesh = pv.read(vtkfile) 
        # Pull the center coordinate (XYZ) for each cell
        cent_avg = lambda coords: [np.mean(coords[i:i+2]) for i in [0,2,4]]
        cell_xyz = np.array([cent_avg(mesh.cell_bounds(n)) for n in range(mesh.n_cells)])
        self.cell_xyz = cell_xyz
        self.vars = list(mesh.cell_data)
        return cell_xyz
       
    def update_df(self, mesh, df, time, varname):
        '''Utility function to store vtk output at points and times into a dataframe'''
        #pdb.set_trace()
        inds = self.wells['Cell_ind'].to_numpy().astype(int)
        df.loc[:,time] = mesh[varname][inds]
        return df
    
    def update_df_layers(self, mesh, df_layers, time):
        '''Utility function to store vtk output at times into a dataframe.
           This pulls Mean Age for entire: soil, saprolite, and bedrock layers -- not at single well points.'''
        #db.set_trace()
        # Read in a porosity field to define layers
        directory = 'wy_2017_2021' # hard code for now
        header    = 'wy_2017_2021' # hard code for now
        fn_porosity = os.path.join(directory, '{}.out.{}.pfb'.format(header,'porosity'))
        por = np.flipud(read_pfb(fn_porosity)[:,0,:])
        
        soil_ind = np.where(por[:,0]==por.max())[0].max()+1
        bed_ind  = np.where(por[:,0]==por.min())[0].min()
        
        age_arr = np.flipud(mesh['Age'].reshape(NZ,559))
        
        df_layers.loc['soil_mu',time] = age_arr[:soil_ind, :-20].mean() # drop the toe of hillslope where BC looks to be factor? Floodplain well is at X528, so keep 559-l>528
        df_layers.loc['soil_sd',time] = age_arr[:soil_ind, :-20].std()
        df_layers.loc['sap_mu',time]  = age_arr[soil_ind:bed_ind, :].mean()
        df_layers.loc['sap_sd',time]  = age_arr[soil_ind:bed_ind, :].std()
        df_layers.loc['bed_mu',time]  = age_arr[bed_ind:, :].mean()
        df_layers.loc['bed_sd',time]  = age_arr[bed_ind:, :].std()
        return df_layers
    
    def update_df_borehole(self, mesh, xinds, df_borehole, time):
        '''Utility function to store vtk output at times into a dataframe.
           This pulls Mean Age along a vertical 'borehole' '''
        #pdb.set_trace()
        age_arr = np.flipud(mesh['Age'].reshape(NZ,559))
        
        for i in range(len(df_borehole)):
            df_borehole[i].loc[np.arange(NZ),time] = np.array(age_arr[:,xinds[i]])
        return df_borehole
    
    def read_vtk(self):
        '''Loop through all vtk files and pull mean age at wells'''
        #pdb.set_trace()
        #inds = self.wells['Cell_ind'].to_numpy().astype(int)
        
        # seperate dataframe for all variables
        gen_df = lambda : pd.DataFrame(index=self.wells.index)
        df0,df1,df2,df3,df4,df5,df6,df7,df8 = [gen_df() for i in range(len(self.vars))]
        self.age_df = [df0,df1,df2,df3,df4,df5,df6,df7,df8]
        
        # holds mean age in soil, saprolite, and bedrock
        self.age_layers_df = pd.DataFrame()
        # borehole mean ages
        self.age_boreX65  = pd.DataFrame(index=np.arange(NZ))
        self.age_boreX265 = pd.DataFrame(index=np.arange(NZ))
        self.age_boreX424 = pd.DataFrame(index=np.arange(NZ))
        self.age_boreX528 = pd.DataFrame(index=np.arange(NZ))
        self.age_boreholes = [self.age_boreX65, self.age_boreX265, self.age_boreX424, self.age_boreX528]
        xinds  = [65, 265, 424, 528]
        xinds_ = ['X65', 'X265', 'X424', 'X528']
        
        
        # populate the dataframes from each timestep
        for i in range(len(self.vtk_files)): # loop through each timestep
            f = self.vtk_files[i]
            mesh_i = pv.read(f)
            ii = int(f.split('.')[-2])
            for j in range(len(self.vars)): # loop through each variable
                self.update_df(mesh_i, self.age_df[j], ii, self.vars[j])
            # layer mean age
            self.update_df_layers(mesh_i, self.age_layers_df, ii)
            # borehole mean age
            self.update_df_borehole(mesh_i, xinds, self.age_boreholes, ii)
        #pdb.set_trace()
        out_dict = {} # Dictionary of Dataframes
        for j in range(len(self.vars)):
            out_dict[self.vars[j]] = self.age_df[j].T.sort_index()
        out_dict['Age_layers'] = self.age_layers_df.T
        for z,r in zip(xinds_, self.age_boreholes):
            out_dict[z] = r
        return out_dict
    

## Moved Below  
## Mean Age Data
## WY 2017-2021
#age           = ecoslim_grid_vtk(well_df)
#vtk_files_c   = age.find_cgrid_vtk('./ecoslim_2017_2021')
#cell_xyzm     = age.read_vtk_grid(vtk_files_c[0])
#age_dict      = age.read_vtk()






#----------------------------------------------------------
#
# The pnts vtk files now
#
#----------------------------------------------------------
class ecoslim_pnts_vtk():
    def __init__(self, wells_df, eco_grid):
        
        # List of vtk files to process
        self.vtk_files = None
        
        # ParFlow-Ecoslim Grid Info
        self.cell_xyz = eco_grid # output from ecoslim_grid_vtk.read_vtk_grid()
        
        # Observation well locations
        self.wells_df = wells_df  # output from well_info_v2.py

        # Age Distribution
        self.rtd = {}
        
        # PID lists
        self.pid_list = []
        

    def find_pnts_vtk(self, dir_loc):
        '''Gathers all *.vtk files in directory with path dir_loc into a list.
           Returns a list of all the vtk files'''
        ff = []
        for file in os.listdir(dir_loc):
            if file.endswith(".vtk"):
                if file.split('.')[0].split('_')[-1] == 'pnts':
                    ff.append(os.path.join(dir_loc, file))
        ff.sort()
        self.vtk_files = ff
        return ff
           
    
    
    def read_vtk(self):
        '''Loop through all vtk files and paticle Time, Mass, and Source around each well.
           Returns a dictionary where first key is the int(model time) and the second key is the str(well location) name.
               This then holds a px3 array where each row is a particle in the control volume for the well and columns are 'Time','Mass','Source'
               -- exmaple: rtd_dict[model_time][well_name] = px3 array'''
        #pdb.set_trace()
        
        # Coordinates of the wells in Parflow space
        x = self.wells_df['X']
        
        # Well volumes from top to bottom of screen -- set up in ER_PLM_ParFlow/utils/well_info_v4.dummywells.py
        z_top = [] # top of screen
        z_bot = []
        for i in range(len(self.wells_df.index)):
            w = self.wells_df.index[i]
            if w == 'PLM7':
                # Sample Depth +/- a couple of meters for borehole PLM7
                zt = self.wells_df.loc[w,'land_surf_cx']  - self.wells_df.loc[w,'smp_depth_m'] + 2.0
                zb = self.wells_df.loc[w,'land_surf_cx']  - self.wells_df.loc[w,'smp_depth_m'] - 2.0
            else:
                # Top of screen then bottom of screen -- for piezos PLM1 and PLM6
                zt = self.wells_df.loc[w,'land_surf_cx']  - self.wells_df.loc[w,'screen_top_m']
                zb = self.wells_df.loc[w,'land_surf_cx']  - self.wells_df.loc[w,'screen_bot_m']
            z_top.append(zt)
            z_bot.append(zb)
        pf_well_coords = np.column_stack((x.to_numpy(), np.zeros(len(x.to_numpy())), z_top, z_bot))

        # Loop through each vtk file
        for i in range(len(self.vtk_files)):
            #pdb.set_trace()
            f = self.vtk_files[i]
            print ('{}/{}'.format(i+1, len(self.vtk_files)))
            
            dd = pv.read(f) 
            # xyz points of all particles -- big list
            pt_xyz = np.array(dd.points)
            
            ii = int(f.split('.')[-2])
            self.rtd[ii] = {}
            
            # Loop through each well
            for w in range(len(pf_well_coords)):
                # Draw a box around well points to collect particles near the well, meters
                xhigh = pt_xyz[:,0] <= pf_well_coords[w,0] + 2.0 #1.5125 #2.0 
                xlow  = pt_xyz[:,0] >= pf_well_coords[w,0] - 2.0 #1.5125 #2.0  
                zhigh = pt_xyz[:,2] <= pf_well_coords[w,2] # top of screen
                zlow  = pt_xyz[:,2] >= pf_well_coords[w,3] # bottom of screen
                
                msk  = np.column_stack((xlow,xhigh,zlow,zhigh))
                msk_ = msk.sum(axis=1)==4

                # Now pull values in the define box
                time   = dd['Time'][msk_]
                mass   = dd['Mass'][msk_]
                source = dd['Source'][msk_]
                Xin    = dd['xInit'][msk_]
                #Yin    = dd['yInit'][msk_]
                #pid    = dd['pid'][msk_]
                
                rtd = np.column_stack((time,mass,source,Xin))
                self.rtd[ii][well_df.index[w]] = rtd
        return self.rtd    
            
    
    
    def read_vtk_boreholes(self):
        #pdb.set_trace()
        rtd_dict = {}
        for i in range(len(times_boreholes)):
            #pdb.set_trace()
            #f = self.vtk_files[times_boreholes[i]] 
            f = self.vtk_files[i]            
            print ('b{}/{}'.format(i+1, len(times_boreholes)))
                
            ii = int(f.split('.')[-2])
            rtd_dict[ii] = {}
                
            dd = pv.read(f) 
            # xyz points of all particles -- big list
            pt_xyz = np.array(dd.points)
                
            # Well volumes from top to bottom of screen -- set up in ER_PLM_ParFlow/utils/well_info_v4.dummywells.py
            pts_z = np.column_stack(len(Xpos)*[pt_xyz[:,2]]) # helps with array comparison to avoid for looping over the different Xpos
            pts_x = np.column_stack(len(Xpos)*[pt_xyz[:,0]])

            for l in range(len(layer_elevs_)-1):
                ll = l_bls[l]
                rtd_dict[ii][ll] = {}
                rtd_df     = pd.DataFrame(index=age_list_, columns=xinds_)
                mass_df    = pd.DataFrame(index=age_list_, columns=xinds_)
                source_df  = pd.DataFrame(index=[1,2,3], columns=xinds_)
                Xin_df     = pd.DataFrame(index=np.arange(len(X_)), columns=xinds_)

                ztop  = layer_elevs_[l]
                zbot  = layer_elevs_[l+1]

                xhigh = pts_x <= Xpos + 1.5125*2 
                xlow  = pts_x >= Xpos - 1.5125*2  
                zhigh = pts_z <= ztop # top of screen
                zlow  = pts_z >= zbot # bottom of screen
                    
                msk   = np.column_stack((xlow,xhigh,zlow,zhigh))
                msk_  = msk.reshape(pts_z.shape[0], len(Xpos), 4) #msk_[:,:,0] will give 4 corresponding to first Xpos
                msk_x = msk_.sum(axis=1)==4
                        
                # Now pull values in the define box
                for j in range(len(Xpos)):
                    w = 'X{}_{}'.format(xinds[j],i)
                    
                    #pdb.set_trace()
                    
                    time     = dd['Time'][msk_x[:,j]] / 8760
                    rtd_mask = [np.where(time < a, True, False) for a in age_list_]
                    rtd_df.loc[:,xinds_[j]] = [q.sum() for q in rtd_mask] # number particles with age less than -- a CDF
                    # now the cumulative mass for each bin
                    mass     = dd['Mass'][msk_x[:,j]]
                    mass_df.loc[:,xinds_[j]] = [mass[z].sum() for z in rtd_mask]
                    
                    wt = mass/mass.sum()
                    
                    # Debugging
                    if j == 0:
                        #print ('layer {}'.format(l))
                        #print ('  {}  Number of particles'.format(len(time)))
                        #print ('  {}  Number of particles less than 1 year'.format([q.sum() for q in rtd_mask][0]))
                        #print ('  {}  Mass of particles less than 1 year'.format([mass[z].sum() for z in rtd_mask][0]))
                        #print ('  {:.4f}  Mean age (years)'.format((time*wt).sum()))
                        pass
        
                    source   = dd['Source'][msk_x[:,j]]
                    source_df.loc[1,xinds_[j]] = (source==1).sum()
                    source_df.loc[2,xinds_[j]] = (source==2).sum()
                    source_df.loc[3,xinds_[j]] = (source==3).sum()
                        
                    Xin     = dd['xInit'][msk_x[:,j]]
                    Xin_df.loc[:,xinds_[j]] = np.array([(Xin <= x).sum() for x in X_])
                        
                    rtd_dict[ii][ll]['rtd_cdf'] = rtd_df.copy()
                    rtd_dict[ii][ll]['mass']    = mass_df.copy()
                    rtd_dict[ii][ll]['source']  = source_df.copy()
                    rtd_dict[ii][ll]['Xin_df']  = Xin_df.copy()
        return rtd_dict            
            
                






#------------------
# WY2017-2021
# WY 2017-2021
#------------------
age           = ecoslim_grid_vtk(well_df)
vtk_files_c   = age.find_cgrid_vtk('./ecoslim_2017_2021')
age.vtk_files = vtk_files_c
cell_xyzm     = age.read_vtk_grid(vtk_files_c[0])
age_dict      = age.read_vtk()

get_rtd            = ecoslim_pnts_vtk(well_df, cell_xyzm)
vtk_files          = get_rtd.find_pnts_vtk('./ecoslim_2017_2021')
get_rtd.vtk_files  = vtk_files[::5]
rtd_dict           = get_rtd.read_vtk()
# Save to a dictionary
with open('./parflow_out/ecoslim_MeanAge.1721.pk', 'wb') as f:
    pickle.dump(age_dict, f) 
with open('./parflow_out/ecoslim_rtd.1721.pk', 'wb') as ff:
    pickle.dump(rtd_dict, ff) 






"""
#---------------------------------------
#
# Vertical Age Dynamics in Boreholes
#
#---------------------------------------

# Parflow variable dz
dz = np.array([1.00, 1.00, 1.00, 1.00, 1.00,       # 52.0 - 102.0
               0.80, 0.80,                         # 36.0 - 52.0
               0.60, 0.60,                         # 24.0 - 36.0
               0.40, 0.40,                         # 16.0 - 24.0
               0.20, 0.20, 0.20,                   # 10.0 - 16.0  -- 2m layers down to 16 m
               0.10, 0.10,                         # 8.0  - 10.0
               0.10, 0.05, 0.05,                   # 6.0  - 8.0   -- 0.5m res possible down to 7.0 m
               0.05, 0.05, 0.05, 0.05,             # 4.0  - 6.0
               0.05, 0.05, 0.05, 0.05,             # 2.0  - 4.0
               0.05, 0.05, 0.05, 0.025, 0.025])    # 0.0  - 2.0  
dz_scale = 10 * dz


directory = 'wy_2017_2021' # hard code for now
header    = 'wy_2017_2021' # hard code for now
fn_porosity = os.path.join(directory, '{}.out.{}.pfb'.format(header,'porosity'))
por = np.flipud(read_pfb(fn_porosity)[:,0,:])

soil_ind = np.where(por[:,0]==por.max())[0].max()+1
bed_ind  = np.where(por[:,0]==por.min())[0].min()

# DEM info
Z_  = np.loadtxt('../ER_PLM_ParFlow/elevation_v4.sa', skiprows=1)
#X_  = np.arange(1., len(Z_))*1.5125 
X_  = np.linspace(1.5125, len(Z_)*1.5125, len(Z_))
dem = np.column_stack((X_,np.zeros_like(X_), Z_))

xinds = [265, 424, 494, 508]
xinds_ = ['X265', 'X424', 'X494', 'X508']
Xpos  = X_[xinds] 

layer_depths = np.flip(dz_scale).cumsum()
layer_elevs  = np.array([np.subtract(Z_[xinds][i], np.flip(dz_scale).cumsum()) for i in range(len(xinds))]).T # depth to bottom of each layer in altitude (m)

soil_elevs = layer_elevs[:soil_ind, :] # note, this does not include land surface elevation
sap_elevs  = layer_elevs[soil_ind:bed_ind, :]
bed_elevs  = layer_elevs[bed_ind:, :]

# The above is very slow with N=32 layers. Manually define some layers instead
#l_bls = np.array([0.5, 2.0, 4.5, 6.5, 9.0, 14.0, 20.0, 30.0, 50.0, 80.0, 102.0])
#layer_elevs_ = np.array([np.subtract(Z_[xinds][i], l_bls) for i in range(len(xinds))]).T
l_bls = layer_depths.copy()
layer_elevs_ = layer_elevs.copy()


# Create a list of ages to calculate CDFs
#age_list_ = np.concatenate((np.arange(1,10,1), np.arange(10,50,5), np.arange(50,1000,10), np.arange(1000, 35000, 100)))
age_list_ = np.concatenate((np.arange(0.5,5.0,0.5), np.arange(5,20,1), np.arange(20,100,5), np.arange(100,1000,10), np.arange(1000, 35000, 100)))



# Create modified timelists to speed up process
def pf_2_dates(startdate, enddate, f):
    '''Assumes ParFlow outputs every 24 hours'''
    s = pd.to_datetime(startdate)
    e = pd.to_datetime(enddate)
    d_list = pd.date_range(start=s, end=e, freq=f)
    # Drop Leap years again
    d_list_ = d_list[~((d_list.month == 2) & (d_list.day == 29))]
    return d_list_

def set_wy(df):
    dates     = df.copy().index
    yrs       = dates.year
    yrs_      = np.unique(yrs)[1:]
    wy_inds_  = [np.where((dates > '{}-09-30'.format(i-1)) & (dates < '{}-10-01'.format(i)), True, False) for i in yrs_]
    wy_inds   = np.array([wy_inds_[i]*yrs_[i] for i in range(len(yrs_))]).sum(axis=0)
    first_yrs = [(wy_inds==i).argmax() for i in yrs_]
    return list(yrs_), list(first_yrs)


dates     = pf_2_dates('2016-09-30', '2021-08-29', '24H')
yrs       = dates.year
yrs_      = np.unique(yrs)[1:]
wy_inds_  = [np.where((dates > '{}-09-30'.format(i-1)) & (dates < '{}-10-01'.format(i)), True, False) for i in yrs_]
wy_inds   = np.array([wy_inds_[i]*yrs_[i] for i in range(len(yrs_))]).sum(axis=0)
#times_boreholes = np.where((dates.day==1) & (wy_inds==2017))[0]
# more targeted dates
date_boreholes = ['2016-10-01', '2016-11-01', '2017-01-01', '2017-03-01',
                  '2017-05-01', '2017-05-10', '2017-05-20', 
                  '2017-06-01', '2017-06-10', '2017-06-20',
                  '2017-07-01', '2017-07-15',
                  '2017-08-01', '2017-09-01']
times_boreholes = np.array([np.where(dates==db_)[0][0] for db_ in date_boreholes])




# Finally, run it
get_rtd            = ecoslim_pnts_vtk(well_df, cell_xyzm)
vtk_files          = get_rtd.find_pnts_vtk('./ecoslim_2017_2021')
rtd_dict_b         = get_rtd.read_vtk_boreholes()
# Save to a dictionary
with open('./parflow_out/ecoslim_rtd.1721.bores.pk', 'wb') as ff:
    pickle.dump(rtd_dict_b, ff) 

"""











        

        






##------
## WY2000-2016
#age           = ecoslim_grid_vtk(well_df)
#vtk_files_c   = age.find_cgrid_vtk('./ecoslim_2000_2016')
#age.vtk_files = vtk_files_c
#cell_xyzm     = age.read_vtk_grid(vtk_files_c[0])
#age_dict      = age.read_vtk()
#
#get_rtd            = ecoslim_pnts_vtk(well_df, cell_xyzm)
#vtk_files          = get_rtd.find_pnts_vtk('./ecoslim_2000_2016')
#get_rtd.vtk_files  = vtk_files[::5]
#rtd_dict           = get_rtd.read_vtk()
## Save to a dictionary
#with open('./parflow_out/ecoslim_MeanAge.0016.pk', 'wb') as f:
#    pickle.dump(age_dict, f)            
#with open('./parflow_out/ecoslim_rtd.0016.pk', 'wb') as ff:
#    pickle.dump(rtd_dict, ff)  







#def flux_wt_rtd(rtd_dict_unsort, model_time, well_name, nbins):
#    '''Convert particle ages at a single model time and well location to a residence time distribution.
#    Returns:
#        - rtd_df:   Dataframe with mass weighted ages for all paricles (sorted by age)
#        - rtd_dfs:  Dataframe similar to above but bins age distribution into discrete intervals  
#    Inputs: 
#        - rtd_dict_unsort: output from ecoslim_pnts_vtk.read_vtk() above
#        - model_time: model time to consider. Must be in rtd_dict_unsort.keys()
#        - well_name: observation well to consider. Must be in rtd_dict_unsort.keys()
#        - nbins: Number of intervals to bin the ages into.'''
#    # Flux Weighted RTD
#    # Info regarding particles at a single timestep and single point (box) of the domain
#    rtd    = rtd_dict_unsort[model_time][well_name] 
#    rtd_df = pd.DataFrame(data=rtd,columns=['Time','Mass','Source'])
#    rtd_df['wt'] = rtd_df['Mass']/rtd_df['Mass'].sum()
#    rtd_df.sort_values('Time', inplace=True)
#    rtd_df['Time'] /= 8760
#    
#    # Now some binning
#    #nbins = 10
#    gb =  rtd_df.groupby(pd.cut(rtd_df['Time'], nbins))
#    rtd_dfs = gb.agg(dict(Time='mean',Mass='sum',Source='mean',wt='sum'))
#    rtd_dfs['count'] = gb.count()['Time']
#
#    return rtd_df, rtd_dfs
#
#rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, 202, 'PLM6', 10)


## Plot the RTD
#fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
##ax[0].plot(rtd_df1['Time'], rtd_df1['wt'], marker='.')
#ax[0].plot(rtd_dfs['Time'], rtd_dfs['wt'], marker='.', color='red')
#ax[0].set_ylabel('Probability')
#ax[0].set_xlabel('Travel Time (years)')
#
#ax[1].plot(rtd_df['Time'], np.cumsum(rtd_df['wt']), marker='.')
#ax[1].plot(rtd_dfs['Time'], np.cumsum(rtd_dfs['wt']), marker='.', color='red')
#ax[1].set_ylabel('Probability CDF')
#ax[1].set_xlabel('Travel Time (years)')
#fig.tight_layout()
#plt.show()





    
 
