# Low K bedrock version of Run127

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker


#----------
# Read in top surface DEM files
#head  = np.loadtxt('./elevation.sa', max_rows=1)
#dem   = np.loadtxt('./elevation.sa', skiprows=1)


savename = 'subsurface_00127' 

#-----------
# ParFlow Grid
LowerX = 0.0
LowerY = 0.0
LowerZ = 0.0

NX = 559
NY = 1
NZ = 32

DX = 1.5125
DY = 1.5125
DZ = 10.0



#------------
# Variable Z scaling
# 							                              # Depth below land surface
dz_scale = np.array([1.00, 1.00, 1.00, 1.00, 1.00,       # 52.0 - 102.0
                     0.80, 0.80,                         # 36.0 - 52.0
                     0.60, 0.60,                         # 24.0 - 36.0
                     0.40, 0.40,                         # 16.0 - 24.0
                     0.20, 0.20, 0.20,                   # 10.0 - 16.0  -- 2m layers down to 16 m
                     0.10, 0.10,                         # 8.0  - 10.0
                     0.10, 0.05, 0.05,                   # 6.0  - 8.0   -- 0.5m res possible down to 7.0 m
                     0.05, 0.05, 0.05, 0.05,             # 4.0  - 6.0
                     0.05, 0.05, 0.05, 0.05,             # 2.0  - 4.0
                     0.05, 0.05, 0.05, 0.025, 0.025])    # 0.0  - 2.0  
NZ = len(dz_scale)


# Layer thickness, counting from bottom up
dz_scaled = np.array(dz_scale) * DZ
# Cumulative thickness, counting from bottom up
dz_cumsum = dz_scaled.cumsum()

# Depth below land surface, cumulative thickness
# Counts top down, because makes more sense to me
dbs = np.flip(dz_scaled).cumsum() 



#---------------------------
# Update 10/15/21
# Treat all 32 layers as own indicator index
# Keep track of permeability per layer

# Update 10/19/21
# Indicator layers through 16m bls are 2 m thick
# These can include multiple mesh layers


#-------------------------------
# Set Statigraphy Properties
# CHANGE ME
#------------------------------

# Write depth to bottom of layer below land surface
# top starts at zero 
layers_dbs = {}
layers_dbs['soil1']   = 1.0
layers_dbs['wshale']  = 4.0
layers_dbs['fshale']  = 102.0 # Bottom of domain

# Saturated hydraulic conductivity in m/s
layers_K = {}
layers_K['soil1']   = 2.5523e-04  # middle between top soil and sub soil in tokonaga 2022
layers_K['wshale']  = 5.070914231103997e-05  
#layers_K['fshale']  = 8.130408338859933e-08  # tokonaga 2022 only gives a tranmisivity value
layers_K['fshale']  = 2.5710608658023473e-08  # half order of magnitude smaller

# Porosity
layers_por = {}
layers_por['soil1']   = 0.40
layers_por['wshale']  = 0.20
layers_por['fshale']  = 0.10

# Van Genutchen relative perm: [alpha,n]
layers_rel = {}
layers_rel['soil1']   = [4.415024485, 4.817367116]
layers_rel['wshale']  = [4.127707749, 3.156116411]
layers_rel['fshale']  = [0.519, 1.595]

# Van Genutchen saturation: [salpha,sn,sres,ssat]
layers_srel = {}
layers_srel['soil1']   = layers_rel['soil1']  + [0.131, 1.0]  #[1.820, 1.789, 0.131, 1.0]
layers_srel['wshale']  = layers_rel['wshale'] + [0.04, 1.0]   #[0.519, 1.595, 0.04,  1.0]
layers_srel['fshale']  = layers_rel['fshale'] + [0.005, 1.0]  #[0.519, 1.595, 0.005, 1.0]



#--------------------------
# Dataframe holding row->indicators and material properties
# Indicator thickness are 2m down to 16 m depth
layers_dbs_nam = list(layers_dbs.keys())
layers_dbs_arr = np.array(list(layers_dbs.values()))

# Map index to layers
ind_2_layer = pd.DataFrame()
ind_2_layer['bls'] = dbs

# Where to make the cuts
#c = (dbs%2 == 0.).astype(float)
#inds = []
#ind_cc = 1.0
#for i in range(len(c)):
#    if c[i] == 0.0:
#        inds.append(ind_cc)
#    elif c[i] == 1.0:
#        inds.append(ind_cc)
#        ind_cc += 1.0

ind_2_layer['ind_num'] = ind_2_layer.index+1

# Assing a lithology to the layers
cc = 0
for i in range(len(layers_dbs_arr)):
    lo_ind = cc
    hi_ind = (dbs == layers_dbs_arr[i]).argmax()+1
    cc = hi_ind
    #ind_2_layer[lo_ind:hi_ind] = layers_dbs_nam[i]
    ind_2_layer.loc[lo_ind:hi_ind, 'layer'] = layers_dbs_nam[i]


# Add in K and porosity
for i in list(layers_K.keys()):
    # m/s
    ind_2_layer.loc[ind_2_layer['layer']==i, 'K_ms'] = layers_K[i]
    # m/hr
    ind_2_layer.loc[ind_2_layer['layer']==i, 'K_mhr'] = layers_K[i]*3600
    # porosity
    ind_2_layer.loc[ind_2_layer['layer']==i, 'porosity'] = layers_por[i]


# Add in Van Genutchen parameters
for i in list(layers_rel.keys()):
    # Relative Perm
    ind_2_layer.loc[ind_2_layer['layer']==i, 'vg_alpha'] = layers_rel[i][0]
    ind_2_layer.loc[ind_2_layer['layer']==i, 'vg_n'] = layers_rel[i][1]
    # Sautration
    ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_alpha'] = layers_srel[i][0]
    ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_n'] = layers_srel[i][1]
    ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_res'] = layers_srel[i][2]
    ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_sat'] = layers_srel[i][3]




#-------------
# Post process for exponentially decreasing K in fshale
fshale_inds = ind_2_layer['layer']=='fshale' # df indices to make life easier
# depth below whatever layer of interest
dbs_shale = ind_2_layer.loc[fshale_inds,'bls'].to_numpy()


def K_exp_fun(K0, K1, depth_prof):
    # solve system k0=A*e^(lam*z0) and k1=A*e^(lam*z1) for lam and A
    z0,z1 = depth_prof[0], depth_prof[-1]
    lam = 1/(z0-z1) *  np.log(K0/K1)
    A = K1*np.exp(-lam*z1)

    K_exp = A*np.exp(lam*dbs_shale)
    return K_exp

"""
# quick plot
fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10,6))
for Kl in [1.e-12, 1.e-10, 1.e-8, 1.e-6]:
    K_exp = K_exp_fun(1.e-6, Kl, dbs_shale)
    ax[0].plot(K_exp/K_exp.max(), dbs_shale, marker='.')
    ax[1].plot(K_exp, dbs_shale, marker='.')
ax[0].invert_yaxis()
ax[1].invert_yaxis()
ax[1].set_xscale('log')

ax[0].set_ylabel('Depth (m)')
ax[0].set_xlabel(r'K/K$_{max}$')
ax[1].set_xlabel('K (m/s)')
[ax[i].yaxis.set_major_locator(ticker.MultipleLocator(10)) for i in [0,1]]
[ax[i].yaxis.set_minor_locator(ticker.MultipleLocator(5)) for i in [0,1]]
plt.show()
"""



#------------
# CHANGE ME
#------------
# pick one to model
Kshale_top      = layers_K['fshale']
Kshale_bottom   = layers_K['fshale']
#Kshale_bottom   = 1.e-8

K_exp = K_exp_fun(Kshale_top, Kshale_bottom, dbs_shale)


# now change dataframe
ind_2_layer['K_mhr_dec'] = ind_2_layer['K_ms'].copy()*3600
ind_2_layer.loc[fshale_inds, 'K_mhr_dec'] = K_exp*3600

# clean up
ind_2_layer.drop(columns=['K_ms','K_mhr'], inplace=True)

#-------
# Finally, group by ind_num
#ind_2_layer_s = ind_2_layer.groupby('ind_num').agg({'bls':'max',
#                                    'porosity':'mean',
#                                    'vg_alpha':'mean',
#                                    'vg_n':'mean',
#                                    'vgSat_alpha':'mean',
#                                    'vgSat_n':'mean',
#                                    'vgSat_res':'mean',
#                                    'vgSat_sat':'mean',
#                                    'K_mhr_dec':'mean'})    





#------------------
# make a grid
x,z = np.meshgrid(np.arange(NX), ind_2_layer['ind_num'].to_numpy())

# Make a layer-cake permeability field
# back to index 0 being bottom left of domain
z_ = np.flip(z).ravel()

#---------------
# Write an SA file
# Each row has own indicator number
np.savetxt('{}.sa'.format(savename), z_, fmt='%.4e', header='559 1 {}'.format(NZ), comments='')




#-----------------------
# Start to write ParFlow tcl 
# Includes modifying the template
# and writing to external file

fname = '{}.mat'.format(savename)

# some meta data first
f = open(fname, 'w')
f.writelines('#Material Property Information\n\n')


### variable dz assignments
f.writelines('#-------------------------\n')
f.writelines('# Variable dz Assignment\n')
f.writelines('#-------------------------\n')
f.writelines('pfset Solver.Nonlinear.VariableDz       True\n')
f.writelines('pfset dzScale.GeomNames                 domain\n')
f.writelines('pfset dzScale.Type                      nzList\n')
f.writelines('pfset dzScale.nzListNumber\t\t{}\n\n'.format(len(dz_scale)))
for i in range(len(dz_scale)):
    l = 'pfset Cell.{}.dzScale.Value\t\t{}\n'.format(i, dz_scale[i])
    f.writelines(l)


### indicator geometry assigments
f.writelines('\n#-------------------------\n')
f.writelines('# Indicator Geometry Input\n')
f.writelines('#-------------------------\n')
# indicator list
indlist = ['{}{}'.format('i', zz) for zz in ind_2_layer['ind_num']]
f.writelines('set indlist "{}"\n\n'.format(' '.join(indlist)))
# some more stuff
f.writelines('pfset GeomInput.indi_input.InputType      IndicatorField\n')
f.writelines('pfset GeomInput.indi_input.GeomNames      $indlist\n')
f.writelines('pfset Geom.indi_input.FileName            $ind_fname\n\n')
# assign indicator numbers 
for i in range(1,len(ind_2_layer)+1):
    f.writelines('pfset GeomInput.i{}.Value\t\t{}\n'.format(i,i))


### permeability assinments, in m/hr for each layer
f.writelines('\n#-------------------------\n')
f.writelines('# Permeability\n')
f.writelines('#-------------------------\n')
f.writelines('pfset Geom.Perm.Names                  "domain $indlist"\n\n')
f.writelines('pfset Geom.domain.Perm.Type             Constant\n')
f.writelines('pfset Geom.domain.Perm.Value            0.01\n\n')
# Loop through permeability dictionary
for i in range(len(ind_2_layer)):
    f.writelines('pfset Geom.i{}.Perm.Type \t\t Constant\n'.format(i+1))
    f.writelines('pfset Geom.i{}.Perm.Value \t\t {:.10f}\n\n'.format(i+1, ind_2_layer.loc[i,'K_mhr_dec']))


### porosity
f.writelines('\n#-------------------------\n')
f.writelines('# Porosity\n')
f.writelines('#-------------------------\n')
f.writelines('pfset Geom.Porosity.GeomNames          "domain $indlist"\n\n')
f.writelines('pfset Geom.domain.Porosity.Type  \t\t Constant\n')
f.writelines('pfset Geom.domain.Porosity.Value \t\t 0.3\n\n')
# Loop through porosoity dictionary
for i in range(len(ind_2_layer)):
    f.writelines('pfset Geom.i{}.Porosity.Type \t\t Constant\n'.format(i+1))
    f.writelines('pfset Geom.i{}.Porosity.Value \t\t {:.4f}\n\n'.format(i+1, ind_2_layer.loc[i,'porosity']))


### relative perm
f.writelines('\n#-------------------------\n')
f.writelines('# Relative Permeability\n')
f.writelines('#-------------------------\n')
f.writelines('pfset Phase.RelPerm.Type                  VanGenuchten\n')
f.writelines('pfset Phase.RelPerm.GeomNames             "domain $indlist"\n\n')
f.writelines('pfset Geom.domain.RelPerm.Alpha \t\t 3.5\n')
f.writelines('pfset Geom.domain.RelPerm.N     \t\t 2.0\n\n')
# Loop through porosoity dictionary
for i in range(len(ind_2_layer)):
    f.writelines('pfset Geom.i{}.RelPerm.Alpha \t\t {:.5f}\n'.format(i+1, ind_2_layer.loc[i,'vg_alpha']))
    f.writelines('pfset Geom.i{}.RelPerm.N \t\t {:.5f}\n\n'.format(i+1, ind_2_layer.loc[i,'vg_n']))

### relative saturation
f.writelines('\n#-------------------------\n')
f.writelines('# Saturation Functionn\n')
f.writelines('#-------------------------\n')
f.writelines('pfset Phase.Saturation.Type                  VanGenuchten\n')
f.writelines('pfset Phase.Saturation.GeomNames             "domain $indlist"\n\n')
f.writelines('pfset Geom.domain.Saturation.Alpha           3.5\n')
f.writelines('pfset Geom.domain.Saturation.N               2.0\n')
f.writelines('pfset Geom.domain.Saturation.SRes            0.2\n')
f.writelines('pfset Geom.domain.Saturation.SSat            1.0\n\n')
# Loop through porosoity dictionary
for i in range(len(ind_2_layer)):
    f.writelines('pfset Geom.i{}.Saturation.Alpha \t\t {:.5f}\n'.format(i+1, ind_2_layer.loc[i,'vgSat_alpha']))
    f.writelines('pfset Geom.i{}.Saturation.N \t\t {:.5f}\n'.format(i+1, ind_2_layer.loc[i,'vgSat_n']))
    f.writelines('pfset Geom.i{}.Saturation.SRes \t\t {:.5f}\n'.format(i+1, ind_2_layer.loc[i,'vgSat_res']))
    f.writelines('pfset Geom.i{}.Saturation.SSat \t\t {:.5f}\n\n'.format(i+1, ind_2_layer.loc[i,'vgSat_sat']))
    
f.writelines('#END\n\n')
f.close()









"""#-----------------------    
# Summarize into some files and dataframe
sub_map = pd.DataFrame()

# indicator numbers
sub_map['ind_ID'] = p_map_.astype(int)
# indicator names from dz_layers
l = [np.flip(list(dz_layers.keys()))[int(i)-1] for i in sub_map['ind_ID']]
sub_map['layer'] = l
# thickness of each layer
sub_map['dz_m'] = np.flip(dz_scale)*DZ
# depth of each layer bls
sub_map['bls_m']  = dbs
# total depth for indicator ID #
ld = [np.flip(list(dz_layers.values()))[int(i)-1] for i in sub_map['ind_ID']]
sub_map['layer_td'] = ld

# Save full to csv
sub_map.to_csv('z_map.txt', index=True, sep='\t')

# something easier to read
sub_map_ = sub_map.groupby('ind_ID').sum()['dz_m'].cumsum()
sub_map_.to_csv('z_maps.txt', index=True, sep='\t')




# Indicator Info For ParFlow Script
ind_id   = np.arange(1,len(dz_layers)+1)
ind_name = list(dz_layers.keys())"""











"""
#----------
# Read in original layering
# Use this for comparison only
p_og  = np.loadtxt('./subsurface_transect_regrid9_extended.sa', skiprows=1)
p_og = np.flipud(p_og.reshape(NZ,NX))
p_og_cake = np.repeat(p_og[:,0],p_og.shape[1]).reshape(NZ,NX)
p_og_cake_ = np.flip(p_og_cake.ravel())

np.savetxt('subsurface_lcog.sa', p_og_cake_, fmt='%.4e', header='559 1 26', comments='')
"""

