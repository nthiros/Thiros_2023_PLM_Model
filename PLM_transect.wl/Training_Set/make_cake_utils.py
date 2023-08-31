import numpy as np
import pandas as pd
import sys
import subprocess
import pdb


class pf_properties():
    def __init__(self, layer_bls, dz_scale):
        self.dbs = layer_bls  # depth to bottom of all layers. Starts at land surface then counts down
                              # np.flip(np.array(dz_scale)*DZ).cumsum()
        self.dz_scale = dz_scale
        
        # Uncertain Parameters to sample over
        self.layers_dbs    =  None  # depth below land surface for each layer
        self.layers_K      =  None  # saturated hydraulic conductivity (m/s) for the soil, w.shale, and f.shale
        self.layers_por    =  None  # porosity as decimal
        self.layers_krel   =  None  # van genutchen relative permeability func [alpha, n]
        self.layers_srel   =  None  # van genutchen relative saturation func [salpha, sn, sres, ssat]

        # Holds dataframe to written at the end after monte carlo sampling
        self.ind_2_layer   =  None  # dataframe that holds all the material properties map to the correct layers


    #-------------------------------------
    # Set Statigraphy and Layer Properties
    #-------------------------------------
    def set_layer_depth(self, d_soil, d_wshale):
        ''' Write depth to bottom of layer below land surface. Top starts at zero 
            Assumes a soil, weathered shale, and fractured shale right now.
            Bottom of fractured shale is bottom of domain'''
        layers_dbs = {}
        layers_dbs['soil']    = d_soil
        layers_dbs['wshale']  = d_wshale
        layers_dbs['fshale']  = self.dbs.max() # Bottom of domain
        self.layers_dbs = layers_dbs
        
    def set_layer_K(self, K_soil, K_wshale, K_fshale):
        '''Saturated hydraulic conductivity in m/s'''
        layers_K = {}
        layers_K['soil']    = K_soil
        layers_K['wshale']  = K_wshale 
        layers_K['fshale']  = K_fshale
        self.layers_K = layers_K
        
    def set_layer_por(self, por_soil, por_wshale, por_fshale):
        '''Porosity'''
        layers_por = {}
        layers_por['soil']    = por_soil
        layers_por['wshale']  = por_wshale
        layers_por['fshale']  = por_fshale
        self.layers_por = layers_por
        
    def set_vg_perm(self, soil_krel, wshale_krel, fshale_krel):
        '''Van Genutchen relative perm
           Takes a list of [alpha, n]'''
        layers_krel = {}
        layers_krel['soil']    = soil_krel
        layers_krel['wshale']  = wshale_krel
        layers_krel['fshale']  = fshale_krel
        self.layers_krel = layers_krel
    
    def set_vg_sat(self, soil_srel, whsale_srel, fshale_srel):
        '''Van Genutchen saturation parameters
           Takes list of :[salpha,sn,sres,ssat]'''
        layers_srel = {}
        layers_srel['soil']    = soil_srel
        layers_srel['wshale']  = whsale_srel
        layers_srel['fshale']  = fshale_srel
        self.layers_srel = layers_srel


    #---------------------------------------------
    # Dataframe holding row->indicators mapping 
    # and material properties
    #---------------------------------------------
    def build_mats_df(self):
        # Each row in df is a layer
        #pdb.set_trace()
        layers_dbs_nam = list(self.layers_dbs.keys())
        layers_dbs_arr = np.array(list(self.layers_dbs.values()))
        
        # Map index to layers
        ind_2_layer = pd.DataFrame()
        ind_2_layer['bls'] = self.dbs
        ind_2_layer['ind_num'] = ind_2_layer.index+1
        
        # Assigning a lithology to the layers
        cc = 0
        for i in range(len(layers_dbs_arr)):
            lo_ind = cc
            hi_ind = (self.dbs == layers_dbs_arr[i]).argmax()+1
            cc = hi_ind
            ind_2_layer.loc[lo_ind:hi_ind, 'layer'] = layers_dbs_nam[i]
        
        # Add in K and porosity
        for i in list(self.layers_K.keys()):
            # m/s
            ind_2_layer.loc[ind_2_layer['layer']==i, 'K_ms']  = self.layers_K[i]
            # m/hr
            ind_2_layer.loc[ind_2_layer['layer']==i, 'K_mhr'] = self.layers_K[i]*3600
            # porosity
            ind_2_layer.loc[ind_2_layer['layer']==i, 'porosity'] = self.layers_por[i]
        
        # Add in Van Genutchen parameters
        for i in list(self.layers_krel.keys()):
            # Relative Perm
            ind_2_layer.loc[ind_2_layer['layer']==i, 'vg_alpha'] = self.layers_krel[i][0]
            ind_2_layer.loc[ind_2_layer['layer']==i, 'vg_n']     = self.layers_krel[i][1]
            # Sautration
            ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_alpha'] = self.layers_srel[i][0]
            ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_n']     = self.layers_srel[i][1]
            ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_res']   = self.layers_srel[i][2]
            ind_2_layer.loc[ind_2_layer['layer']==i, 'vgSat_sat']   = self.layers_srel[i][3]
    
        self.ind_2_layer = ind_2_layer
        return ind_2_layer


    #---------------------------------
    # Post process for exponentially 
    # decreasing K in fshale
    #---------------------------------
    def K_exp_fun(self, K0, K1, depth_prof):
        '''solve system k0=A*e^(lam*z0) and k1=A*e^(lam*z1) for lam and A.
           K0 is conductivity at top, K1 is conductivity at bottom, depth_prof is layering dbs'''
        z0,z1 = depth_prof[0], depth_prof[-1]
        lam = 1/(z0-z1) * np.log(K0/K1)
        A = K1*np.exp(-lam*z1)
        K_exp = A*np.exp(lam*depth_prof)
        return K_exp

    def apply_exp_K(self, Ktop, Kbot):
        '''Ktop is K at the wshale:fshale boundary, Kbot is K at bottom of domain.
           Set Ktop=Kbot for no decay.
           Setting Ktop and Kbot controls the exponential K decay rate.'''
        # This only considers exponential K in fshale layer
        fshale_inds = self.ind_2_layer['layer'] == 'fshale' 
        dbs_shale   = self.ind_2_layer.loc[fshale_inds,'bls'].to_numpy()
        
        #Kshale_top      = Ktop
        #Kshale_bottom   = Kbot
        
        K_exp = self.K_exp_fun(Ktop, Kbot, dbs_shale)
        
        # Update dataframe
        self.ind_2_layer['K_mhr_dec'] = self.ind_2_layer['K_ms'].copy()*3600
        self.ind_2_layer.loc[fshale_inds, 'K_mhr_dec'] = K_exp*3600
        self.ind_2_layer.drop(columns=['K_ms','K_mhr'], inplace=True) # cleanup

    #----------------------------------
    # Write ParFlow tcl materials
    # to external file
    #----------------------------------
    def write_mat_pfb(self, savenum):
        fname = 'subsurface_{}.mat'.format(savenum)
        #
        # some meta data first
        f = open(fname, 'w')
        f.writelines('#Material Property Information\n\n')
        #
        # variable dz assignments
        f.writelines('#-------------------------\n')
        f.writelines('# Variable dz Assignment\n')
        f.writelines('#-------------------------\n')
        f.writelines('pfset Solver.Nonlinear.VariableDz       True\n')
        f.writelines('pfset dzScale.GeomNames                 domain\n')
        f.writelines('pfset dzScale.Type                      nzList\n')
        f.writelines('pfset dzScale.nzListNumber\t\t{}\n\n'.format(len(self.dz_scale)))
        for i in range(len(self.dz_scale)):
            l = 'pfset Cell.{}.dzScale.Value\t\t{}\n'.format(i, self.dz_scale[i])
            f.writelines(l)
        #
        # indicator geometry assigments
        f.writelines('\n#-------------------------\n')
        f.writelines('# Indicator Geometry Input\n')
        f.writelines('#-------------------------\n')
        # indicator list
        indlist = ['{}{}'.format('i', zz) for zz in self.ind_2_layer['ind_num']]
        f.writelines('set indlist "{}"\n\n'.format(' '.join(indlist)))
        # some more stuff
        f.writelines('pfset GeomInput.indi_input.InputType      IndicatorField\n')
        f.writelines('pfset GeomInput.indi_input.GeomNames      $indlist\n')
        f.writelines('pfset Geom.indi_input.FileName            $ind_fname\n\n')
        # assign indicator numbers 
        for i in range(1,len(self.ind_2_layer)+1):
            f.writelines('pfset GeomInput.i{}.Value\t\t{}\n'.format(i,i))
        #
        # permeability assinments, in m/hr for each layer
        f.writelines('\n#-------------------------\n')
        f.writelines('# Permeability\n')
        f.writelines('#-------------------------\n')
        f.writelines('pfset Geom.Perm.Names                  "domain $indlist"\n\n')
        f.writelines('pfset Geom.domain.Perm.Type             Constant\n')
        f.writelines('pfset Geom.domain.Perm.Value            0.01\n\n')
        # Loop through permeability dictionary
        for i in range(len(self.ind_2_layer)):
            f.writelines('pfset Geom.i{}.Perm.Type \t\t Constant\n'.format(i+1))
            f.writelines('pfset Geom.i{}.Perm.Value \t\t {:.10f}\n\n'.format(i+1, self.ind_2_layer.loc[i,'K_mhr_dec']))
        #
        # porosity
        f.writelines('\n#-------------------------\n')
        f.writelines('# Porosity\n')
        f.writelines('#-------------------------\n')
        f.writelines('pfset Geom.Porosity.GeomNames          "domain $indlist"\n\n')
        f.writelines('pfset Geom.domain.Porosity.Type  \t\t Constant\n')
        f.writelines('pfset Geom.domain.Porosity.Value \t\t 0.3\n\n')
        # Loop through porosoity dictionary
        for i in range(len(self.ind_2_layer)):
            f.writelines('pfset Geom.i{}.Porosity.Type \t\t Constant\n'.format(i+1))
            f.writelines('pfset Geom.i{}.Porosity.Value \t\t {:.5f}\n\n'.format(i+1, self.ind_2_layer.loc[i,'porosity']))
        #
        # relative perm
        f.writelines('\n#-------------------------\n')
        f.writelines('# Relative Permeability\n')
        f.writelines('#-------------------------\n')
        f.writelines('pfset Phase.RelPerm.Type                  VanGenuchten\n')
        f.writelines('pfset Phase.RelPerm.GeomNames             "domain $indlist"\n\n')
        f.writelines('pfset Geom.domain.RelPerm.Alpha \t\t 3.5\n')
        f.writelines('pfset Geom.domain.RelPerm.N     \t\t 2.0\n\n')
        # Loop through porosoity dictionary
        for i in range(len(self.ind_2_layer)):
            f.writelines('pfset Geom.i{}.RelPerm.Alpha \t\t {:.5f}\n'.format(i+1, self.ind_2_layer.loc[i,'vg_alpha']))
            f.writelines('pfset Geom.i{}.RelPerm.N \t\t {:.5f}\n\n'.format(i+1, self.ind_2_layer.loc[i,'vg_n']))
        #
        # relative saturation
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
        for i in range(len(self.ind_2_layer)):
            f.writelines('pfset Geom.i{}.Saturation.Alpha \t\t {:.5f}\n'.format(i+1, self.ind_2_layer.loc[i,'vgSat_alpha']))
            f.writelines('pfset Geom.i{}.Saturation.N \t\t {:.5f}\n'.format(i+1, self.ind_2_layer.loc[i,'vgSat_n']))
            f.writelines('pfset Geom.i{}.Saturation.SRes \t\t {:.5f}\n'.format(i+1, self.ind_2_layer.loc[i,'vgSat_res']))
            f.writelines('pfset Geom.i{}.Saturation.SSat \t\t {:.5f}\n\n'.format(i+1, self.ind_2_layer.loc[i,'vgSat_sat']))
            
        f.writelines('#END\n\n')
        f.close()


    #------------------
    # make a grid
    #------------------
    def make_grid(self, NX, savenum):
        '''Writes an indicator.sa file.
           NX is number of dimensions in x direction
           savenum is a unique file number'''
        x,z = np.meshgrid(np.arange(NX), self.ind_2_layer['ind_num'].to_numpy())
        # Make a layer-cake permeability field
        # back to index 0 being bottom left of domain
        z_ = np.flip(z).ravel()
        #
        # Write an SA file
        # Save each row has own indicator number
        np.savetxt('subsurface_{}.sa'.format(savenum), z_, fmt='%.4e', header='{} 1 {}'.format(int(NX), int(len(self.dbs))), comments='')
        # Write to pfb file
        call = ['tclsh','file_conversion.sa2pfb.tcl','subsurface_{}.sa'.format(savenum)]
        subprocess.run(['/usr/bin/tclsh8.5','file_conversion.sa2pfb.tcl','subsurface_{}.sa'.format(savenum)])
        # Write material properties  
        self.write_mat_pfb(savenum)






"""
Moved to make_cake_sobol.py script
#----------------------------------------
# ParFlow Grid
#
LowerX = 0.0
LowerY = 0.0
LowerZ = 0.0

NX = 559
NY = 1
NZ = 32

DX = 1.5125
DY = 1.5125
DZ = 10.0

#----------------------------------------
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
# Depth below land surface, cumulative thickness
# Counts top down, because makes more sense to me
dbs = np.flip(dz_scaled).cumsum() 


#----------------------------------------
# Set hillslope properties
# 
savenum = 1

pf = pf_properties(dbs, dz_scale)
# Setup Properties -- Vary these with LHS sampling
d_soil, d_wshale                    = 5.0, 9.0
K_soil, K_wshale, K_fshale          = 5.605e-05, 5.605e-06, 5.605e-08
por_soil, por_wshale, por_fshale    = 0.40, 0.10, 0.05
soil_krel, wshale_krel, fshale_krel = [1.820,1.789], [0.519,1.595], [0.519,1.595]
soil_srel, whsale_srel, fshale_srel = [1.820,1.789,0.131,1.0], [0.519,1.595,0.04,1.0], [0.519,1.595,0.005,1.0]
#
# Now write it all to files for ParFlow
pf.set_layer_depth(d_soil, d_wshale)
pf.set_layer_K(K_soil, K_wshale, K_fshale)
pf.set_layer_por(por_soil, por_wshale, por_fshale)
pf.set_vg_perm(soil_krel, wshale_krel, fshale_krel)
pf.set_vg_sat(soil_srel, whsale_srel, fshale_srel)
#
pf.build_mats_df()
#
pf.apply_exp_K(Ktop=K_fshale, Kbot=K_fshale) # no decay
#pf.apply_exp_K(Ktop=K_fshale, Kbot=K_fshale*1.e-3)  # decay by 3 orders of magnitude
#
pf.make_grid(NX, '{:05d}'.format(savenum))
"""


