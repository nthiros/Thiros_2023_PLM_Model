import numpy as np
import pandas as pd
import sys
import subprocess
import pdb

import skopt
from scipy.stats import qmc

import make_cake_utils as d_utils



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
# 							 # Depth below land surface
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
# Sobol Sampling of parameters


# Define prior means for parameters
par_priors =   {'dz_soil':    1.0,                 # tickness of soil layer in meters
                'dz_wshale':  3.0,                 # thicknes of w.shale layer in meters
                
                'K_soil':     np.log10(1.0e-04),    # log10 K (m/s)
                
                #'K_wshale':  np.log10(5.605e-06),   # log10 K
                'K_wshale_f': -1.0,                  # log-factor decrease of wshale compared to soil K; 0 means log(K_wshale)=log(K_soil), -1 means log(K_wshale)=log(K_soil)-1
                                                     # This is a way to enforce K_soil >= K_wshale >= K_fshale
                                                     
                #'K_fshale':  np.log10(1.0e-08),     # log10 K
                'K_fshale_f': -1.0,                  # log-factor decrease on K_wshale for top of K_fshale
                'K_exp':       0.0,                  # log factor to decrease fshale w/ exp decay
                                                     # K_fshale_f=-2 means log(K_fshale)=log(K_wshale)-2; 
                                                     # K_exp=0.0 is no decay
                
                'por_soil':    0.4,               # porosity
                'por_wshale':  0.2,               # porosity
                'por_fshale':  0.1 ,              # porosity
                
                'soil_rel_alpha':   1.820,        # vg alpha
                'wshale_rel_alpha': 1.092,        # vg alpha
                'fshale_rel_alpha': 0.519,        # vg alpha
                
                'soil_rel_n':       1.789,        # vg n
                'wshale_rel_n':     1.363,        # vg n
                'fshale_rel_n':     1.595,        # vg n
                
                'soil_rel_rsat':    0.131,        # vg residual saturation
                'wshale_rel_rsat':  0.040,        # vg residual saturation
                'fshale_rel_rsat':  0.005,        # vg residual saturation
                }
                

# Define the ranges for parameters 
# Only these parameters will be varied, rest well be set at priors (above)
# Intergers will result in integer discretization
par_bounds =   {#'dz_soil':    [2,  4],       # thickness of soil layer
                #'dz_wshale':  [2,  6],       # thickness of wshale layer
                # thickness of fshale will be the remainder of the domain
                
                'K_soil':     [-5.301, -3.301],  # log10 K (m/s)
               
                #'K_wshale':   [-7.0, -5.0],     # log10 K -- using factor decreases (below) instead of absolute values
                'K_wshale_f': [-2.0, 0.0],       # factor decrease (log-space) of K_wshale from K_soil. ie. log_kwshale=log_ksoil-1.0
                
                #'K_fshale':   [-9.0, -7.0],     # log10 K -- using factor decreases (below) instead of absolute values
                'K_fshale_f': [-3.0,  -1.0],      # factor decrease (log-space) of K_fshale from K_wshale. ie. log_fwshale=log_wshale-0.0 (same)
                #'K_exp':      [-2.0, 0.0],      # factor (log-space) to decrease fshale w/ exp decay, 0 is no decay, -2 means bottom is two orders of magnitude smaller K than top
                
                #'por_soil':   [0.3, 0.4],    # porosity
                #'por_wshale': [0.1, 0.3],    # porosity
                #'por_fshale': [0.01, 0.1],   # porosity
                
                'soil_rel_alpha':   [0.1,  5.5],   # vg alpha
                'wshale_rel_alpha': [0.1,  5.5],   # vg alpha
                #'fshale_rel_alpha': [0.05, 0.5],  # vg alpha
                
                'soil_rel_n':       [1.0, 5.0],        # vg n
                'wshale_rel_n':     [1.0, 5.0],        # vg n
                #'fshale_rel_n':     [1.0, 5.0],       # vg n
                
                #'soil_rel_rsat':    0.0,       # vg residual saturation
                #'wshale_rel_rsat':  0.0,       # vg residual saturation
                #'fshale_rel_rsat':  0.0,       # vg residual saturation
                }
                
            


class sample_priors():
    def __init__(self, pars_bnds_pr, num_runs, random_state):
        self.par_bnd = pars_bnds_pr
        self.num_runs = num_runs
        self.keys = list(pars_bnds_pr.keys())
        self.rnd_st = random_state
    
    def random(self):
        '''Simple, random sampling from the bounds'''
        np.random.seed(self.rnd_st)
        arr = np.ones((self.num_runs, len(self.keys)))
        for i in range(len(self.keys)):
            arr[:,i] = np.random.uniform(self.par_bnd[self.keys[i]][0], self.par_bnd[self.keys[i]][1], self.num_runs)
        df = pd.DataFrame(data=arr, columns=self.keys)
        return df
        
    def LHS(self):
        '''Latin Hypercube Sampling'''
        #bnds_arr = pd.DataFrame((self.par_bnd)).to_numpy().T
        bnds_arr = skopt.Space(list(par_bounds.values()))
        space = skopt.Space(bnds_arr)
        lhs_arr = np.array(skopt.sampler.Lhs(lhs_type="classic").generate(space.dimensions, self.num_runs, self.rnd_st))
        lhs_df = pd.DataFrame(data=lhs_arr, columns=self.keys)
        return lhs_df
    
    def sobol(self):
        '''Generate a Sobol Sequence'''
        #bnds_arr = pd.DataFrame((self.par_bnd)).to_numpy().T
        bnds_arr = skopt.Space(list(par_bounds.values()))
        space = skopt.Space(bnds_arr)
        sobol_arr = np.array(skopt.sampler.Sobol().generate(space.dimensions, self.num_runs, self.rnd_st))
        sobol_df = pd.DataFrame(data=sobol_arr, columns=self.keys)
        return sobol_df



num_runs = 64  # the number of training instances to generate - should be a power of two e.g. 2**6
random_state = 22133


lhs_samps   = sample_priors(par_bounds, num_runs, random_state).LHS()
sobol_samps = sample_priors(par_bounds, num_runs, random_state).sobol()

#
# Mangle the parameters for ParFlow
#

# pick LHS or sobol set
#samps = lhs_samps.copy()
samps = sobol_samps.copy()

# Set parameters that were not varied to the prior mean
#samps['soil_rel_rsat']   = par_priors['soil_rel_rsat']
#samps['wshale_rel_rsat'] = par_priors['wshale_rel_rsat']
#samps['fshale_rel_rsat'] = par_priors['fshale_rel_rsat']
#samps['fshale_rel_alpha'] = par_priors['fshale_rel_alpha']
#samps['fshale_rel_n']     = par_priors['fshale_rel_n']
for i in par_priors.keys():
    if i not in samps.columns:
        samps[i] = par_priors[i]



# Tack-on prior means as first row of dataframe - this will be the first run
s_pr = samps.iloc[0,:].copy()
for i in list(par_priors.keys()):
    if i in s_pr.index.to_list():
        s_pr[i] = par_priors[i]
samps = pd.concat((s_pr, pd.DataFrame(samps).T), axis=1).T
samps.reset_index(drop=True, inplace=True)


# Calculate K values based on factor decrase parameters and soil K
#   -- log(K_wshale) = log(K_soil) + factor decrease (the K_wshale_f parameter)
#   -- log(K_fshale) = log(K_wshale) + factor decrease (the K_fshale_f parameter)
# comment out if defining log K directly and not defining relative decreases
samps['K_wshale'] = samps.loc[:,'K_soil']   + samps.loc[:,'K_wshale_f']
samps['K_fshale'] = samps.loc[:,'K_wshale'] + samps.loc[:,'K_fshale_f']


# depths bls to bottom of layer
samps['d_soil']   = (samps.loc[:,'dz_soil'].copy()).astype(int)  
samps['d_wshale'] = (samps.loc[:,'dz_soil'] + samps.loc[:,'dz_wshale']).astype(int)



# Write to files
for n in range(num_runs):
    #----------------------------------------
    # Set hillslope properties
    #
    s = samps.copy()
    savenum = n
    
    #
    # Setup Properties -- Vary these with LHS sampling
    d_soil, d_wshale  = s.loc[n,'d_soil'], s.loc[n,'d_wshale']
    
    K_soil, K_wshale, K_fshale = 10**s.loc[n,'K_soil'], 10**s.loc[n,'K_wshale'], 10**s.loc[n,'K_fshale']   
    
    exp_fact = 10**s.loc[n,'K_exp'] #K_fshale at bottom is this factor smaller than K_fshale at top
    
    por_soil, por_wshale, por_fshale = s.loc[n,'por_soil'], s.loc[n,'por_wshale'], s.loc[n,'por_fshale']   
    
    soil_krel   = [s.loc[n,'soil_rel_alpha'], s.loc[n,'soil_rel_n']]
    soil_srel   =  soil_krel   + [s.loc[n,'soil_rel_rsat']]  + [1.0]
    
    wshale_krel = [s.loc[n,'wshale_rel_alpha'],s.loc[n,'wshale_rel_n']]
    wshale_srel =  wshale_krel + [s.loc[n,'wshale_rel_rsat']] + [1.0] 
    
    fshale_krel = [s.loc[n,'fshale_rel_alpha'],s.loc[n,'fshale_rel_n']]       
    fshale_srel =  fshale_krel + [s.loc[n,'fshale_rel_rsat']] + [1.0]
    
    
    #
    # Now write it all to files for ParFlow
    pf = d_utils.pf_properties(dbs, dz_scale)
    pf.set_layer_depth(d_soil, d_wshale)
    pf.set_layer_K(K_soil, K_wshale, K_fshale)
    pf.set_layer_por(por_soil, por_wshale, por_fshale)
    pf.set_vg_perm(soil_krel, wshale_krel, fshale_krel)
    pf.set_vg_sat(soil_srel, wshale_srel, fshale_srel)
    #
    pf.build_mats_df()
    #
    #.apply_exp_K(Ktop=K_fshale, Kbot=K_fshale) # no decay
    pf.apply_exp_K(Ktop=K_fshale, Kbot=K_fshale*exp_fact)  # decay by 3 orders of magnitude
    #
    pf.make_grid(NX,'{:05d}'.format(savenum))



"""
Single Realiztion
#----------------------------------------
# Set hillslope properties
#
s = samps.copy()
n  = 0 
savenum = n

#
# Setup Properties -- Vary these with LHS sampling
d_soil, d_wshale  = s.loc[n,'d_soil'], s.loc[n,'d_wshale']

K_soil, K_wshale, K_fshale = 10**s.loc[n,'K_soil'], 10**s.loc[n,'K_wshale'], 10**s.loc[n,'K_fshale']   

exp_fact = 10**s.loc[n,'K_exp'] #K_fshale at bottom is this factor smaller than K_fshale at top

por_soil, por_wshale, por_fshale = s.loc[n,'por_soil'], s.loc[n,'por_wshale'], s.loc[n,'por_fshale']   

soil_krel   = [s.loc[n,'soil_rel_alpha'],s.loc[n,'soil_rel_n']]
soil_srel   =  soil_krel   + [s.loc[n,'soil_rel_rsat']]  + [1.0]

wshale_krel = [s.loc[n,'wshale_rel_alpha'],s.loc[n,'wshale_rel_n']]
wshale_srel =  wshale_krel + [s.loc[n,'wshale_rel_rsat']] + [1.0] 

fshale_krel = [s.loc[n,'fshale_rel_alpha'],s.loc[n,'fshale_rel_n']]       
fshale_srel =  fshale_krel + [s.loc[n,'fshale_rel_rsat']] + [1.0]


#
# Now write it all to files for ParFlow
pf = d_utils.pf_properties(dbs, dz_scale)
pf.set_layer_depth(d_soil, d_wshale)
pf.set_layer_K(K_soil, K_wshale, K_fshale)
pf.set_layer_por(por_soil, por_wshale, por_fshale)
pf.set_vg_perm(soil_krel, wshale_krel, fshale_krel)
pf.set_vg_sat(soil_srel, wshale_srel, fshale_srel)
#
pf.build_mats_df()
#
#.apply_exp_K(Ktop=K_fshale, Kbot=K_fshale) # no decay
pf.apply_exp_K(Ktop=K_fshale, Kbot=K_fshale*exp_fact)  # decay by 3 orders of magnitude
#
pf.make_grid(NX, '{:05d}'.format(savenum))
"""













