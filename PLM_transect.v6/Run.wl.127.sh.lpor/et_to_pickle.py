# update 10/13/22

# Script to pull simulated outputs from ParFlow*.pfb files and CLM*.pfb files


import numpy as np
import pandas as pd
import pickle as pk
import os
import glob
from parflowio.pyParflowio import PFData




def read_pfb(fname):
    pfdata = PFData(fname)
    pfdata.loadHeader()
    pfdata.loadData()
    return pfdata.copyDataArray()



#
# ET output files
#
#nsteps = 1793 + 1

def pull_pf_et(dirname, filename, nsteps):
    et_out_dict = {}
    for i in range(nsteps):
        # Read in the file
        nn = '{0:05d}'.format(i)
        et_f    = './{}/{}.out.evaptrans.{}.pfb'.format(dirname, filename, nn)
        if os.path.exists(et_f):
            et_out  = read_pfb(et_f)[:,0,:]
            et_out_dict[i] = et_out
        else:
            print ('timstep {} does not exists'.format(i))
    # Save it
    savename = '_'.join(dirname.split('_')[1:])
    with open('parflow_out/et_out_dict.{}.pk'.format(savename), 'wb') as f:
        pk.dump(et_out_dict, f)        



#
# CLM output files
#

# pg. 147 of parflow manual
clm_var_base = ['eflx_lh_tot',  'eflx_lwrad_out','eflx_sh_tot',   'eflx_soil_grnd',
                'qflx_evap_tot','qflx_evap_grnd','qflx_evap_soil','qflx_evap_veg',
                'qflx_tran_veg','qflx_infl',     'swe_out',       't_grnd',        't_soil']

# have clm set to 4 soil layers
clm_var = clm_var_base + ['t_soil_{}'.format(i) for i in range(4)]


def pull_clm(dirname, filename, nsteps):
    clm_out_dict = {}
    for i in range(nsteps):
        # Read in the file
        nn = '{0:05d}'.format(i)
        clm_f    = './{}/{}.out.clm_output.{}.C.pfb'.format(dirname, filename, nn)
        if os.path.exists(clm_f):
            clm_out  = read_pfb(clm_f)
            clm_df   = pd.DataFrame(clm_out[:,0,:].T, columns=clm_var)
            clm_out_dict[i] = clm_df
        else:
            print ('timstep {} does not exists'.format(i))
    # Save it
    savename = '_'.join(dirname.split('_')[1:])
    with open('parflow_out/clm_out_dict.{}.pk'.format(savename), 'wb') as f:
        pk.dump(clm_out_dict, f)      



def find_nsteps(dirname):
    l = glob.glob(os.path.join('.', dirname,'*evaptrans*'))
    nl = []
    for ll in l:
        try:
             nl.append(int(ll.split('.')[-2]))
        except ValueError:
            pass
    nsteps = len(nl)+1  
    return nsteps


#
# WY 2000 - 2016
print ('Working on WY 2000 to 2016 ET files')
dirname, filename = 'wy_2000_2016', 'wy_2000_2016'
nsteps = find_nsteps(dirname)
pull_clm(dirname, filename, nsteps)

#
# WY 2017 - 2021
print ('Working on wy 2017 to 2021 ET files')
dirname, filename = 'wy_2017_2021', 'wy_2017_2021'
nsteps = find_nsteps(dirname)
pull_clm(dirname, filename, nsteps)













