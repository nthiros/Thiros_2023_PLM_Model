# EcoSLIM needs ParFlow outputs in pfb of:
# velx, vely, velz, saturation, porosity
# out.evaptrans, parflow C.pfb
# 
# NOTE: ecoslim assumes 10 soil layers (can change in Ecoslim.f90)


import os
import numpy as np
import shutil


# directory to pull from
dr = os.path.join(os.getcwd(), 'wy_spinup')
rn = 'wy_spinup' # runname

# new directory to copy into
dr_new = os.path.join(os.getcwd(), 'wy_spinup_lng')

# Total number of output parflow fields, likely set to days=24hr
ft = 3650

# file number to start copying at
fts = ft-365*2  # want to repeat the last year of data
                # using two years because 730/10=73 is and integer value

# pull outputs every tp days, remember original output frequeency is 1 day
tp = 10  # 240 hours

# list of file numbers to use
# note, ts for clm stuff is +1
ts = np.arange(fts, ft+1, tp)
# need to renumber files starting from 0
ts_reset = np.arange(len(ts))


# list of files I want
filename = ['velx','vely','velz','satur','porosity','evaptrans','clm_output']
f  = [] # original file names
fn = [] # updated file names
for i,n in zip(ts, ts_reset):
    for j in filename:
        if j== 'clm_output':
            if i>0:
                ff  = '{}/{}.out.{}.{:05d}.C.pfb'.format(dr, rn, j, i)
                ffn = '{}/{}.out.{}.{:05d}.C.pfb'.format(dr_new, rn, j, n)
        elif j == 'evaptrans':
            if i>0:
                ff  = '{}/{}.out.{}.{:05d}.pfb'.format(dr, rn, j, i)
                ffn = '{}/{}.out.{}.{:05d}.pfb'.format(dr_new, rn, j, n)
        elif j == 'porosity':
            ff  = '{}/{}.out.{}.pfb'.format(dr, rn, j)
            ffn = '{}/{}.out.{}.pfb'.format(dr_new, rn, j)
        else:
            ff  = '{}/{}.out.{}.{:05d}.pfb'.format(dr, rn, j, i)
            ffn = '{}/{}.out.{}.{:05d}.pfb'.format(dr_new, rn, j, n)	
        f.append(ff)
        fn.append(ffn)


# Create a new directory
if os.path.exists(dr_new) and os.path.isdir(dr_new):
    shutil.rmtree(dr_new)
os.makedirs(dr_new)

for i,j in zip(f,fn):
#    fnew = os.path.join(dr_new, f[i].split('/')[-1])
#    #print (fnew) 
    shutil.copyfile(i, j) 




