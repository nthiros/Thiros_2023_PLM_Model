

import pandas as pd
import numpy as np
import pickle
import os
import sys
import glob
from scipy import stats

from parflowio.pyParflowio import PFData
import pyvista as pv
import parflow.tools as pftools


import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator
plt.rcParams['font.size']=14
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.dates as mdates


import pdb


#------------------------------------
#
# Field Observations -- Water Levels
#
#------------------------------------
wells = pd.read_csv('../PLM_transect.v6/ER_PLM_ParFlow/utils/wells_2_pf_v4.dummy.csv', index_col='well') 


# Transducer Data
# From ESS-dive: Faybishenko 2022,  https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1866836
plm1_obs = pd.read_csv('../PLM_transect.v6/ER_PLM_ParFlow/Field_Data/er_plm1_waterlevel_daily_2016-2022.csv', skiprows=5, index_col='Date_time', parse_dates=True)
plm6_obs = pd.read_csv('../PLM_transect.v6/ER_PLM_ParFlow/Field_Data/er_plm6_waterlevel_daily_2016-2022.csv', skiprows=5, index_col='Date_time', parse_dates=True)

# Water Table Depth below land surface
plm1_obs.loc[:,'bls'] = 2787.16 - plm1_obs.copy().loc[:,'Waterlevel_1d_trend']
plm6_obs.loc[:,'bls'] = 2759.65 - plm6_obs.copy().loc[:,'Waterlevel_1d_trend']

plm1_obs.index = plm1_obs.index.normalize()
plm6_obs.index = plm6_obs.index.normalize()

plm1_obs = plm1_obs[np.where((plm1_obs.index>'2016-09-30') & (plm1_obs.index<'2021-10-01'), True, False)]
plm6_obs = plm6_obs[np.where((plm6_obs.index>'2016-09-30') & (plm6_obs.index<'2021-10-01'), True, False)]

# Clean plm6 offset 
#plm6_obs.iloc[:282,:] = plm6_obs.iloc[:282,:] + plm6_obs[plm6_obs.index < '2017-10-01'].diff(periods=2).max()
#_diff = plm6_obs.loc['2017-09-06','bls'] - plm6_obs.loc['2017-09-08','bls']
_diff = plm6_obs.loc['2017-09-03','bls'] - plm6_obs.loc['2017-09-08','bls']
plm6_obs.loc[plm6_obs.index < '2017-09-03','bls'] = plm6_obs.loc[plm6_obs.index < '2017-09-03','bls'] - _diff
msk1 = np.where((plm6_obs.index>='2017-09-03')&(plm6_obs.index<'2017-09-08'), True, False)
plm6_obs.loc[msk1,'bls'] = np.nan
plm6_obs['bls'].interpolate(inplace=True)







#------------------------------------------
#
# Forcing MET data
#
#------------------------------------------
# 2012-2021
met       = pd.read_csv('./MET/met.2012-2021.1hr.txt', delim_whitespace=True, names=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])
tstart    = pd.to_datetime('2012-10-01 00', format='%Y-%m-%d %H')
tend      = pd.to_datetime('2021-09-30 23', format='%Y-%m-%d %H') # Water year 21 is not over yet
hours     = pd.DatetimeIndex(pd.Series(pd.date_range(tstart, tend, freq='1H')))
met.index = hours





#
# Utility Funcitons
#
def pf_2_dates(startdate, enddate, f, dropleap=False):
    '''Assumes ParFlow outputs every 24 hours'''
    s = pd.to_datetime(startdate)
    e = pd.to_datetime(enddate)
    d_list = pd.date_range(start=s, end=e, freq=f)
    if dropleap:
        d_list = d_list[~((d_list.month == 2) & (d_list.day == 29))]
    return d_list

dates    = pf_2_dates('2012-10-01', '2021-09-30', '24H')
date_map = pd.DataFrame(dates, columns=['Date'])


def set_wy(df):
    #pdb.set_trace()
    try:
        dates = pd.Index(df.copy()['Date'])
    except KeyError:
        dates = pd.Index(df.copy().index)
    yrs       = dates.year
    yrs_      = np.unique(yrs)[1:]
    wy_inds_  = [np.where((dates > '{}-09-30'.format(i-1)) & (dates < '{}-10-01'.format(i)), True, False) for i in yrs_]
    wy_inds   = np.array([wy_inds_[i]*yrs_[i] for i in range(len(yrs_))]).sum(axis=0)
    first_yrs = [(wy_inds==i).argmax() for i in yrs_]
    return list(wy_inds), list(first_yrs)



# Create a new directory for figures
if os.path.exists('figures') and os.path.isdir('figures'):
    pass
else:
    os.makedirs('figures')




#-------------------------------------------------------
#
# Water Level Plots
#
#-------------------------------------------------------


wells  = {'PLM1':404, 'PLM6':494}
w  = list(wells.keys())
w_ = list(wells.values())


#
# WL 2017-2021 Plot
#
#dirlist = glob.glob('/Volumes/Extreme SSD/ParFlow_Runs/*Run0*')
dirlist = glob.glob('../wl_montecarlo/*Run0*')
labs    = [int(i.split('/')[-1][3:]) for i in dirlist]
dmap    = dict(zip(labs,dirlist))


cm = plt.cm.tab20(np.linspace(0,1,len(dirlist)))
#cm = plt.cm.tab10(np.linspace(0,1,10))


#
# Using saturation field and Parflow Hydrylogy Utils
#
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9.0, 2.2))
fig.subplots_adjust(hspace=0.07, wspace=0.3, top=0.88, right=0.95, bottom=0.12, left=0.14)
# Plot water levels
for i in range(len(w)):
    for j in range(len(dirlist)):
        try:
            # Using saturation field and Parflow Hydrylogy Utils
            wll = pd.read_pickle('{}/parflow_out/pf_out_dict_1221.pk'.format(dirlist[j]))['wtd']
            ax[i].plot(dates, wll[:,w_[i]], color=cm[j], alpha=0.5)#, label=labs[j])
        except ValueError:
            print ('value:', dirlist[j])
            pass
        except IndexError:
            print ('index:', dirlist[j])
            pass
        except FileNotFoundError:
            print ('file:', dirlist[j])
            pass
    ax[i].invert_yaxis()
    # Now add field observations
    if w[i] == 'PLM1':
        ax[i].plot(plm1_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    elif w[i] == 'PLM6':
        ax[i].plot(plm6_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    # Fill between years
    for ii in np.arange(2017,2021).reshape(2,2):
        ax[i].axvspan(pd.to_datetime('{}-10-01'.format(ii[0])), pd.to_datetime('{}-09-30'.format(ii[1])), alpha=0.04, color='red')
    ax[i].yaxis.set_major_locator(MaxNLocator(4))
    # Xticks as october 1st
    ax[i].xaxis.set_major_locator(mdates.YearLocator(base=1, month=10, day=1))
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    #ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7]))
    ax[i].xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    #ax[i].yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[i].set_ylabel(w[i])
    ax[i].tick_params(axis='y', which='both', right=False)
    ax[i].margins(x=0.01)
    ax[i].grid()
    ax[i].set_xlim(pd.to_datetime('2016-09-30'), pd.to_datetime('2021-10-02'))
# Clean up
fig.text(0.53, 0.95, 'All 128 Monte Carlo Runs', va='center', ha='center')
fig.text(0.01, 0.55, 'Water Table\nDepth (mbls)', va='center', rotation='vertical')
#ax[0].tick_params(axis='x', labelbottom=False)
#ax[1].tick_params(axis='x', labelbottom=True, rotation=0)#, pad=0.1)
#ax[0].legend(loc='upper left', bbox_to_anchor=(0.98, 1.08), handlelength=0.9, labelspacing=0.5, handletextpad=0.25)
ax[0].legend(loc='lower right',  handlelength=0.9, labelspacing=0.5, handletextpad=0.25)
plt.savefig('./figures/waterlevels_ens.jpg', dpi=300)
plt.savefig('./figures/waterlevels_ens.svg', format='svg')
plt.show()





# 
# Find Best Models
#
plm1_sse  = pd.DataFrame()  # use sum-squared error metric
plm6_sse  = pd.DataFrame()
plm_rmse  = pd.DataFrame()


date_mask_plm1 = np.array([i in plm1_obs.index for i in dates])
date_mask_plm6 = np.array([i in plm6_obs.index for i in dates])
  

for j in range(len(dirlist)):
    try:
        plm1_nse_den = ((plm1_obs['bls'].to_numpy()-plm1_obs['bls'].to_numpy().mean())**2).sum()
        plm6_nse_den = ((plm6_obs['bls'].to_numpy()-plm6_obs['bls'].to_numpy().mean())**2).sum()    
    
        # Using saturation field and Parflow Hydrylogy Utils
        mod = int(dirlist[j].split('/')[-1][3:])
        wll = pd.read_pickle('{}/parflow_out/pf_out_dict_1221.pk'.format(dirlist[j]))['wtd'] # modeled
        plm1_mod = wll[date_mask_plm1, wells['PLM1']]
        plm6_mod = wll[date_mask_plm6, wells['PLM6']]
        
        plm1_res = np.array(plm1_mod - plm1_obs['bls'])
        plm6_res = np.array(plm6_mod - plm6_obs['bls'])
        
        plm1_sse.loc[mod,'sse']  = (plm1_res**2).sum()
        plm6_sse.loc[mod,'sse']  = (plm6_res**2).sum()
        plm_rmse.loc[mod,'rmse'] = np.sqrt((plm1_res**2 + plm6_res**2).sum())
        
        #plm1_nse.append(1 - (((plm1_mod - plm1_obs['bls'])**2).sum() / plm1_nse_den)) 
        #plm6_nse.append(1 - (((plm6_mod - plm6_obs['bls'])**2).sum() / plm6_nse_den)) 
        
        # Using pressure at screen depths
        #wll = pd.read_csv('{}/parflow_out/wy_2012_2021_wt_bls.csv'.format(dirlist[j]))[1:]
        #plm1_mod = wll[date_mask_plm1]['PLM1']
        #plm6_mod = wll[date_mask_plm6]['PLM6']
        #plm1_sse.append(((plm1_mod.to_numpy() - plm1_obs['bls'].to_numpy())**2).sum())
        #plm6_sse.append(((plm6_mod.to_numpy() - plm6_obs['bls'].to_numpy())**2).sum())
        #plm1_nse.append(1 - (((plm1_mod.to_numpy() - plm1_obs['bls'].to_numpy())**2).sum() / plm1_nse_den)) 
        #plm6_nse.append(1 - (((plm6_mod.to_numpy() - plm6_obs['bls'].to_numpy())**2).sum() / plm6_nse_den)) 
        
    except ValueError:
        print ('value:', dirlist[j])
        plm1_sse.loc[mod, 'sse'] = np.NaN
        plm6_sse.loc[mod, 'sse'] = np.NaN
        plm_rmse.loc[mod, 'rmse'] = np.NaN
        #plm1_sse.append(np.NaN)
        #plm6_sse.append(np.NaN)
        #plm_rmse.append(np.NaN)
        pass
    except IndexError:
        print ('index:', dirlist[j])
        plm1_sse.loc[mod, 'sse'] = np.NaN
        plm6_sse.loc[mod, 'sse'] = np.NaN
        plm_rmse.loc[mod, 'rmse'] = np.NaN
        #plm1_sse.append(np.NaN)
        #plm6_sse.append(np.NaN)
        #plm_rmse.append(np.NaN)
        pass
    except FileNotFoundError:
        print ('file:', dirlist[j])
        plm1_sse.loc[mod, 'sse'] = np.NaN
        plm6_sse.loc[mod, 'sse'] = np.NaN
        plm_rmse.loc[mod, 'rmse'] = np.NaN
        #plm1_sse.append(np.NaN)
        #plm6_sse.append(np.NaN)
        #plm_rmse.append(np.NaN)
        pass

plm1_opt_ = plm1_sse.sort_values(by='sse')
plm6_opt_ = plm6_sse.sort_values(by='sse')

# consider both wells together
plm_opt_  = (plm1_sse + plm6_sse).sort_values(by='sse')
plm_rmse_ = plm_rmse.sort_values(by='rmse')


# Just the indices
plm1_opt = plm1_opt_.index.to_numpy()
plm6_opt = plm6_opt_.index.to_numpy()
plm_opt  = plm_opt_.index.to_numpy()

plm_rmse = plm_rmse_.index.to_numpy()


#
# --- Plotting ---
#
nmods = 6 # plot the top 4 models
cm = plt.cm.coolwarm(np.linspace(0.0,1.0,nmods))

# --- Plot 1 --- best models according to plm1 only
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9.0, 2.2))
fig.subplots_adjust(hspace=0.07, wspace=0.3, top=0.88, right=0.95, bottom=0.12, left=0.14)
# Plot water levels
for i in range(len(w)):
    for j in range(nmods):  # 10 best models
        try:
            # Using saturation field and Parflow Hydrylogy Utils
            wll = pd.read_pickle('{}/parflow_out/pf_out_dict_1221.pk'.format(dmap[plm1_opt[j]]))['wtd'] 
            ax[i].plot(dates, wll[:,w_[i]], color=cm[j], alpha=0.9, label=labs[j])
            #wll = pd.read_csv('{}/parflow_out/wy_2012_2021_wt_bls.csv'.format(dmap[plm1_opt[j]]))[1:]             
            #ax[i].plot(dates, wll.loc[:,w[i]], color=cm[j], alpha=0.9, label=labs[j])
        except ValueError:
            print ('value:', dirlist[j])
            pass
        except IndexError:
            print ('index:', dirlist[j])
            pass
        except FileNotFoundError:
            print ('file:', dirlist[j])
            pass
    ax[i].invert_yaxis()
    # Now add field observations
    if w[i] == 'PLM1':
        ax[i].plot(plm1_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    elif w[i] == 'PLM6':
        ax[i].plot(plm6_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    # Fill between years
    for ii in np.arange(2017,2021).reshape(2,2):
        ax[i].axvspan(pd.to_datetime('{}-10-01'.format(ii[0])), pd.to_datetime('{}-09-30'.format(ii[1])), alpha=0.04, color='red')
    ax[i].yaxis.set_major_locator(MaxNLocator(4))
    # Xticks as october 1st
    ax[i].xaxis.set_major_locator(mdates.YearLocator(base=1, month=10, day=1))
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    #ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7]))
    ax[i].xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    #ax[i].yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[i].set_ylabel(w[i])
    ax[i].tick_params(axis='y', which='both', right=False)
    ax[i].margins(x=0.01)
    ax[i].grid()
    ax[i].set_xlim(pd.to_datetime('2016-09-30'), pd.to_datetime('2021-10-02'))
# Clean up
fig.text(0.53, 0.95, 'Best Models Considering PLM1 Only', va='center', ha='center')
fig.text(0.01, 0.55, 'Water Table\nDepth (mbls)', va='center', rotation='vertical')
#ax[0].legend(loc='upper left', bbox_to_anchor=(0.98, 1.08), handlelength=0.9, labelspacing=0.5, handletextpad=0.25)
plt.savefig('./figures/waterlevels_ens.plm1.jpg', dpi=300)
plt.savefig('./figures/waterlevels_ens.plm1.svg', format='svg')
plt.show()



# --- Plot 2 --- best models according to plm6 only
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9.0, 2.2))
fig.subplots_adjust(hspace=0.07, wspace=0.3, top=0.88, right=0.95, bottom=0.12, left=0.14)
# Plot water levels
for i in range(len(w)):
    for j in range(nmods):  # 10 best models
        try:
            # Using saturation field and Parflow Hydrylogy Utils
            wll = pd.read_pickle('{}/parflow_out/pf_out_dict_1221.pk'.format(dmap[plm6_opt[j]]))['wtd']              
            ax[i].plot(dates, wll[:,w_[i]], color=cm[j], alpha=0.9, label=labs[j])
            #wll = pd.read_csv('{}/parflow_out/wy_2012_2021_wt_bls.csv'.format(dmap[plm6_opt[j]]))[1:]             
            #ax[i].plot(dates, wll.loc[:,w[i]], color=cm[j], alpha=0.9, label=labs[j])
        except ValueError:
            print ('value:', dirlist[j])
            pass
        except IndexError:
            print ('index:', dirlist[j])
            pass
        except FileNotFoundError:
            print ('file:', dirlist[j])
            pass
    ax[i].invert_yaxis()
    # Now add field observations
    if w[i] == 'PLM1':
        ax[i].plot(plm1_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    elif w[i] == 'PLM6':
        ax[i].plot(plm6_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    # Fill between years
    for ii in np.arange(2017,2021).reshape(2,2):
        ax[i].axvspan(pd.to_datetime('{}-10-01'.format(ii[0])), pd.to_datetime('{}-09-30'.format(ii[1])), alpha=0.04, color='red')
    ax[i].yaxis.set_major_locator(MaxNLocator(4))
    # Xticks as october 1st
    ax[i].xaxis.set_major_locator(mdates.YearLocator(base=1, month=10, day=1))
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    #ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7]))
    ax[i].xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    #ax[i].yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[i].set_ylabel(w[i])
    ax[i].tick_params(axis='y', which='both', right=False)
    ax[i].margins(x=0.01)
    ax[i].grid()
    ax[i].set_xlim(pd.to_datetime('2016-09-30'), pd.to_datetime('2021-10-02'))
#ax[0].set_title('Best Models Considering PLM6 Only', fontsize=14)
# Clean up
fig.text(0.53, 0.95, 'Best Models Considering PLM6 Only', va='center', ha='center')
fig.text(0.01, 0.55, 'Water Table\nDepth (mbls)', va='center', rotation='vertical')
#ax[0].legend(loc='upper left', bbox_to_anchor=(0.98, 1.08), handlelength=0.9, labelspacing=0.5, handletextpad=0.25)
plt.savefig('./figures/waterlevels_ens.plm6.jpg', dpi=300)
plt.savefig('./figures/waterlevels_ens.plm6.svg', format='svg')
plt.show()



# --- Plot 3 --- best models according to plm6 and plm1 together
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9.0, 2.2))
fig.subplots_adjust(hspace=0.07, wspace=0.3, top=0.88, right=0.95, bottom=0.12, left=0.14)
# Plot water levels
for i in range(len(w)):
    for j in range(nmods):  # 10 best models
        try:
            wll =  pd.read_pickle('{}/parflow_out/pf_out_dict_1221.pk'.format(dmap[plm_opt[j]]))['wtd']
            ax[i].plot(dates, wll[:,w_[i]], color=cm[j], alpha=0.9, label='A{}'.format(labs[j]+1))
            #wll = pd.read_csv('{}/parflow_out/wy_2012_2021_wt_bls.csv'.format(dmap[plm_opt[j]]))[1:]             
            #ax[i].plot(dates, wll.loc[:,w[i]], color=cm[j], alpha=0.9, label=labs[j])
        except ValueError:
            print ('value:', dirlist[j])
            pass
        except IndexError:
            print ('index:', dirlist[j])
            pass
        except FileNotFoundError:
            print ('file:', dirlist[j])
            pass
    ax[i].invert_yaxis()
    # Now add field observations
    if w[i] == 'PLM1':
        ax[i].plot(plm1_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    elif w[i] == 'PLM6':
        ax[i].plot(plm6_obs['bls'], color='black', alpha=1.0, linewidth=1.5, linestyle='--', label='Obs.')
    # Fill between years
    for ii in np.arange(2017,2021).reshape(2,2):
        ax[i].axvspan(pd.to_datetime('{}-10-01'.format(ii[0])), pd.to_datetime('{}-09-30'.format(ii[1])), alpha=0.04, color='red')
    ax[i].yaxis.set_major_locator(MaxNLocator(4))
    # Xticks as october 1st
    ax[i].xaxis.set_major_locator(mdates.YearLocator(base=1, month=10, day=1))
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    #ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7]))
    ax[i].xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    #ax[i].yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[i].set_ylabel(w[i])
    ax[i].tick_params(axis='y', which='both', right=False)
    ax[i].margins(x=0.01)
    ax[i].grid()
    ax[i].set_xlim(pd.to_datetime('2016-09-30'), pd.to_datetime('2021-10-02'))
#ax[0].set_title('Best Models Considering PLM1 and PLM6', fontsize=14)
# Clean up
fig.text(0.53, 0.95, 'Best Models Considering PLM1 and PLM6 Simultaneously', va='center', ha='center')
fig.text(0.01, 0.55, 'Water Table\nDepth (mbls)', va='center', rotation='vertical')
#ax[0].legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), handlelength=0.9, labelspacing=0.5, handletextpad=0.25)
plt.savefig('./figures/waterlevels_ens.both.jpg', dpi=300)
plt.savefig('./figures/waterlevels_ens.both.svg', format='svg')
plt.show()






# --- Table With Parameter Values ---
#samps = pd.read_csv('/Volumes/Extreme SSD/ParFlow_Runs/Training_Set/sobol_parameters.csv', index_col='Run')
samps = pd.read_csv('../wl_montecarlo/Training_Set/sobol_parameters.csv', index_col='Run')
#samps.loc[:, ['K_soil','K_wshale','K_fshale']] = 10**(samps.loc[:, ['K_soil','K_wshale','K_fshale']])

samps_prior = pd.DataFrame([samps.loc[0,:]])
samps_prior.insert(loc=0, column='well',value='prior')

samps_plm1 = samps.loc[plm1_opt[:4],:]
samps_plm1.insert(loc=0, column='well',value='plm1')

samps_plm6 = samps.loc[plm6_opt[:4],:]
samps_plm6.insert(loc=0, column='well',value='plm6')

samps_both = samps.loc[plm_opt[:5],:]
samps_both.insert(loc=0, column='well',value='both')


samps_opt  = pd.concat((samps_prior, samps_plm1, samps_plm6, samps_both))
samps_opt_ = samps_opt.copy()[['well','K_soil','K_wshale','K_fshale','soil_rel_alpha','wshale_rel_alpha','soil_rel_n','wshale_rel_n']]
#samps_opt_.to_csv('opt_models_wl.csv', index=True)

samps_all  = samps.copy().loc[plm_rmse]
samps_all_ = samps_all.copy()[['K_soil','K_wshale','K_fshale','soil_rel_alpha','wshale_rel_alpha','soil_rel_n','wshale_rel_n']]
samps_all_['rmse'] = plm_rmse_['rmse'].to_numpy()



# --- Bi-Variate distributions with Parameter Values ---

import itertools


_labmap = {'K_soil' : r'$log_{10} K_{soil}$',
           'K_wshale' : r'$log_{10} K_{wshale}$',
           'K_fshale' : r'$log_{10} K_{fshale}$',
           'soil_rel_alpha' : r'$\alpha_{soil}$',
           'wshale_rel_alpha' : r'$\alpha_{wshale}$',
           'soil_rel_n' : r'$n_{soil}$',
           'wshale_rel_n' : r'$n_{wshale}$'}



samps_both = samps.loc[plm_opt[:5],:]
samps_both = samps_both.copy()[['K_soil','K_wshale','K_fshale','soil_rel_alpha','wshale_rel_alpha','soil_rel_n','wshale_rel_n']]

base = samps_both.iloc[0,:].copy()
#base['well']
base['K_soil']   = np.log10(2.5e-5)
base['K_wshale'] = np.log10(4.0e-6)
base['K_fshale'] = np.log10(2.2e-8)
base['soil_rel_alpha'] = 1.82
base['wshale_rel_alpha'] = 0.52
base['soil_rel_n'] = 1.79
base['wshale_rel_n'] = 1.6
base.name = 'base'
base = pd.DataFrame(base).T



# top models and base
samps_opt  = pd.concat((base, samps_both))
samps_opt_ = samps_opt.copy()[['K_soil','K_wshale','K_fshale','soil_rel_alpha','wshale_rel_alpha','soil_rel_n','wshale_rel_n']]
samps_opt_.index.name = 'model'
samps_opt_.to_csv('opt_models_wl.csv', index=True)




subplot_vals = np.array(list(itertools.combinations(np.arange(len(samps_both.columns)), 2)))
#subplot_vals[:,1] = np.subtract(subplot_vals[:,1], subplot_vals[:,0]+1)

pars = list(itertools.combinations(samps_both.columns, 2))


cm = plt.cm.coolwarm(np.linspace(0.0,1.0, len(samps_both)))


fig, axes = plt.subplots(ncols=len(samps_both.columns), nrows=len(samps_both.columns), figsize=(10,8))
fig.subplots_adjust(top=0.9, bottom=-0.1, right=0.9, left=-0.12, hspace=0.1, wspace=0.1)
for i in range(len(samps_both.index)):
    for ii in range(len(subplot_vals)):
        r  = subplot_vals[ii][0]
        c  = subplot_vals[ii][1]
        vr = pars[ii][1]
        vc = pars[ii][0]
        
        ax = axes[r,c]
        _i = samps_both.index[i]
        #ax.scatter(samps_both.loc[_i,vr], samps_both.loc[_i,vc], color=cm[i], label='A{}'.format(i+1))
        ax.scatter(samps_both.loc[_i,vr], samps_both.loc[_i,vc], color=cm[i], s=50, label='A{} ({:.1f})'.format(i+1, plm_rmse_.iloc[i,0]))
        
        ax.minorticks_on()
        ax.tick_params(axis='y', which='both', right=True, labelright=True, labelleft=False)
        ax.tick_params(axis='x', which='both', top=True, labeltop=True, labelbottom=False, rotation=0, pad=0.001)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
        ax.xaxis.set_major_locator(ticker.MaxNLocator(3))
        if r == 0:
            ax.set_xlabel(_labmap[vr])
            ax.xaxis.set_label_position("top")
        else: 
            ax.tick_params(axis='x', top=True, labeltop=False, labelbottom=False)
            
        if c == len(samps_opt_.columns[1:])-1:
            ax.set_ylabel(_labmap[vc])
            ax.yaxis.set_label_position("right")
        else:
            ax.tick_params(axis='y', right=True, labelright=False, labelleft=False)    
        
        ax.set_xlim(samps[vr].min(), samps[vr].max())
        ax.set_ylim(samps[vc].min(), samps[vc].max())
        
        if i == 0:
            # parameter bounds from MC sampling
            #ax.axhline(samps[vc].min(), linestyle=':', color='grey', alpha=0.9)
            #ax.axhline(samps[vc].max(), linestyle=':', color='grey', alpha=0.9)
            #ax.axvline(samps[vr].min(), linestyle=':', color='grey', alpha=0.9)
            #ax.axvline(samps[vr].max(), linestyle=':', color='grey', alpha=0.9)
    
            # best-fit line
            m, intr, r_value, p_value, std_err = stats.linregress(samps_both[vr].to_numpy(), samps_both[vc].to_numpy())
            xx = np.linspace(samps[vr].min(), samps[vr].max(), 100)
            yy = m*xx + intr
            ax.plot(xx, yy, linestyle='--', color='black', alpha=0.6)
            if (r,c) == (0,1):
                ax.text(0.4, 0.09, r'R$^{{2}}=${:.1f}'.format(r_value**2), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            else:
                ax.text(0.85, 0.09, '{:.1f}'.format(r_value**2), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        if i == len(samps_both.index)-1:
            # plot the base
            ax.scatter(base.loc['base',vr], base.loc['base',vc], color='black', marker='s', label='base ({:.1f})'.format(plm_rmse_.loc[0,'rmse']))
# turn off all subplots that are not being used
_negsubplots = np.array(list(itertools.product(np.arange(len(samps_both.columns)), repeat=2)))
for zz in _negsubplots:
    r,c = zz
    if  (zz==subplot_vals).all(axis=1).sum()==1:
        pass
    else:
        try:
            axes[r,c].axis('off')
        except IndexError:
            pass
axes[0,1].legend(ncol=1, loc='upper left', bbox_to_anchor=(-0.05, 0.0), handlelength=0.9, labelspacing=0.5, columnspacing=0.8, handletextpad=0.2, frameon=False)
plt.savefig('./figures/parameter_cov.jpg', dpi=300)
plt.savefig('./figures/parameter_cov.svg', format='svg')
plt.show()
    



#
# --- Parameter values versus RMSE for top (all?) models
#
fig, axes = plt.subplots(ncols=len(samps_all_.columns)-1, nrows=1, figsize=(10,2))
fig.subplots_adjust(wspace=0.14, left=0.1, right=0.98, bottom=0.3)
for i in range(len(samps_all_.columns)-1):
    ax = axes[i]
    _p = samps_all_.iloc[:,i].to_numpy()
    _rmse = samps_all_['rmse'].to_numpy()
    ax.scatter(_p, _rmse, marker='.', color='black')
    ax.set_yscale('log')
    ax.set_xlabel(_labmap[samps_all_.columns[i]], labelpad=0.1)
    ax.minorticks_on()
    ax.tick_params(axis='x', rotation=25, pad=0.1)
    if i != 0:
        ax.tick_params(axis='y', labelleft=False)
    else:
        ax.set_ylabel('Water Level\nRMSE (m)')
    # replot top with color
    for j in range(5):
        ax.scatter(_p[j], _rmse[j], marker='.', s=90, color=cm[j])
plt.savefig('./figures/parameter_rmse.jpg', dpi=300)
plt.savefig('./figures/parameter_rmse.svg', format='svg')
plt.show()









"""
#
# --- ET disributions ---
#
def read_pf_out(_dir, _var):
    #pdb.set_trace()
    clm_out_ = pd.read_pickle(os.path.join(_dir,'parflow_out','clm_out_dict.2012_2021.pk')) # contains ET from CLM output files
    clm_out  = {i: v for i, v in enumerate(clm_out_.values())}  # Reset timestep keys to start with zero index
    clm_keys = list(clm_out[0].keys())
    
    df_    = pd.DataFrame((clm_out[i][_var]) for i in list(clm_out.keys())).T
    df_.columns = dates
    df = pd.DataFrame(df_)
    return df


_clm_swe  = []
_clm_et   = []

xpos = 559 - 10
dx, dy = 1.5125*1000, 1.5125*1000

# Put everyting in m3/day
for j in range(len(dirlist)):
    try:
        # SWE from CLM -- mm
        #swe = read_pf_out(dirlist[j],'swe_out').sum(axis=0) # mm of SWE on the hillslope per day
        swe = (read_pf_out(dirlist[j],'swe_out')*dx*dy).sum(axis=0) / 1000**3 # m3/day
        _clm_swe.append(swe)
    
        # ET from CLM -- mm/s
        #clm_et = pd.DataFrame((read_pf_out(dirlist[j],'qflx_evap_tot').loc[:xpos,:]*dx*dy).sum(axis=0)*3600*24)/1000**3 # m3/day
        clm_et = pd.DataFrame((read_pf_out(dirlist[j],'qflx_evap_tot')*dx*dy).sum(axis=0)*3600*24)/1000**3 # m3/day
        clm_et.rename(columns={0:'ET'}, inplace=True)
        clm_et['wy'] = set_wy(clm_et)[0]
        #clm_et['ET'][clm_et['ET']<0.0] = 0.0
        clm_et_cs = clm_et.groupby(by='wy').cumsum()
        _clm_et.append(clm_et_cs)
    

    except ValueError:
        print ('value:', dirlist[j])
        pass
    except KeyError:
        print ('key:', dirlist[j])
        pass
    except IndexError:
        print ('index:', dirlist[j])
        pass
    except FileNotFoundError:
        print ('file:', dirlist[j])
        pass



fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8.0,2.7))
fig.subplots_adjust(hspace=0.15, top=0.92, bottom=0.1, right=0.78, left=0.15)
# Total Precip from NLDAS
dx, dy = 1.5125*1000, 1.5125*1000
#ax.plot(met_cs*dx*dy*559/(1000**3), color='black', alpha=0.8, linestyle='--', label='Precip.')
#ax.fill_between(x=met_cs.index, y2=0, y1=met_cs.to_numpy().ravel()*dx*dy*559/(1000**3), color='grey', alpha=0.1)  
ax.fill_between(x=swe.index, y2=np.array(_clm_swe).min(axis=0),          y1=np.array(_clm_swe).max(axis=0), color='C0', alpha=0.65, label='CLM SWE')
ax.fill_between(x=swe.index, y2=np.array(_clm_et)[:,:,0].min(axis=0),    y1=np.array(_clm_et)[:,:,0].max(axis=0), color='C1', alpha=0.65, label='CLM ET')
ax.yaxis.set_major_locator(ticker.MultipleLocator(200))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.set_ylabel('Annual Cumulative\nFlux (m$^{3}$/yr)', labelpad=0.5)
ax.grid()
ax.margins(x=0.01)
ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
for z in np.arange(2017,2021).reshape(2,2):
    ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.04, color='red')
ax.legend(ncol=1, loc='upper center', bbox_to_anchor=(1.13, 1.03), handlelength=1.0, labelspacing=0.25, handletextpad=0.5, columnspacing=0.5, fontsize=13)
#plt.savefig('./figures/clm_comps.png', dpi=300)
plt.show()

"""



