
import pandas as pd
import numpy as np
import pickle
import os
import sys
import pdb

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
from matplotlib import lines






#------------------------------------
#
# Field Observations -- Water Levels
#
#------------------------------------
wells = pd.read_csv('./ER_PLM_ParFlow/utils/wells_2_pf_v4.dummy.csv', index_col='well') 

# Water level observations from ESS-Dive
plm1_obs = pd.read_csv('./ER_PLM_ParFlow/Field_Data//er_plm1_waterlevel_daily_2016-2022.csv', skiprows=5, index_col='Date_time', parse_dates=True)
plm6_obs = pd.read_csv('./ER_PLM_ParFlow/Field_Data//er_plm6_waterlevel_daily_2016-2022.csv', skiprows=5, index_col='Date_time', parse_dates=True)

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






#-------------------------------------------------------------------
#
# Field Observations -- Tracers
#
#-------------------------------------------------------------------
#
# Use the recharge temp and excess air corrected mixing ratio ensembles
# These are output from Age_Modeling directory from age_modeling_mcmc.prep.py
# CFC and SF6 in pptv, 3H in TU, and 4He_ter in cm3/g
#
map_dict   = pd.read_pickle('./ER_PLM_ParFLow/Field_Data/map_dict.pk')
ens_dict   = pd.read_pickle('./ER_PLM_ParFLow/Field_Data/ens_dict.pk')



#
# Utility Funcitons
#
def pf_2_dates(startdate, enddate, f):
    '''Assumes ParFlow outputs every 24 hours'''
    s = pd.to_datetime(startdate)
    e = pd.to_datetime(enddate)
    d_list = pd.date_range(start=s, end=e, freq=f)
    # Drop Leap years again
    #d_list_ = d_list[~((d_list.month == 2) & (d_list.day == 29))]
    return d_list

dates    = pf_2_dates('2016-10-01', '2021-10-01', '24H')
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



#------------------------------------------
#
# Forcing MET data
#
#------------------------------------------
# spinup
met_sp    = pd.read_csv('./MET/met.avg.10yr.1hr.txt', delim_whitespace=True, names=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])
tstart    = pd.to_datetime('1989-10-01 00', format='%Y-%m-%d %H')
tend      = pd.to_datetime('1999-10-08 23', format='%Y-%m-%d %H') # Water year 21 is not over yet
hours     = pd.DatetimeIndex(pd.Series(pd.date_range(tstart, tend, freq='1H')))
#hours     = hours[~((hours.month == 2) & (hours.day == 29))] # No leap years
met_sp.index = hours
# 2017-2021
met_17_21 = pd.read_csv('./MET/met.2017-2021.1hr.txt', delim_whitespace=True, names=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])
tstart    = pd.to_datetime('2016-10-01 00', format='%Y-%m-%d %H')
tend      = pd.to_datetime('2021-09-30 23', format='%Y-%m-%d %H') # Water year 21 is not over yet
hours     = pd.DatetimeIndex(pd.Series(pd.date_range(tstart, tend, freq='1H')))
met_17_21.index = hours
# 2000-2016
met_00_16 = pd.read_csv('./MET/met.2000-2016.1hr.txt', delim_whitespace=True, names=['rad_s','rad_l','prcp','temp','wnd_u','wnd_v','press','vap'])
tstart    = pd.to_datetime('1999-10-01 00', format='%Y-%m-%d %H')
tend      = pd.to_datetime('2016-09-30 23', format='%Y-%m-%d %H') # Water year 21 is not over yet
hours     = pd.DatetimeIndex(pd.Series(pd.date_range(tstart, tend, freq='1H')))
met_00_16.index = hours
# 1979-2021
met_comp = pd.concat((met_00_16, met_17_21))


#----------------------------------------
#
# Summarize Precipitation
#
#----------------------------------------
# monthly (mm) of precipitation
# spinup
prcp0_summ       = (met_sp['prcp']*3600*1).groupby(pd.Grouper(freq='M')).sum() # mm/month
prcp0_summ.index = prcp0_summ.index.map(lambda x: x.replace(day=1))
# 2000 - 2016
prcp1_summ       = (met_00_16['prcp']*3600*1).groupby(pd.Grouper(freq='M')).sum() # mm/month
prcp1_summ.index = prcp1_summ.index.map(lambda x: x.replace(day=1))
# 2017-2021
prcp2_summ       = (met_17_21['prcp']*3600*1).groupby(pd.Grouper(freq='M')).sum() # mm/month
prcp2_summ.index = prcp2_summ.index.map(lambda x: x.replace(day=1))

# yearly accumulated precip
prcp0_sumy = (met_sp['prcp']*3600*1).groupby(pd.Grouper(freq='Y')).sum() # mm/yr
prcp1_sumy = (met_00_16['prcp']*3600*1).groupby(pd.Grouper(freq='Y')).sum() # mm/yr
prcp2_sumy = (met_17_21['prcp']*3600*1).groupby(pd.Grouper(freq='Y')).sum() # mm/yr

# Use the yearly one and clip; prevents skewed months because water year starts october
prcp_summ       = (met_comp['prcp']*3600*1).groupby(pd.Grouper(freq='M')).sum() # mm/month
prcp_summ.index = prcp_summ.index.map(lambda x: x.replace(day=1))
prcp_sumy       = (met_comp['prcp']*3600*1).groupby(pd.Grouper(freq='Y')).sum() # mm/yr

# annual wy precip totals
prcp_wy = pd.DataFrame((met_17_21['prcp']*3600*1).groupby(pd.Grouper(freq='D')).sum()) # mm/day
prcp_wy['wy'] = set_wy(prcp_wy)[0]
prcp_wy = prcp_wy.groupby(by='wy').sum()






met_daily = pd.DataFrame(met_17_21.copy()['prcp']*3600)
met_daily = met_daily.groupby(pd.Grouper(freq='D')).sum()
met_daily['wy'] = set_wy(met_daily)[0]
met_cs = met_daily.groupby(by='wy').cumsum()




#
# --- Pick sets to plot ---
#
dirlist = ['Run.wl.tetsu22.sh', 'Run.wl.127.sh', 'Run.wl.082.sh', 'Run.wl.104.sh', 'Run.wl.037.sh', 'Run.wl.101.sh']
labs    = ['Base', 'A1', 'A2', 'A3', 'A4', 'A5']
fdir = 'figures.opt.manuscript'

#dirlist = ['Run.wl.127.sh', 'Run.wl.127.sh.hk', 'Run.wl.127.sh.lk', 'Run.wl.127.sh.expk1', 'Run.wl.127.sh.lpor', 'Run.wl.127.sh.hpor']
#labs    = ['A1', 'High K', 'Low K', 'Exp. D. K', 'Low por', 'High por']
#fdir = 'figures.127.manuscript'


# Create a new directory for figures
if os.path.exists(fdir) and os.path.isdir(fdir):
    pass
else:
    os.makedirs(fdir)




#-------------------------------------------------------
#
# Water Level Plots
#
#-------------------------------------------------------
w = ['PLM1','PLM7','PLM6']

cm = plt.cm.tab10(np.linspace(0,1,10))

# --- WL 2017-2021 Plot ---
wells  = {'PLM1':404, 'PLM6':494}
w  = list(wells.keys())
w_ = list(wells.values())

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7.0, 2.2))
fig.subplots_adjust(hspace=0.15, top=0.96, right=0.75, bottom=0.12, left=0.15)
#fig, axes = plt.subplots(nrows=et_.shape[0]+1, ncols=et_.shape[1], figsize=(7.3,6))
#fig.subplots_adjust(top=0.86, bottom=0.12, left=0.08, right=0.97, hspace=0.15, wspace=0.1)
w = ['PLM1','PLM6']
# Plot water levels
for i in range(len(w)):
    for j in range(len(dirlist)):
        #wll = pd.read_csv('./{}/parflow_out/wy_2017_2021_wt_bls.csv'.format(dirlist[j]))
        #ax[i].plot(dates, wll[w[i]], color=cm[j], alpha=0.75, label=labs[j])
        wll = pd.read_pickle('{}/parflow_out/pf_out_dict_0021.pk'.format(dirlist[j]))['wtd']
        ax[i].plot(dates, wll[-len(dates):,w_[i]], color=cm[j], alpha=0.7, label=labs[j])
    ax[i].invert_yaxis()
    # Now add field observations
    if w[i] == 'PLM1':
        ax[i].plot(plm1_obs['bls'], color='black', alpha=1.0, linewidth=1.75, linestyle='--', label='Obs.')
    elif w[i] == 'PLM6':
        ax[i].plot(plm6_obs['bls'], color='black', alpha=1.0, linewidth=1.75, linestyle='--', label='Obs.')
    # Fill between years
    for ii in np.arange(2017,2021).reshape(2,2):
        ax[i].axvspan(pd.to_datetime('{}-10-01'.format(ii[0])), pd.to_datetime('{}-09-30'.format(ii[1])), alpha=0.04, color='red')
    ax[i].set_ylim(5.0, 0.0)
    #ax[i].yaxis.set_major_locator(MaxNLocator(4))
    ax[i].yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax[i].yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    #ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    # Xticks as october 1st
    ax[i].xaxis.set_major_locator(mdates.YearLocator(base=1, month=10, day=1))
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    #ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7]))
    ax[i].xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    #    label.set(horizontalalignment='right', rotation_mode="anchor")#, horizontalalignment='right')
    ax[i].set_ylabel(w[i])
    #ax[i].text(0.1, 0.15, w[i], horizontalalignment='center', verticalalignment='center', transform=ax[i].transAxes,bbox=dict(facecolor='white', edgecolor='white', alpha=0.6))
    #
    ax[i].tick_params(axis='x', bottom=True, top=False)
    ax[i].tick_params(axis='y', which='both', right=False)
    ax[i].margins(x=0.01)
    ax[i].grid()
    ax[i].set_xlim(pd.to_datetime('2016-09-30'), pd.to_datetime('2021-10-02'))
# Clean up
fig.text(0.035, 0.55, 'Water Table\nDepth (m)', va='center', ha='center', rotation='vertical')
ax[0].tick_params(axis='x', labelbottom=False)
ax[1].tick_params(axis='x', labelbottom=True, rotation=0)#, pad=0.1)
ax[0].legend(loc='upper left', bbox_to_anchor=(0.99, 1.08), handlelength=0.9, labelspacing=0.25, handletextpad=0.25)
#plt.savefig('./figures/waterlevels_comp.1721.jpg', dpi=300)
plt.savefig(os.path.join(fdir,'waterlevels_comp.1721.jpg'), dpi=300)
plt.savefig(os.path.join(fdir,'waterlevels_comp.1721.svg'), format='svg')
plt.show()








#--------------------------------------------
#
# A1-A5 Parameter Values
#
#--------------------------------------------
pars = pd.read_csv('../PLM_transect.wl/opt_models_wl.csv', index_col='model')
pars.index = ['Base','A1','A2','A3','A4','A5']

_labmap = {'K_soil' : r'$log_{10}$''\n''$K_{soil}$',
           'K_wshale' : r'$log_{10} K_{wshale}$',
           'K_fshale' : r'$log_{10} K_{fshale}$',
           'soil_rel_alpha' : r'$\alpha_{soil}$',
           'wshale_rel_alpha' : r'$\alpha_{wshale}$',
           'soil_rel_n' : r'$n_{soil}$',
           'wshale_rel_n' : r'$n_{wshale}$'}


fig, axes = plt.subplots(ncols=1, nrows=pars.shape[1], figsize=(3.0, 2.2))
for i,p in enumerate(pars.columns):
    ax  = axes[i]
    _dd = pars.loc[:,p] 
    for j,n in enumerate(_dd.index):
        ax.scatter(_dd[n], 1.0, color='C{}'.format(j), marker='o', label=n)
    
    ax.set_ylim(0.9, 1.1)
    ax.tick_params(axis='y', labelleft=False, left=False)
    if p in ['K_soil', 'K_wshale', 'K_fshale']:
        ax.set_xlim(-8, -3)
        ax.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False)
        if p != 'K_soil':
            ax.tick_params(axis='x', top=True, labeltop=False, bottom=False, labelbottom=False)
    else:
        ax.set_xlim(0, 5)
        ax.tick_params(axis='x', top=False, labeltop=False, bottom=True, labelbottom=True)
        if p != 'wshale_rel_n':
            ax.tick_params(axis='x', top=False, labeltop=False, bottom=True, labelbottom=False)
    
    ax.set_ylabel(_labmap[p], rotation=70)
fig.tight_layout()
plt.show()










#--------------------------------------------
#
# Ecoslim Age Distributions Plotting
#
#--------------------------------------------
def flux_wt_rtd(rtd_dict_unsort, model_time, well_name, nbins):
    '''Convert particle ages at a single model time and well location to a residence time distribution.
    Returns:
        - rtd_df:   Dataframe with mass weighted ages for all paricles (sorted by age)
        - rtd_dfs:  Dataframe similar to above but bins age distribution into discrete intervals  
    Inputs: 
        - rtd_dict_unsort: output from ecoslim_pnts_vtk.read_vtk() above
        - model_time: model time to consider. Must be in rtd_dict_unsort.keys()
        - well_name: observation well to consider. Must be in rtd_dict_unsort.keys()
        - nbins: Number of intervals to bin the ages into.'''
    # Flux Weighted RTD
    # Info regarding particles at a single timestep and single point (box) of the domain
    #pdb.set_trace()
    rtd    = rtd_dict_unsort[model_time][well_name] 
    rtd_df = pd.DataFrame(data=rtd,columns=['Time','Mass','Source','Xin'])
    rtd_df['wt'] = rtd_df['Mass']/rtd_df['Mass'].sum()
    rtd_df.sort_values('Time', inplace=True)
    rtd_df['Time'] /= 8760
    
    # Now some binning
    #nbins = 10
    gb =  rtd_df.groupby(pd.cut(rtd_df['Time'], nbins))
    rtd_dfs = gb.agg(dict(Time='mean',Mass='sum',Source='mean',Xin='mean',wt='sum'))
    rtd_dfs['count'] = gb.count()['Time']
    rtd_dfs['Time'] = rtd_dfs['Time'].interpolate(method='linear')
    return rtd_df, rtd_dfs



#------------
def build_timelist(yrs_list):
    #pdb.set_trace()
    #yr = np.arange(2017,2022)
    yr = yrs_list
    wy_inds_   = [np.where((date_map['Date'] > '{}-09-30'.format(i-1)) & (date_map['Date'] < '{}-10-01'.format(i)), True, False) for i in yr]
    wy_inds    = [date_map.index[wy_inds_[i]] for i in range(len(yr))]
    wy_map     = dict(zip(yr, wy_inds))
    wy_inds    = np.concatenate(((wy_inds)))
    time_list  = list(wy_inds[np.isin(wy_inds, list(rtd_dict.keys()))]) # model times that I want
    return time_list



# Pull Data For Single Conceptual model
fpath = dirlist[0]
rtd_dict = pd.read_pickle(os.path.join(fpath,'parflow_out', 'ecoslim_rtd.1721.pk'))


samp_date = '2021-05-11'
model_time = date_map[date_map['Date'] == samp_date].index[0]
model_time = list(rtd_dict.keys())[abs(list(rtd_dict.keys()) - model_time).argmin()]
model_time_samp = model_time


#
# Ensemble CDFs in seperate subplots
#
time_list = build_timelist([2021])
# Generate a Dictionary with all the RTDS and the taus for the chosen years and wells
def build_rtd_trans(rtd_dict, wells, time_list):
    rtd_trans = {}
    for i in range(len(wells)):
        w       = wells[i]
        tau_mu  = []
        tau_med = []
        rtd_df_ = []
        keep_list = [] # needed for when there are no particles for a timestep
        rtd_trans[w] = {}
        for t in range(len(time_list)):
            model_time = time_list[t]
            try:
                rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, model_time, w, 10)
                rtd_df_.append(rtd_df)
                tau_mu.append((rtd_df['Time'] * rtd_df['wt']).sum())
                tau_med.append(np.median(rtd_df['Time']))
                # not sure here, some NANs where there are zero particles
                #rtd_dfs['Time'] = rtd_dfs['Time'].interpolate(method='linear')
                keep_list.append(t)
            except ValueError:
                print (w, model_time)
                pass
        rtd_trans[w]['tau_mu']    = tau_mu
        rtd_trans[w]['tau_med']   = tau_med
        rtd_trans[w]['rtd_df']    = rtd_df_
        rtd_trans[w]['keep_list'] = keep_list
    return rtd_trans



# Ensemble Plot -- includes different timestamps
wells = ['PLM1','PLM7','PLM6']
colors = plt.cm.coolwarm(np.linspace(0,1,12+1))
fig, axes = plt.subplots(nrows=3, ncols=len(dirlist), figsize=(9.8, 5.2))
fig.subplots_adjust(wspace=0.20, hspace=0.38, top=0.94, bottom=0.13, right=0.88, left=0.1)
for i in range(3):
    w = wells[i]
    for j in range(len(dirlist)):
        ax = axes[i,j]
        rtd_dict  = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
        rtd_trans = build_rtd_trans(rtd_dict, wells, time_list)
        for t in range(len(time_list)):
            _date = date_map.loc[time_list[t],'Date']
            try:
                _rtd_df = rtd_trans[w]['rtd_df'][t]
                _tau    = rtd_trans[w]['tau_mu'][t]
                #ax.plot(_rtd_df['Time'], np.cumsum(_rtd_df['wt']), color='black', alpha=0.75, zorder=4)
                #ax.axvline(_tau, color='grey', alpha=0.3, linestyle='-', zorder=2)
                ax.plot(_rtd_df['Time'], np.cumsum(_rtd_df['wt']), color=colors[_date.month], alpha=0.7, zorder=4, label=_date.month)
            except IndexError:
                print ('no {} {}'.format(w, time_list[t]))
                pass
        # Plot mean ages as vertical lines
        #ax.axvspan(min(rtd_trans[w]['tau_mu']), max(rtd_trans[w]['tau_mu']), color='grey', alpha=0.3)
        # Replot the median 
        mid_ind  = np.where(np.array(rtd_trans[w]['tau_mu'])==np.sort(rtd_trans[w]['tau_mu'])[len(rtd_trans[w]['tau_mu'])//2])[0][0]
        _rtd_df_ = rtd_trans[w]['rtd_df'][mid_ind]
        _tau_    = rtd_trans[w]['tau_mu'][mid_ind]
        ax.plot(_rtd_df_['Time'], np.cumsum(_rtd_df_['wt']), color='black', lw=2.0, linestyle='-', alpha=0.90, zorder=6, label='mean')
        #ax.axvline(_tau_, ymin=0.0, ymax=0.5, color='black', alpha=0.65, linestyle='--', zorder=5)
        ax.axvline(_rtd_df_['Time'][abs(np.cumsum(_rtd_df_['wt'])-0.5).idxmin()], ymin=0.0, ymax=0.5, color='black', alpha=0.65, linestyle='--', zorder=5)
        # Cleanup 
        ax.minorticks_on()
        ax.grid()
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.xaxis.set_major_locator(MaxNLocator(2))
        ax.tick_params(axis='x', rotation=20, pad=0.1)
        if i == 0:
            ax.set_title(labs[j], fontsize=14, pad=0.1)
        if j == 0:     
            ax.set_ylabel('{}'.format(w))
        else:
            ax.tick_params(axis='y', labelleft=False)
fig.text(0.5, 0.015, 'Particle Age (years)', ha='center')
fig.text(0.0075, 0.53, 'CDF', va='center', rotation='vertical')
# legend
handles, labels = axes[0,len(dirlist)-1].get_legend_handles_labels()
unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
axes[0,len(dirlist)-1].legend(*zip(*unique), bbox_to_anchor=(0.98, 1.1), handlelength=1.0, labelspacing=0.25, handletextpad=0.5, title='month')
#plt.savefig('./figures/ecoslim_rtd_ens.jpg',dpi=300)
plt.savefig(os.path.join(fdir,'ecoslim_rtd_ens.jpg'), dpi=300)
plt.savefig(os.path.join(fdir,'ecoslim_rtd_ens.svg'), format='svg')
plt.show()




# Box-Plot of age dyanmics at the groundwater wells
# Big array of all particle ages over the entire year
fig_labs = {'PLM1':'A','PLM7':'B','PLM6':'C'}
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5, 3.0))
fig.subplots_adjust(hspace=0.16, top=0.96, bottom=0.2, right=0.94, left=0.22)
age_wted_glob = {}
off_ = -0.3
for i in range(3):
    w = wells[i]
    age_wted = []
    for j in range(len(dirlist)):
        rtd_dict  = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
        rtd_trans = build_rtd_trans(rtd_dict, wells, time_list)
        _age_arr  = []
        _mass_arr = []
        for t in range(len(time_list)):
            try:
                _age_arr.append(rtd_trans[w]['rtd_df'][t]['Time'].to_numpy())
                _mass_arr.append(rtd_trans[w]['rtd_df'][t]['Mass'].to_numpy())
            except IndexError:
                pass
        _age_arr  = np.concatenate(_age_arr)
        _mass_arr = np.concatenate(_mass_arr)
        #_age_arr  = np.concatenate(([rtd_trans[w]['rtd_df'][t]['Time'].to_numpy() for t in range(len(time_list))]))
        #_mass_arr = np.concatenate(([rtd_trans[w]['rtd_df'][t]['Mass'].to_numpy() for t in range(len(time_list))]))
        # Mass-weighted ages
        _age_df   = pd.DataFrame(np.column_stack((_age_arr,_mass_arr)), columns=['Age','Mass'])
        _age_df.sort_values(by='Age', inplace=True)
        # Bin Ages into intervals and accumulate Mass for each of the bins
        _age_bin = _age_df.groupby(by=pd.cut(_age_df['Age'], bins=60)).agg(dict(Age='mean', Mass='sum'))
        # Normalize mass between 1000 and 1
        _age_bin['Mass_norm'] = (_age_bin['Mass'] - _age_bin['Mass'].min()) / (_age_bin['Mass'].max()-_age_bin['Mass'].min())
        _age_bin['Mass_norm'] = _age_bin['Mass_norm']* (1000-1)+1
        _age_bin.dropna(inplace=True)
        # Now make array with age values occuring proportional to mass in that bin
        _age_wted = np.concatenate(([[_age_bin['Age'].iloc[i]]*int(_age_bin['Mass_norm'].iloc[i]) for i in range(len(_age_bin))]))
        age_wted.append(_age_wted)
    age_wted_glob[w] = age_wted
    print ('')
    print('{} min tau = {}'.format(w, np.concatenate(age_wted).min()))
    print('{} max tau = {}'.format(w, np.concatenate(age_wted).max()))
    print ('')
    for zz in range(len(age_wted)):
        print ('{} Percent SD {}'.format(w, 100*(age_wted[zz].std()/age_wted[zz].mean())))
    # Box-plot
    #ax = axes[i]
    xx = np.arange(len(age_wted))
    ax.boxplot(x=age_wted, positions=xx+off_, widths=0.25, showfliers=True, 
               flierprops={'color':'grey','alpha':0.2,'markersize':4}, medianprops={'color':'C{}'.format(i),'linewidth':1.5})#, labels=[w]*len(xx))
    off_ += 0.3
    ax.set_yscale('log')
    ax.set_xticks(ticks=np.arange(len(dirlist)), labels=labs)
    ax.tick_params(axis='x', rotation=25, pad=0.025)
    #ax.tick_params(axis='x', rotation=0, pad=0.5) #pad=0.025)
    ax.grid()
    ax.set_ylim(10, 2000)
    ax.set_xlim(-0.5, len(dirlist)-0.5)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=20))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(0.1,1,0.1), numticks=20))
    ax.set_ylabel('Simulated Particle Ages\n(years)')
    # Annoying, but plot lines off yscale for legend
    ax.scatter(x=-10, y=-10, marker='_', color='C{}'.format(i), s=100, label=w)
    for ii in np.arange(0,6).reshape(3,2):
        ax.axvspan(ii[0]-0.5, ii[1]-0.5, alpha=0.02, color='red')
ax.legend(handlelength=1.0, labelspacing=0.25, handletextpad=0.1)
plt.savefig(os.path.join(fdir,'ecoslim_rtd_ens.simple.jpg'), dpi=300)
plt.show()


# percent change from PLM1 to PLM6
pc = []
for i in range(len(age_wted)):
    pchange = (age_wted_glob['PLM6'][i].mean()-age_wted_glob['PLM1'][i].mean())/age_wted_glob['PLM1'][i].mean()*100
    print ('{} -- {:.3f} % increase from PLM1 to PLM6 '.format(labs[i], pchange))
    pc.append(pchange)










#
# Soil Wells
#
wells = ['X404', 'X494', 'X508']
names = ['PLM1 Soil', 'PLM6 Soil', 'Floodplain']


# Ensemble Plot with different timestamps
colors = plt.cm.coolwarm(np.linspace(0,1,12+1))
fig, axes = plt.subplots(nrows=3, ncols=len(dirlist), figsize=(9.8, 5.2))
fig.subplots_adjust(wspace=0.20, hspace=0.38, top=0.94, bottom=0.13, right=0.88, left=0.1)
for i in range(3):
    w = wells[i]
    for j in range(len(dirlist)):
        ax = axes[i,j]
        rtd_dict  = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
        rtd_trans = build_rtd_trans(rtd_dict, wells, time_list)
        for t in range(len(time_list)):
            _date = date_map.loc[time_list[t],'Date']
            try:
                _rtd_df = rtd_trans[w]['rtd_df'][t]
                _tau    = rtd_trans[w]['tau_mu'][t]
                #ax.plot(_rtd_df['Time'], np.cumsum(_rtd_df['wt']), color='black', alpha=0.75, zorder=4)
                #ax.axvline(_tau, color='grey', alpha=0.3, linestyle='-', zorder=2)
                ax.plot(_rtd_df['Time'], np.cumsum(_rtd_df['wt']), color=colors[_date.month], alpha=0.7, zorder=4, label=_date.month)
            except IndexError:
                print ('no {} {}'.format(w, time_list[t]))
                pass
        # Plot mean ages as vertical lines
        #ax.axvspan(min(rtd_trans[w]['tau_mu']), max(rtd_trans[w]['tau_mu']), color='grey', alpha=0.3)
        # Replot the median 
        mid_ind  = np.where(np.array(rtd_trans[w]['tau_mu'])==np.sort(rtd_trans[w]['tau_mu'])[len(rtd_trans[w]['tau_mu'])//2])[0][0]
        _rtd_df_ = rtd_trans[w]['rtd_df'][mid_ind]
        _tau_    = rtd_trans[w]['tau_mu'][mid_ind]
        ax.plot(_rtd_df_['Time'], np.cumsum(_rtd_df_['wt']), color='black', lw=2.0, linestyle='-', alpha=0.90, zorder=6, label='mean')
        #ax.axvline(_tau_, ymin=0.0, ymax=0.5, color='black', alpha=0.65, linestyle='--', zorder=5)
        ax.axvline(_rtd_df_['Time'][abs(np.cumsum(_rtd_df_['wt'])-0.5).idxmin()], ymin=0.0, ymax=0.5, color='black', alpha=0.65, linestyle='--', zorder=5)
        # Cleanup 
        ax.minorticks_on()
        ax.grid()
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.xaxis.set_major_locator(MaxNLocator(2))
        ax.tick_params(axis='x', rotation=20, pad=0.1)
        if i == 0:
            ax.set_title(labs[j], fontsize=14, pad=0.1)
        if j == 0:     
            #ax.set_ylabel('{}'.format(w))
            ax.set_ylabel('{}'.format(names[i]))
        else:
            ax.tick_params(axis='y', labelleft=False)
fig.text(0.5, 0.015, 'Particle Age (years)', ha='center')
fig.text(0.0075, 0.53, 'CDF', va='center', rotation='vertical')
# legend
handles, labels = axes[0,len(dirlist)-1].get_legend_handles_labels()
unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
axes[0,len(dirlist)-1].legend(*zip(*unique), bbox_to_anchor=(0.98, 1.1), handlelength=1.0, labelspacing=0.25, handletextpad=0.5, title='month')
#plt.savefig('./figures/ecoslim_rtd_ens.jpg',dpi=300)
plt.savefig(os.path.join(fdir,'ecoslim_rtd_ens.soil.jpg') ,dpi=300)
plt.show()



















#-------------------------------------------------------
#
# Fraction Younger in Soil  Plots
#
#-------------------------------------------------------
date_map['wy'] = set_wy(date_map)[0]
date_map.loc[0,'wy'] = 2017

yr = np.arange(2017,2022)

wy_inds_   = [np.where((date_map['Date'] > '{}-09-30'.format(i-1)) & (date_map['Date'] < '{}-10-01'.format(i)), True, False) for i in yr]
wy_inds    = [date_map.index[wy_inds_[i]] for i in range(len(yr))]
wy_map     = dict(zip(yr, wy_inds))
wy_inds    = np.concatenate(((wy_inds)))

time_list       = list(wy_inds[np.isin(wy_inds, list(rtd_dict.keys()))]) #+ [model_time_samp] # model times that I want
time_list_dates = date_map.loc[time_list,'Date'].to_list() # dates
time_list_yr    = date_map.loc[time_list,'wy'].to_list() # corresponding water year
time_list_map_  = [np.where(np.array(time_list_yr)==y)[0] for y in np.unique(np.array(time_list_yr))]
time_list_map   = dict(zip(np.unique(np.array(time_list_yr)), time_list_map_))

# Want first day of the month for x labels -- pick a full year like 2019
month_map  = {1:'O',2:'N',3:'D',4:'J',5:'F',6:'M',7:'A',8:'M',9:'J',10:'J',11:'A',12:'S'}
months     = list(month_map.values())
days_ = [date_map['Date'][date_map['wy']==2019].iloc[i].month for i in range(len(date_map['Date'][date_map['wy']==2019]))]
first_month  = [(np.array(days_)==i).argmax() for i in [10,11,12,1,2,3,4,5,6,7,8,9]]

# Sample date to model time
samp_date   = '2021-05-11' # date to consider and extract RTD info for
samp_date_  = np.array([(i-pd.to_datetime(samp_date)).days for i in time_list_dates])
_model_time = abs(samp_date_).argmin()

# ---
wells = ['X404', 'X494', 'X508']
names = ['PLM1 Soil', 'PLM6 Soil', 'Floodplain']
names_dict = dict(zip(wells,names))



#
# --- Fraction Younger and Snow Fraction --- 
#
w = 'X508'

fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(6.5, 7.5))
fig.subplots_adjust(hspace=0.14, top=0.98, bottom=0.08, right=0.7, left=0.18)
# Get the data
replotbase   = False
snow_base    = None
fy_diff_base = None
for j in range(len(dirlist)):
    ax = axes[j]
    frac_young_ = {0.25:[],1:[],10:[],100:[],1000:[]}
    dates_      = []
    tau_        = []
    snow_       = [] # source=1 is intitial particles, source=2 is rain, source=3 is snow
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for t in range(len(time_list)):
        rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
        tau_.append((rtd_df['Time'] * rtd_df['wt']).sum())
        dates_.append(time_list_dates[t])
        snow_.append((rtd_df['Source']==3).sum() / len(rtd_df))
        for k in list(frac_young_.keys()):
            frac_young_[k].append(rtd_df[rtd_df['Time']<k]['wt'].sum())
    fkeys = list(frac_young_.keys())
    # Find cumulative fraction younger
    fy = 100*pd.DataFrame.from_dict(frac_young_).T
    fy_diff = pd.DataFrame([fy[i].sub(fy[i].shift()) for i in fy.columns])
    fy_diff.iloc[:,0] = fy.iloc[0,:]
    fy_diff[1001.0] = 100-fy_diff.sum(axis=1) # add greater than 1000
    fy_diff_ = fy_diff.to_numpy().T
    ax.stackplot(dates_, fy_diff_, colors=plt.cm.coolwarm(np.linspace(0,1,len(frac_young_)+1)), labels=list(frac_young_.keys())+ ['>{}'.format(list(frac_young_.keys())[-1])], alpha=0.7)
    #
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.tick_params(axis='x', rotation=20, pad=0.1)
    ax.text(0.015, 0.88, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.0))
    ax.set_ylim(-0.05,150)
    ax.set_yticks(ticks=[0,50,100], labels=[0,50,100])
    ax.set_yticks(ticks=[10,20,30,40,60,70,80,90], minor=True)
    #ax.grid()
    ax.margins(x=0.01)
    # For legend...
    ax.plot(dates_, np.ones(len(dates_))+180, color='black', alpha=0.9, linestyle='-', label='Snow' if j==0 else '')
    # --- plot snow versus rain ---
    ax2 = ax.twinx()
    ax2.plot(dates_, np.array(snow_)*100, color='black', alpha=0.9, linestyle='-', label='Snow Fraction' if j==0 else '')
    ax2.set_ylim(-5,175)
    ax2.set_yticks(ticks=[0,50,100], labels=[0,50,100])
    ax2.set_yticks(ticks=[10,20,30,40,60,70,80,90], minor=True)
    ax2.invert_yaxis()
    ax2.axhline(50, color='grey', lw=1.0, alpha=0.5)
    #    
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='grey')
    # Re-plot base model for comparisons?
    if replotbase == True:
        if j == 0:
            snow_base = np.array(snow_).copy()
            fy_diff_base = fy_diff_.copy() 
            ax.plot(dates_, fy_diff_base[:2,:].sum(axis=0), color='black', linestyle=':', alpha=0.65) 
            ax2.plot(dates_, snow_base*100, color='grey', alpha=0.65, linestyle='--')
        if j !=0 :
            # First two columns in fy_diff_base are the 0.25 and 1 fractions -- replot these
            ax.plot(dates_, fy_diff_base[:2,:].sum(axis=0), color='black', linestyle=':', alpha=0.65)  
            ax2.plot(dates_, snow_base*100, color='grey', alpha=0.75, linestyle='--')  
#fig.text(0.96, 0.50, 'Mean Age (yrs)'.format(names_dict[w]), va='center', rotation='vertical')
fig.text(0.06, 0.50, 'Fraction Younger (%)'.format(list(frac_young_.keys())[0]), va='center', rotation='vertical', color='black')
fig.text(0.79, 0.47, 'Fraction from Snow (%)', va='center', rotation='vertical', color='black')
axes[0].legend(ncol=1, loc='upper left', bbox_to_anchor=(1.15, 1.15),
          handlelength=0.8, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13, title='Fraction\nYounger\n(years)') #fig.tight_layout()
plt.savefig(os.path.join(fdir,'fraction_younger.png'), dpi=300)
plt.savefig(os.path.join(fdir,'fraction_younger.svg'), format='svg')
plt.show()


# --- Scatter plot to compare snow fraction across models ---
snow = []
for j in range(len(dirlist)):
    _snow = []
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for t in range(len(time_list)):
        rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
        _snow.append((rtd_df['Source']==3).sum() / len(rtd_df))
    snow.append(_snow)
snow = np.array(snow).T
#
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.2))
fig.subplots_adjust(hspace=0.14, top=0.95, bottom=0.36, right=0.98, left=0.28)
bplot1 = ax.boxplot(snow*100, patch_artist=True)
colors = ['C{}'.format(i) for i in range(len(dirlist))]
for i in range(len(bplot1)):
    bplot1['boxes'][i].set_facecolor(colors[i])
    bplot1['medians'][i].set_color('black')   
ax.set_xticks(ticks=np.arange(len(dirlist))+1, labels=labs)
ax.tick_params(axis='x', rotation=30, pad=0.1)
ax.grid(axis='y')
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.set_ylabel('Fraction from\nSnow (%)', labelpad=1.2)
#ax.yaxis.set_label_coords(-0.2, 0.4)
#for label in ax.get_xticklabels(which='major'):
#    label.set(horizontalalignment='right', rotation_mode="anchor")
plt.savefig(os.path.join(fdir, 'fraction_younger_snow_simple.png'), dpi=300)
plt.savefig(os.path.join(fdir, 'fraction_younger_snow_simple.svg'), format='svg')
plt.show()


# --- Updated Snow vs. SWE plot ---
snow = []
for j in range(len(dirlist)):
    _snow = []
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for t in range(len(time_list)):
        rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
        _snow.append((rtd_df['Source']==3).sum() / len(rtd_df))
    snow.append(_snow)
snow = np.array(snow).T
# Seperate snow into each WY
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5, 2.7))
fig.subplots_adjust(hspace=0.14, top=0.95, bottom=0.36, right=0.75, left=0.15)
_p = prcp_wy.copy().to_numpy()
_offset = np.linspace(-25,25,6)
_wy_marker = ['o','^','X','s','v']
_wy_marker = ['o']*6 
for j,m in enumerate(labs):
    for i,y in enumerate(list(time_list_map.keys())):
        _snow = snow[time_list_map[y],j] * 100
        _snow_unc = abs(_snow.mean() - np.array([np.percentile(_snow,10), np.percentile(_snow,90)]))[:,np.newaxis] 
        ax.errorbar(_p[i]+_offset[j], y=_snow.mean(), yerr=_snow_unc, alpha=0.6, color='C{}'.format(j), zorder=9)
        ax.scatter(_p[i]+_offset[j], _snow.mean(), color='C{}'.format(j), marker=_wy_marker[i], zorder=10)
# Generate manual legends for simplicity
for j,m in enumerate(labs): #colors for each model 
    ax.scatter(-100, -100, marker='o', color='C{}'.format(j), label=m)
for i,y in enumerate(list(time_list_map.keys())): #shapes for WYs
    #ax.scatter(-100, -100, marker=_wy_marker[i], color='black', label=y)
    #ax.vlines(x=_p[i], ymin=0, ymax=30, color='grey', linestyle='--', alpha=0.75, zorder=8)
    ax.axvspan(float(_p[i]+_offset[0]-2), float(_p[i]+_offset[-1]+2), color='grey', alpha=0.20)
    if y == 2018:
        ax.text(x=_p[i], y=20, s='WY{}'.format(str(y)[2:]), fontsize=12, va='bottom', ha='center')
    else:
        ax.text(x=_p[i], y=20, s=str(y)[2:], fontsize=12, va='bottom', ha='center')
ax.grid(axis='y')
ax.minorticks_on()
ax.set_ylim(19,90)
ax.set_xlim(340, 770)
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax.set_xticks(ticks=_p.ravel(), labels=_p.ravel().astype(int))
ax.set_ylabel('Fraction from\nSnow (%)', labelpad=1.2)
ax.set_xlabel('Annual Precip. (mm/year)')
#ax.yaxis.set_label_coords(-0.2, 0.4)
#for label in ax.get_xticklabels(which='major'):
#    label.set(horizontalalignment='right', rotation_mode="anchor")
ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.0, 1.1), handlelength=1.5, labelspacing=0.25, columnspacing=0.25, handletextpad=0.05, frameon=False)
plt.savefig(os.path.join(fdir, 'eco_snow_vs_p.png'), dpi=300)
plt.savefig(os.path.join(fdir, 'eco_snow_vs_p.svg'), format='svg')
plt.show()




# --- Age CDF to easily compare models --- 
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5, 2.7))
fig.subplots_adjust(hspace=0.14, top=0.95, bottom=0.36, right=0.75, left=0.15)
for j in range(len(dirlist)):
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, model_time, w, 30)
    tau = (rtd_df['Time'] * rtd_df['wt']).sum()
    tau_med = np.median(rtd_df['Time'])
    # not sure here, some NANs where there are zero particles
    rtd_dfs['Time'] = rtd_dfs['Time'].interpolate(method='linear')
    # CDF plot
    ax.plot(rtd_df['Time'],  np.cumsum(rtd_df['wt']), linewidth=1.5, alpha=1.0, color=cm[j], label=labs[j])
    #ax.axvline(tau, color='grey', linestyle='--')
ax.set_ylabel('CDF')
# Clean up a bit
#matcks = ax.get_xticks()
#ax.xaxis.set_minor_locator(MultipleLocator((matcks[1]-matcks[0])/2))
#matcks = ax.get_xticks()
#ax.xaxis.set_minor_locator(MultipleLocator((matcks[1]-matcks[0])/2))
ax.minorticks_on()
ax.set_xscale('log')
#ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs='all', numticks=6))
ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=20))
ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(0.1,1,0.1), numticks=20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.tick_params(axis='x', which='both', rotation=0)#, pad=0.025)
ax.grid()
ax.set_xlim(0.01, 2000)
ax.set_xlabel('Particle Age (years)')
#ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.1), handlelength=1.5, labelspacing=0.25, handletextpad=0.5)
plt.savefig(os.path.join(fdir,'fraction_younger_cdf_simple.png'),dpi=300)
plt.savefig(os.path.join(fdir,'fraction_younger_cdf_simple.svg'),format='svg')
plt.show()



#
# --- Fraction Younger as Masses - not fractions ---
#
fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(6.5, 7.5))
fig.subplots_adjust(hspace=0.14, top=0.98, bottom=0.08, right=0.7, left=0.18)
# Get the data
for j in range(len(dirlist)):
    ax = axes[j]
    frac_young_ = {0.25:[],1:[],10:[],100:[],100:[],1000:[]}
    dates_      = []
    tau_        = []
    mass_       = []  # total particle mass in the control volume
    snow_       = [] # source=1 is intitial particles, source=2 is rain, source=3 is snow
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for t in range(len(time_list)):
        rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
        tau_.append((rtd_df['Time'] * rtd_df['wt']).sum())
        dates_.append(time_list_dates[t])
        mass_.append(rtd_df['Mass'].sum())
        snow_.append((rtd_df['Source']==3).sum() / len(rtd_df))
        for k in list(frac_young_.keys()):
            frac_young_[k].append(rtd_df[rtd_df['Time']<k]['Mass'].sum())
            
    # Find cumulative fraction younger
    fy = pd.DataFrame.from_dict(frac_young_).T
    fy.loc[list(frac_young_.keys())[-1]+1, :] = mass_
    fy_diff = pd.DataFrame([fy[i].sub(fy[i].shift()) for i in fy.columns])
    fy_diff.iloc[:,0] = fy.iloc[0,:]
    fy_diff_ = fy_diff.to_numpy().T
    ax.stackplot(dates_, fy_diff_/1000, colors=plt.cm.coolwarm(np.linspace(0,1,len(frac_young_)+1)), labels=list(frac_young_.keys()) + ['>{}'.format(list(frac_young_.keys())[-1])], alpha=0.7)
        
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.tick_params(axis='x', rotation=20, pad=0.1)
    ax.text(0.015, 0.88, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.0))
    #if j in [0,1,4]:
    #    ax.set_ylim(0,14.9)
    #    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    #    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    #elif j == 2:
    #    ax.set_ylim(0,6.9)
    #    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    #    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    #elif j in [3,5]:
    #    ax.set_ylim(0,65)
    #    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    #    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_ylim(0, fy.to_numpy().max()/1000 * 1.05)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.margins(x=0.01)
    # For legend...
    ax.plot(dates_, np.ones(len(dates_))+180, color='black', alpha=0.9, linestyle='-', label='Snow' if j==0 else '')
    # --- plot snow versus rain ---
    ax2 = ax.twinx()
    ax2.plot(dates_, np.array(snow_)*100, color='black', alpha=0.9, linestyle='-', label='Snow Fraction' if j==0 else '')
    ax2.set_ylim(-5,175)
    ax2.set_yticks(ticks=[0,50,100], labels=[0,50,100])
    ax2.set_yticks(ticks=[10,20,30,40,60,70,80,90], minor=True)
    ax2.invert_yaxis()
    #ax2.grid()
    ax2.axhline(50, color='grey', lw=1.0, alpha=0.5)
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='grey')
#fig.text(0.96, 0.50, 'Mean Age (yrs)'.format(names_dict[w]), va='center', rotation='vertical')
fig.text(0.08, 0.50, 'Total Mass (kg)'.format(list(frac_young_.keys())[0]), va='center', rotation='vertical', color='black')
fig.text(0.79, 0.47, 'Fraction from Snow (%)', va='center', rotation='vertical', color='black')
axes[0].legend(ncol=1, loc='upper left', bbox_to_anchor=(1.15, 1.15),
          handlelength=0.8, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13, title='Fraction\nYounger\n(years)') #fig.tight_layout()
plt.savefig(os.path.join(fdir,'fraction_younger_mass_q.png'), dpi=300)
plt.savefig(os.path.join(fdir,'fraction_younger_mass_q.svg'), format='svg')
plt.show()



# --- Scatter plot to compare the total masses ---
mass = []
for j in range(len(dirlist)):
    _mass  = []  # total particle mass in the control volume
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for t in range(len(time_list)):
        rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
        _mass.append(rtd_df['Mass'].sum())
    mass.append(_mass)    
mass = np.array(mass).T
#
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.2))
fig.subplots_adjust(hspace=0.14, top=0.95, bottom=0.36, right=0.98, left=0.28)
bplot1 = ax.boxplot(mass/1000, patch_artist=True)
colors = ['C{}'.format(i) for i in range(len(dirlist))]
for i in range(len(bplot1)):
    bplot1['boxes'][i].set_facecolor(colors[i])
    bplot1['medians'][i].set_color('black')
ax.set_xticks(ticks=np.arange(len(dirlist))+1, labels=labs)
ax.tick_params(axis='x', rotation=30, pad=0.025)
ax.grid(axis='y')
ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.set_ylabel('Total Mass (kg)')
#ax.yaxis.set_label_coords(-0.2, 0.45)
for label in ax.get_xticklabels(which='major'):
    label.set(horizontalalignment='right', rotation_mode="anchor")
plt.savefig(os.path.join(fdir,'fraction_younger_mass_simple.png'), dpi=300)
plt.show()







#
# --- Mean age plots ---
#
fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(6.5, 7.5))
fig.subplots_adjust(hspace=0.14, top=0.98, bottom=0.08, right=0.7, left=0.18)
cc = ['tab:grey', 'tab:blue', 'tab:brown']
ll = ['-', '--', ':']
# Get the data
wells = ['X508', 'PLM1', 'PLM6']
names = ['Soil Well', 'PLM1', 'PLM6']
for j in range(len(dirlist)):
    ax = axes[j]
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for i,w in enumerate(wells):
        tau_   = []
        dates_ = []
        for t in range(len(time_list)):
            rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
            tau_.append((rtd_df['Time'] * rtd_df['wt']).sum())
            dates_.append(time_list_dates[t])
        # plot mean age simulations
        ax.plot(dates_, tau_, color=cc[i], linestyle=ll[i], alpha=1.0, linewidth=2.0, label=names[i]if j==0 else '')
        # add PLM1 and PLM6 field obs -- mean ages from BMM?
    # Cleanup
    #ax.set_ylim(ax.get_ylim()[0]*0.9, ax.get_ylim()[1]*1.1)
    # want somewhat consistent axis
    if labs[j] in ['Base','A2']:
        ax.set_ylim(200, 1100)
    elif labs[j] in ['A5']:
        ax.set_ylim(10, 60)
    else:
        ax.set_ylim(20,200)
    #else:
    #    pass
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.tick_params(axis='x', rotation=20, pad=0.1)
    #ax.set_ylabel(labs[j], labelpad=1.0)
    ax.text(0.015, 0.88, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.0))
    ax.grid()
    ax.margins(x=0.01)
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='grey')
fig.text(0.055, 0.50, 'Mean Age (years)', va='center', rotation='vertical', color='black')
axes[0].legend(ncol=1, loc='upper left', bbox_to_anchor=(1.005, 1.12),
          handlelength=1.2, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13)
#fig.tight_layout()
plt.savefig(os.path.join(fdir,'fraction_tau.png'), dpi=300)
plt.savefig(os.path.join(fdir,'fraction_tau.svg'), format='svg')
plt.show()



# --- Mean Age Plot update ---
# add the tracer-derived mean ages
from mpl_toolkits.axes_grid1 import make_axes_locatable
epm = pd.read_csv('epm_post_comp.csv', index_col=0)
bmm = pd.read_csv('bmm_post_comp.csv', index_col=0)

fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(6.5, 7.5))
fig.subplots_adjust(hspace=0.14, top=0.98, bottom=0.08, right=0.7, left=0.18)
cc = ['tab:grey', 'tab:blue', 'tab:brown']
ll = ['-', '--', ':']
# Get the data
wells = ['X508', 'PLM1', 'PLM6']
names = ['Soil Well', 'PLM1', 'PLM6']
for j,n in enumerate(dirlist):
    ax = axes[j]
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for i,w in enumerate(wells):
        tau_   = []
        dates_ = []
        for t in range(len(time_list)):
            rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
            tau_.append((rtd_df['Time'] * rtd_df['wt']).sum())
            dates_.append(time_list_dates[t])
        # plot mean age simulations
        ax.plot(dates_, tau_, color=cc[i], linestyle=ll[i], alpha=1.0, linewidth=2.0, label=names[i]if j==0 else '')
    
    # add PLM1 and PLM6 field obs -- mean ages from BMM?
    if labs[j] in ['A1','A3']:
        ax.scatter(pd.to_datetime('2019-05-01'),  epm.loc['tau1','PLM1'], marker='^', color='tab:blue', s=60, zorder=10, label='PLM1 EPM')
        ax.errorbar(pd.to_datetime('2019-05-01'), epm.loc['tau1','PLM1'], yerr=abs(epm.loc['tau1','PLM1'] - np.array([epm.loc['tau1_10th','PLM1'], epm.loc['tau1_90th','PLM1']])[:,np.newaxis]), color='tab:blue', zorder=9)
        
        ax.scatter(pd.to_datetime('2019-05-01'),  bmm.loc['tau','PLM1'], marker='^', facecolors='white', s=60, edgecolors='tab:blue', zorder=10, label='PLM1 BMM')
        ax.errorbar(pd.to_datetime('2019-05-01'), bmm.loc['tau','PLM1'], yerr=abs(bmm.loc['tau','PLM1'] - np.array([bmm.loc['tau_10th','PLM1'], bmm.loc['tau_90th','PLM1']])[:,np.newaxis]), color='tab:blue', zorder=9)
        
        ax.scatter(pd.to_datetime('2019-06-30'),  epm.loc['tau1','PLM6'], marker='o', color='tab:brown', zorder=10, label='PLM6 EPM')
        ax.errorbar(pd.to_datetime('2019-06-30'), epm.loc['tau1','PLM6'], yerr=abs(epm.loc['tau1','PLM6'] - np.array([epm.loc['tau1_10th','PLM6'], epm.loc['tau1_90th','PLM6']])[:,np.newaxis]), color='tab:brown', zorder=9)
    
    if labs[j] in ['A4', 'A5']:
        ax.scatter(pd.to_datetime('2019-05-01'),  epm.loc['tau1','PLM1'], marker='^', color='tab:blue', s=60, zorder=10, label='PLM1 EPM')
        ax.errorbar(pd.to_datetime('2019-05-01'), epm.loc['tau1','PLM1'], yerr=abs(epm.loc['tau1','PLM1'] - np.array([epm.loc['tau1_10th','PLM1'], epm.loc['tau1_90th','PLM1']])[:,np.newaxis]), color='tab:blue', zorder=9)
        
    if labs[j] in ['Base', 'A2']:
        ax.scatter(pd.to_datetime('2019-05-01'),  bmm.loc['tau','PLM1'], marker='^', facecolors='white', s=60, edgecolors='tab:blue', zorder=10, label='PLM1 BMM')
        ax.errorbar(pd.to_datetime('2019-05-01'), bmm.loc['tau','PLM1'], yerr=abs(bmm.loc['tau','PLM1'] - np.array([bmm.loc['tau_10th','PLM1'], bmm.loc['tau_90th','PLM1']])[:,np.newaxis]), color='tab:blue', zorder=9)
        
    if labs[j] in ['Base','A2']:
        # try to add y-axis break
        ax.spines['top'].set_visible(False)
        divider = make_axes_locatable(ax)
        ax2 = divider.new_vertical(size="40%", pad=0.05)
        fig.add_axes(ax2)
        
        #BMM
        #ax2.scatter(pd.to_datetime('2019-05-01'),  bmm.loc['tau','PLM1'], marker='^', facecolors='white', edgecolors='tab:blue', zorder=10, label='PLM1 BMM')
        #ax2.errorbar(pd.to_datetime('2019-05-01'), bmm.loc['tau','PLM1'], yerr=abs(bmm.loc['tau','PLM1'] - np.array([bmm.loc['tau_10th','PLM1'], bmm.loc['tau_90th','PLM1']])[:,np.newaxis]), color='tab:blue', zorder=9)
        
        ax2.scatter(pd.to_datetime('2019-06-30'),  bmm.loc['tau','PLM6'], marker='o', facecolors='white', edgecolors='tab:brown', zorder=10, label='PLM6 BMM')
        ax2.errorbar(pd.to_datetime('2019-06-30'), bmm.loc['tau','PLM6'], yerr=abs(bmm.loc['tau','PLM6'] - np.array([bmm.loc['tau_10th','PLM6'], bmm.loc['tau_90th','PLM6']])[:,np.newaxis]), color='tab:brown', zorder=9)

        ax2.plot(dates_, np.array(tau_)+1000, color='black', alpha=0.0) # dummy line to maintain same axis
        
        ax2.set_ylim(999,10001)
        ax2.set_yscale('log')
        
        d = .01  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
        ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
        
        kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
        ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
                
        ax2.grid()
        ax2.margins(x=0.01)
        ax2.tick_params(which='both', bottom=False, labelbottom=False)
        ax2.spines['bottom'].set_visible(False)
        ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
        ax2.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
        for z in np.arange(2017,2021).reshape(2,2):
            ax2.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='grey')
    # Cleanup
    #ax.set_ylim(ax.get_ylim()[0]*0.9, ax.get_ylim()[1]*1.1)
    # want somewhat consistent axis
    if labs[j] in ['Base','A2']:
        ax.set_ylim(100, 999)
    elif labs[j] in ['A4','A5']:
        ax.set_ylim(10, 115)
    else:
        ax.set_ylim(20,260)
    #else:
    #    pass
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.tick_params(axis='x', rotation=20, pad=0.1)
    #ax.set_ylabel(labs[j], labelpad=1.0)
    ax.text(0.015, 0.88, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.0))
    ax.grid()
    ax.margins(x=0.01)
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='grey')
fig.text(0.055, 0.50, 'Mean Age (years)', va='center', rotation='vertical', color='black')
axes[0].legend(ncol=1, loc='upper left', bbox_to_anchor=(1.005, 1.10),
          handlelength=1.2, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13)
axes[1].legend(ncol=1, loc='upper left', bbox_to_anchor=(1.005, 1.12),
          handlelength=1.2, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13)
ax2.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.005, 1.12),
          handlelength=1.2, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13)
#fig.tight_layout()
plt.savefig(os.path.join(fdir,'fraction_tau_obs.png'), dpi=300)
plt.savefig(os.path.join(fdir,'fraction_tau_obs.svg'), format='svg')
plt.show()



# --- plot mean age as anonomolies from base ---
fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(6.5, 7.5))
fig.subplots_adjust(hspace=0.14, top=0.98, bottom=0.08, right=0.7, left=0.18)
# Get the data
tau_base = {}
for j in range(len(dirlist)):
    ax = axes[j]
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for i,w in enumerate(wells):
        tau_   = []
        dates_ = []
        for t in range(len(time_list)):
            rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
            tau_.append((rtd_df['Time'] * rtd_df['wt']).sum())
            dates_.append(time_list_dates[t])
        # plot mean age simulations
        if j == 0:
            ax.plot(dates_, tau_, color=cc[i], linestyle=ll[i], alpha=1.0, linewidth=2.0, label=names[i]if j==0 else '')
            tau_base[i] = np.array(tau_).copy()
        else:
            ax.plot(dates_, np.array(tau_)-tau_base[i], color=cc[i], linestyle=ll[i], alpha=1.0, linewidth=2.0, label=names[i]if j==0 else '') 
    # Cleanup
    #ax.set_ylim(ax.get_ylim()[0]-(abs(ax.get_ylim()[0]*0.5)), ax.get_ylim()[1]+(abs(ax.get_ylim()[1]*0.5)))
    if labs[j] in ['High K','Low por']:
        ax.set_ylim(-130, -10)
    elif labs[j] in ['Low K', 'Exp. D. K']:
        ax.set_ylim(50, 800)
    else:
        ax.set_ylim(ax.get_ylim()[0]-(abs(ax.get_ylim()[0]*0.2)), ax.get_ylim()[1]+(abs(ax.get_ylim()[1]*0.3)))
    #else:
    #    ax.set_ylim(20,200)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.tick_params(axis='x', rotation=20, pad=0.1)
    #ax.set_ylabel(labs[j], labelpad=1.0)
    ax.text(0.015, 0.88, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.0))
    ax.grid()
    ax.margins(x=0.01)
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='grey')
axes[0].set_ylabel('A1 Mean\nAge (years)')
fig.text(0.055, 0.50, 'Mean Age Anomaly from A1 (years)', va='center', ha='center', rotation='vertical', color='black')
axes[0].legend(ncol=1, loc='upper left', bbox_to_anchor=(1.005, 1.12),
          handlelength=1.2, labelspacing=0.25, columnspacing=0.5, handletextpad=0.25, fontsize=13)
plt.savefig(os.path.join(fdir,'fraction_tau_anomaly.png'), dpi=300)
plt.savefig(os.path.join(fdir,'fraction_tau_anomaly.svg'), format='svg')
plt.show()








# --- Mean Age Biases ---
for j in range(len(dirlist)):
    ax = axes[j]
    rtd_dict = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_rtd.1721.pk'))
    for i,w in enumerate(wells):
        tau_       = []
        tau_trunc  = []

        for t in range(len(time_list)):
            rtd_df, rtd_dfs = flux_wt_rtd(rtd_dict, time_list[t], w, 30)
            tau_.append((rtd_df['Time'] * rtd_df['wt']).sum())

            rtd_df_trunc = rtd_df.copy()[rtd_df.copy()['Time']<10]
            tau_trunc.append((rtd_df_trunc['Time'] * rtd_df_trunc['wt']).sum())
        tau1 = np.array(tau_)
        tau2 = np.array(tau_trunc)
        
        if w=='X508':
            print ('{} {} abs diff = {:.3f} years'.format(labs[j], w, (tau1-tau2).mean()))
            print ('{} {} %   diff = {:.3f} %'.format(labs[j], w, (100*(tau2-tau1)/tau1).mean()))












#------------------------------------------------------
#
# EcoSlim Mass Balance Plots
#
#------------------------------------------------------

_f1 = 'wy2017_2021_eco_ET_output.txt'
_f2 = 'wy2017_2021_eco_flow_output.txt'
_f3 = 'wy2017_2021_eco_PET_balance.txt'

A  = 1.5125*559*1.5125 # transect surface area

dates = pf_2_dates('2016-10-01', '2021-08-29', '24H')

def read_mass_files(_dir):
    et    = pd.read_csv(os.path.join(_dir,'ecoslim_2017_2021',_f1), delim_whitespace=True)
    flow  = pd.read_csv(os.path.join(_dir,'ecoslim_2017_2021',_f2), delim_whitespace=True)
    pet   = pd.read_csv(os.path.join(_dir,'ecoslim_2017_2021',_f3), delim_whitespace=True)
    
    et['Date']   = date_map['Date']
    flow['Date'] = date_map['Date']
    pet['Date']  = date_map['Date']
    
    et['wy']    = set_wy(et)[0]
    flow['wy']  = set_wy(flow)[0]
    pet['wy']   = set_wy(pet)[0]
    
    return et, flow, pet
    
#
# CLM ET to compare
#
# Generated using 'et_to_pickle.py'
def read_pf_out(_dir, _var):
    #pdb.set_trace()
    clm_out_ = pd.read_pickle(os.path.join(_dir,'parflow_out','clm_out_dict.2017_2021.pk')) # contains ET from CLM output files
    clm_out  = {i: v for i, v in enumerate(clm_out_.values())}  # Reset timestep keys to start with zero index
    clm_keys = list(clm_out[0].keys())
    
    df_    = pd.DataFrame((clm_out[i][_var]) for i in list(clm_out.keys())).T
    df_.columns = dates
    df = pd.DataFrame(df_)
    return df










#--------------------------------------------------------
#
# Mass-balance Plots
#
#--------------------------------------------------------

xpos = 559 - 20

#
# Subsurface Storage And Surface Storage
#
pf       = pd.read_pickle(os.path.join(dirlist[0], 'parflow_out/pf_out_dict_0021.pk'))
dates    = pf_2_dates('1999-10-01', '2021-09-30', '24H') # should be 8/29/20
_dates   = dates[dates>'2016-09-30']
_dates_inds = np.where(dates>'2016-09-30')[0]

sub_ss_base  = pf['specific_storage'].sum(axis=1).sum(axis=1)[_dates_inds[0]:]#m3/day
surf_ss_base = pf['surface_storage'].sum(axis=1)[_dates_inds[0]:] #m3/day

fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(7,8))
fig.subplots_adjust(hspace=0.12, top=0.98, bottom=0.05, right=0.85, left=0.20)
for j in range(len(dirlist)):
    pf    = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/pf_out_dict_0021.pk'))
    #sat_mask = [np.where(pf['sat'][t][:,0,:]==1, True, False) for t in _dates_inds]
    sub_ss   = pf['specific_storage'].sum(axis=1).sum(axis=1)[_dates_inds[0]:]
    #sub_ss   = pf['specific_storage'].sum(axis=1)[:,:xpos].sum(axis=1)[_dates_inds[0]:]
    ax = axes[j]
    ax.plot(_dates, sub_ss, color='C0', alpha=0.8)
    #
    surf_ss  = pf['surface_storage'].sum(axis=1)[_dates_inds[0]:]
    ax2 = ax.twinx()
    ax2.plot(_dates, surf_ss, color='C1', linestyle='--', alpha=0.8)
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(3))
    # Plot the base model for comparision
    if j != 0:
        ax.plot(_dates,  sub_ss_base, color='grey', alpha=0.5)
        ax2.plot(_dates, surf_ss_base, color='grey', linestyle='--', alpha=0.5)
    ax.grid()
    ax.margins(x=0.01)
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.set_ylabel(labs[j], labelpad=1.5)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.04, color='red')
ymax = max([axes[i].get_ylim()[1] for i in range(len(dirlist))])
ymin = min([axes[i].get_ylim()[0] for i in range(len(dirlist))])
for j in range(len(dirlist)):
    axes[j].set_ylim(ymin, ymax)
    axes[j].yaxis.set_minor_locator(ticker.AutoMinorLocator())
fig.text(0.03, 0.50, r'Subsurface Storage (m$^{3}$)', color='C0', va='center', rotation='vertical')
fig.text(0.96, 0.50,  r'Surface Storage (m$^{3}$)', color='C1', va='center', rotation='vertical')
#plt.savefig('./figures/storage_comps.png', dpi=300)
#plt.savefig(os.path.join(fdir,'storage_comps.png'), dpi=300)
plt.show()



#
# ParFlow Overland Flow
#
dates       = pf_2_dates('1999-10-01', '2021-08-29', '24H')
_dates      = dates[dates>'2016-09-30']
_dates_inds = np.where(dates > '2016-09-30')[0]

pf_base = pd.read_pickle(os.path.join(dirlist[0], 'parflow_out/pf_out_dict_0021.pk'))
Q_base  = np.array([pf_base['overland_arr'][i][0].ravel()[xpos] for i in _dates_inds])

fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(7,8))
fig.subplots_adjust(hspace=0.15, top=0.98, bottom=0.05, right=0.85, left=0.20)
for j in range(len(dirlist)):
    pf  = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/pf_out_dict_0021.pk'))
    Q   = np.array([pf['overland_arr'][i][0].ravel()[xpos] for i in _dates_inds])
    
    ax = axes[j]
    ax.plot(_dates, Q, color='C0', alpha=0.9, zorder=8, label='Runoff')
    
    if j != 0:
        ax.plot(_dates, Q_base, color='grey', linestyle='--', alpha=0.5, zorder=7)

    def m3hr_2_mmday(q):
        return q/(1.5125*559*1.5125)*1000*24
    def mmday_2_m3hr(q):
        return q*(1.5125*559*1.5125)/1000/24
    
    ax2 = ax.secondary_yaxis('right', functions=(m3hr_2_mmday, mmday_2_m3hr))
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    
    ax.grid()
    ax.margins(x=0.01)
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.set_ylabel(labs[j], labelpad=1.5)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    for z in np.arange(2017,2021).reshape(2,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.04, color='red')
ymax = max([axes[i].get_ylim()[1] for i in range(len(dirlist))])
ymin = min([axes[i].get_ylim()[0] for i in range(len(dirlist))])
for j in range(len(dirlist)):
    axes[j].set_ylim(ymin, ymax)
    axes[j].yaxis.set_minor_locator(ticker.AutoMinorLocator())
fig.text(0.05, 0.50, r'Q [m$^{3}$/hr]', color='black', va='center', rotation='vertical')
fig.text(0.93, 0.50, r'Q [mm/day]',     color='black', va='center', rotation='vertical')
#plt.savefig('./figures/overland_flow.png', dpi=300)
plt.show()
    







#
# ET, SWE, Precip., Infiltration Cumulative Curves -- mm/day
#
dates    = pf_2_dates('2016-10-01', '2021-09-30', '24H')
_clm_swe = []
_clm_et  = []
_clm_inf = []  
_overland = []
_storage  = []

area = 559*1.5125*1.5125 / 1000 # makes it units of mm/yr

for j in range(len(dirlist)):
    # SWE from CLM
    swe = pd.DataFrame([read_pf_out(dirlist[j],'swe_out').mean(axis=0)]).T
    _clm_swe.append(swe)

    # ET from CLM
    clm_et = pd.DataFrame([read_pf_out(dirlist[j],'qflx_evap_tot').loc[:xpos,:].mean(axis=0)*60*60*24]).T # mm/day - contains ET from CLM output files
    clm_et.rename(columns={0:'ET'}, inplace=True)
    clm_et['wy'] = set_wy(clm_et)[0]
    #clm_et['ET'][clm_et['ET']<0.0] = 0.0
    clm_et_cs = clm_et.groupby(by='wy').cumsum()
    _clm_et.append(clm_et_cs)

    # Infiltration from CLM
    inf = pd.DataFrame([read_pf_out(dirlist[j],'qflx_infl').loc[:xpos,:].mean(axis=0)*60*60*24]).T # mm/day - contains ET from CLM output files
    inf.rename(columns={0:'inf'}, inplace=True)
    inf['wy'] = set_wy(inf)[0]
    clm_inf_cs = inf.groupby(by='wy').cumsum()
    _clm_inf.append(clm_inf_cs)

    # Outflow from ParFlow    
    pf = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/pf_out_dict_0021.pk'))
    _dates_inds = np.where(pf_2_dates('1999-10-01', '2021-09-30', '24H') > '2016-09-30')[0]
    Q = pf['overland_arr'][:,0,:][:,xpos][_dates_inds[0]:]
    Q = pd.DataFrame(Q, index=dates)
    Q.columns = ['Q']
    Q['wy'] = set_wy(Q)[0]
    Q_cs = Q.groupby(by='wy').cumsum()
    _overland.append(Q_cs/(1.5125*559*1.5125)*1000*24)
    
    # Change is storage m3
    sub_ss   = pf['specific_storage'].sum(axis=1).sum(axis=1)[_dates_inds[0:]]   
    surf_ss  = pf['surface_storage'].sum(axis=1)[_dates_inds[0:]]   
    tot_ss   = pd.DataFrame(sub_ss+surf_ss, index=dates) #.diff()
    tot_ss['wy'] = set_wy(tot_ss)[0]
    ss_diff_yr   = [tot_ss.iloc[np.where(tot_ss['wy']==y)[0][[0,-1]], :][0].diff()[1] for y in tot_ss['wy'].unique()] # Annual change  -- m3/day
    ss_diff      = tot_ss.diff().groupby(by='wy').cumsum() # storage_day1 - storage_day0 -- m3/day
    #_storage.append(ss_diff.cumsum())
    _storage.append(pd.DataFrame(ss_diff_yr, index=tot_ss['wy'].unique()))
    
    
# plot all models together
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.0, 2.3))
fig.subplots_adjust(hspace=0.15, top=0.96, right=0.75, bottom=0.12, left=0.15)
# Total Precip from NLDAS
ax.plot(met_cs, color='black', alpha=0.8, linestyle='--', label='Precip.')
ax.fill_between(x=met_cs.index, y2=0, y1=met_cs.to_numpy().ravel(), color='grey', alpha=0.1) 
ax.fill_between(x=swe.index, y2=np.array(_clm_et)[:,:,0].min(axis=0),  y1=np.array(_clm_et)[:,:,0].max(axis=0), color='C2', alpha=0.65, label='ET')
ax.fill_between(x=swe.index, y2=np.array(_overland)[:,:,0].min(axis=0),  y1=np.array(_overland)[:,:,0].max(axis=0), color='C1', alpha=0.5, label='Runoff')
#ax.fill_between(x=swe.index, y2=np.array(_clm_swe)[:,:,0].min(axis=0), y1=np.array(_clm_swe)[:,:,0].max(axis=0), color='C0', alpha=0.65, label='CLM SWE')
#ax.fill_between(x=swe.index, y2=np.array(_clm_inf)[:,:,0].min(axis=0), y1=np.array(_clm_inf)[:,:,0].max(axis=0), color='C2', alpha=0.5, label='CLM Infil.')
# Storage
stor = np.array(_storage)[:,:,0]
_err = np.column_stack((stor.mean(axis=0)-stor.min(axis=0), stor.max(axis=0)-stor.mean(axis=0)))
_xx = [pd.to_datetime('{}-09-30'.format(z)) for z in np.arange(2017,2022)]
#ax.scatter(_xx, stor.mean(axis=0)/area, color='C2', label=r'$\Delta$Storage')
ax.errorbar(_xx, y=stor.mean(axis=0)/area, yerr=_err.T, color='C3', linestyle='None', marker='o', label=r'$\Delta$Storage')
#
ax.grid()
ax.margins(x=0.01)
ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
#ax.set_ylim(0,800)
ax.yaxis.set_major_locator(ticker.MultipleLocator(200))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
for z in np.arange(2017,2021).reshape(2,2):
    ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.04, color='red')
#ax.set_ylabel('Annual Cumulative\nFlux (mm)', labelpad=0.01)
fig.text(0.04, 0.53,  'Annual Cumulative\nFlux (mm)', color='black', ha='center', va='center', rotation='vertical')
#ax.legend(ncol=1, loc='upper center', bbox_to_anchor=(1.14, 1.02), handlelength=1.0, labelspacing=0.25, handletextpad=0.5, columnspacing=0.5, fontsize=13)
ax.legend(ncol=1, loc='upper center', bbox_to_anchor=(1.16, 1.04), handlelength=1.0, labelspacing=0.25, handletextpad=0.5, columnspacing=0.5, fontsize=13)
#plt.savefig('./figures/clm_comps.png', dpi=300)
plt.savefig(os.path.join(fdir,'clm_comps.mmday.png'), dpi=300)
plt.savefig(os.path.join(fdir,'clm_comps.mmday.svg'), format='svg')
plt.show()






#
# ET, SWE, Precip., Infiltration -- m3/day -- avoid averaging
#
dates = pf_2_dates('2016-10-01', '2021-09-30', '24H')
_clm_swe  = []
_clm_et   = []
_clm_inf  = []  
_clm_cond = []
_overland = []
_storage  = []
_resid    = []

dx, dy = 1.5125*1000, 1.5125*1000
# Put everyting in m3/day
for j in range(len(dirlist)):
    # SWE from CLM -- mm
    #swe = read_pf_out(dirlist[j],'swe_out').sum(axis=0) # total mm of SWE on the hillslope per day
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
    
    # Infiltration from CLM -- mm/s
    inf = pd.DataFrame((read_pf_out(dirlist[j],'qflx_infl').loc[:xpos,:]*dx*dy).sum(axis=0)*60*60*24)/1000**3 # m3/day
    inf.rename(columns={0:'inf'}, inplace=True)
    inf['wy'] = set_wy(inf)[0]
    clm_inf_cs = inf.groupby(by='wy').cumsum()
    _clm_inf.append(clm_inf_cs)

    # Condensation as negative ET fluxes
    clm_cond = clm_et.copy()
    clm_cond = clm_et.copy()
    clm_cond.loc[clm_cond['ET']>0.0, 'ET'] = 0.0
    clm_cond_cs = clm_cond.groupby(by='wy').cumsum()
    _clm_cond.append(clm_cond_cs)

    # Outflow from ParFlow -- m3/hr
    pf = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/pf_out_dict_0021.pk'))
    _dates_inds = np.where(pf_2_dates('1999-10-01', '2021-09-30', '24H') > '2016-09-30')[0]
    Q = pf['overland_arr'][:,0,:][:,xpos][_dates_inds[0]:]
    Q = pd.DataFrame(Q, index=dates)
    Q.columns = ['Q']
    Q *= 24 #m3/day
    Q['wy'] = set_wy(Q)[0]
    Q_cs = Q.groupby(by='wy').cumsum()
    _overland.append(Q_cs)
    #_overland.append(Q_cs/(1.5125*559*1.5125)*1000*24)
    
    # Change is storage m3
    sub_ss   = pf['specific_storage'].sum(axis=1).sum(axis=1)[_dates_inds[0:]]   
    surf_ss  = pf['surface_storage'].sum(axis=1)[_dates_inds[0:]]   
    tot_ss   = pd.DataFrame(sub_ss+surf_ss, index=dates) #.diff()
    tot_ss['wy'] = set_wy(tot_ss)[0]
    ss_diff_yr   = [tot_ss.iloc[np.where(tot_ss['wy']==y)[0][[0,-1]], :][0].diff()[1] for y in tot_ss['wy'].unique()] # Annual change  -- m3/day
    ss_diff      = tot_ss.diff().groupby(by='wy').cumsum() # storage_day1 - storage_day0 -- m3/day
    #_storage.append(ss_diff.cumsum())
    _storage.append(pd.DataFrame(ss_diff_yr, index=tot_ss['wy'].unique()))
    
    # Mass balance residuals -- P-ET-Runoff (m3/day)
    _inds = [np.where(clm_et['wy']==y)[0][-1] for y in clm_et['wy'].unique()]
    #_met_cs = met_cs.iloc[:len(Q_cs),:]*dx*dy*559/(1000**3)
    _met_cs = met_cs.iloc[:len(Q_cs),:]*dx*dy*xpos/(1000**3)
    _res = _met_cs.iloc[_inds,:].to_numpy().ravel() - clm_et_cs.iloc[_inds,:].to_numpy().ravel() - Q_cs.iloc[_inds,:].to_numpy().ravel()
    _resid.append(pd.DataFrame(_res, index=clm_et['wy'].unique()))
    
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8.0,2.7))
fig.subplots_adjust(hspace=0.15, top=0.92, bottom=0.1, right=0.78, left=0.15)
# Total Precip from NLDAS
ax.plot(met_cs*dx*dy*559/(1000**3), color='black', alpha=0.8, linestyle='--', label='Precip.')
#ax.fill_between(x=met_cs.index, y2=0, y1=met_cs.to_numpy().ravel()*dx*dy*559/(1000**3), color='grey', alpha=0.1)  
ax.fill_between(x=swe.index, y2=np.array(_clm_swe).min(axis=0),          y1=np.array(_clm_swe).max(axis=0), color='C0', alpha=0.65, label='CLM SWE')
ax.fill_between(x=swe.index, y2=np.array(_clm_et)[:,:,0].min(axis=0),    y1=np.array(_clm_et)[:,:,0].max(axis=0), color='C2', alpha=0.65, label='CLM ET')
#ax.fill_between(x=swe.index, y2=np.array(_clm_inf)[:,:,0].min(axis=0),   y1=np.array(_clm_inf)[:,:,0].max(axis=0), color='C2', alpha=0.5, label='CLM Infil.')
ax.fill_between(x=swe.index, y2=np.array(_overland)[:,:,0].min(axis=0),  y1=np.array(_overland)[:,:,0].max(axis=0), color='C1', alpha=0.5, label='Runoff')
#ax.fill_between(x=swe.index, y2=abs(np.array(_clm_cond))[:,:,0].min(axis=0), y1=abs(np.array(_clm_cond)[:,:,0]).max(axis=0), color='C4', alpha=0.65, label='CLM Cond.')
#ax.set_ylim(0,1000)
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
plt.savefig(os.path.join(fdir,'clm_comps.m3day.png'), dpi=300)
plt.show()





# --- New style - resolve all components ---
xpos = 559 - 10
dates = pf_2_dates('2016-10-01', '2021-09-30', '24H')
_dates_inds = np.where(dates>'2016-09-30')[0]
_clm_et   = []
_overland = []
_storage  = []
_resid    = []

dx, dy = 1.5125*1000, 1.5125*1000
# Put everyting in m3/day
for j in range(len(dirlist)):
    # ET from CLM -- mm/s
    #clm_et = pd.DataFrame((read_pf_out(dirlist[j],'qflx_evap_tot').loc[:xpos,:]*dx*dy).sum(axis=0)*3600*24)/1000**3 # m3/day
    clm_et = pd.DataFrame((read_pf_out(dirlist[j],'qflx_evap_tot')*dx*dy).iloc[:xpos,:].sum(axis=0)*3600*24)/1000**3 # m3/day
    clm_et.rename(columns={0:'ET'}, inplace=True)
    clm_et['wy'] = set_wy(clm_et)[0]
    #clm_et['ET'][clm_et['ET']<0.0] = 0.0
    _clm_et.append(clm_et.groupby(by='wy').sum())

    # Outflow from ParFlow -- m3/hr
    pf = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/pf_out_dict_0021.pk'))
    _dates_inds = np.where(pf_2_dates('1999-10-01', '2021-09-30', '24H') > '2016-09-30')[0]
    Q = pf['overland_arr'][:,0,:][:,xpos][_dates_inds[0]:]
    Q = pd.DataFrame(Q, index=dates)
    Q.columns = ['Q']
    Q *= 24 #m3/day
    Q['wy'] = set_wy(Q)[0]
    _overland.append(Q.groupby(by='wy').sum())
    #_overland.append(Q_cs/(1.5125*559*1.5125)*1000*24)
    
    # Change is storage m3
    sub_ss   = pf['specific_storage'].sum(axis=1).sum(axis=1)[_dates_inds[0:]]   
    #sub_ss   = pf['specific_storage'].sum(axis=1)[:,:xpos].sum(axis=1)[_dates_inds[0:]]
    surf_ss  = pf['surface_storage'].sum(axis=1)[_dates_inds[0:]]   
    tot_ss   = pd.DataFrame(sub_ss+surf_ss, index=dates)
    tot_ss['wy'] = set_wy(tot_ss)[0]
    ss_diff_yr   = [tot_ss.iloc[np.where(tot_ss['wy']==y)[0][[0,-1]], :][0].diff()[1] for y in tot_ss['wy'].unique()] # Annual change  -- m3/day
    ss_diff      = tot_ss.diff().groupby(by='wy').cumsum() # storage_day1 - storage_day0 -- m3/day
    #_storage.append(ss_diff.cumsum())
    _storage.append(pd.DataFrame(ss_diff_yr, index=tot_ss['wy'].unique()))
    
    # Mass balance residuals -- P-ET-Runoff (m3/day)
    _inds = [np.where(clm_et['wy']==y)[0][-1] for y in clm_et['wy'].unique()]
    _met_cs = met_cs*dx*dy*559/(1000**3)
    #_met_cs = met_cs.iloc[:len(Q_cs),:]*dx*dy*xpos/(1000**3)
    _res = _met_cs.iloc[_inds,:].to_numpy().ravel() - clm_et_cs.iloc[_inds,:].to_numpy().ravel() - Q_cs.iloc[_inds,:].to_numpy().ravel()
    _resid.append(pd.DataFrame(_res, index=clm_et['wy'].unique()))

et_       = np.array(_clm_et)[:,:,0]
overland_ = np.array(_overland)[:,:,0]
stor_     = np.array(_storage)[:,:,0]
resid_    = np.array(_resid)[:,:,0]
met_      = _met_cs.iloc[_inds,:].to_numpy().ravel()
_labs = dict(zip(np.arange(len(dirlist)), ['A','B','C','D','E','F']))



#labs = ['A1', 'High\nK', 'Low\nK', 'Exp. D.\nK', 'Low\npor', 'High \npor']
area = 559*1.5125*1.5125 / 1000 # makes it units of mm/yr
fig, axes = plt.subplots(nrows=et_.shape[0], ncols=et_.shape[1], figsize=(6.7, 4.0))
fig.subplots_adjust(top=0.84, bottom=0.12, left=0.08, right=0.98, hspace=0.2, wspace=0.1)
hh = 1.0
for i in range(et_.shape[1]): # the number of years as columns
    for j in range(et_.shape[0]): # the number of models as rows
        ax = axes[j,i]
        ax.axvline(met_[i]/area/10, color='C0', linestyle='--', linewidth=2.0, label='Precip.')
        ax.barh(2,  overland_[j,i]/area/10, height=hh, color='C1', label=r'Runoff')
        ax.barh(1,  et_[j,i]/area/10, height=hh, color='C2', label=r'ET')
        #ax.barh(0,  (met_[i]-overland_[j,i]-et_[j,i])/area, height=hh, color='C4', label=r'P$-$ET$-$Runoff')
        ax.barh(0, stor_[j,i]/area/10, height=hh, color='C3', label=r'$\Delta$Storage')
        if i == 0:
            ax.set_ylabel(labs[j])
        if j == 0:
            ax.set_title('WY{}'.format(str(met_daily['wy'].unique()[i])[-2:]), fontsize=13)
for i in range(axes.shape[0]):
    for j in range(axes.shape[1]):
        axes[i,j].set_ylim(-1.0, 3.5)
        #axes[i,j].set_xlim(-300, 800)
        axes[i,j].set_xlim(-25, 80)
        axes[i,j].grid(axis='x')
        
        axes[i,j].xaxis.set_major_locator(ticker.MultipleLocator(30))
        axes[i,j].xaxis.set_minor_locator(ticker.MultipleLocator(10))
        
        axes[i,j].tick_params(axis='y', labelleft=False, left=False)
        axes[i,j].axvline(0.0, color='black', linestyle='--')
        if i != axes.shape[0]-1:
            axes[i,j].tick_params(axis='x', labelbottom=False)
        else:
            axes[i,j].tick_params(axis='x', pad=1.0, rotation=0)
        axes[i,j].spines[['right', 'top']].set_visible(False)
#fig.text(0.53, 0.01, 'Annual Fluxes (m$^{3}$/yr)', ha='center')
fig.text(0.53, 0.012, 'Annual Fluxes (cm/year)', ha='center')
axes[1,0].legend(ncol=4, loc='upper center', bbox_to_anchor=(2.75, 3.95), handlelength=1.1, labelspacing=0.25, handletextpad=0.5, columnspacing=1.0, fontsize=13)
plt.savefig(os.path.join(fdir, 'clm_comps.mass_balance.sep.svg'), format='svg')
plt.savefig(os.path.join(fdir, 'clm_comps.mass_balance.sep.png'), dpi=300)
plt.show()








#------------------------------------------------------------
#
# Precipitation and SWE Timeseries
#
#------------------------------------------------------------
dates    = pf_2_dates('2016-10-01', '2021-09-30', '24H')
_prcp_sum = (met_17_21['prcp']*3600*1).groupby(pd.Grouper(freq='M')).sum() # mm/week
_prcp_sum.index = _prcp_sum.index.map(lambda x: x.replace(day=15))

_clm_swe = []
for j in range(len(dirlist)):
    # SWE from CLM
    swe = pd.DataFrame([read_pf_out(dirlist[j],'swe_out').mean(axis=0)]).T
    _clm_swe.append(swe)


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 2.5))
fig.subplots_adjust(hspace=0.15, top=0.83, bottom=0.14, right=0.83, left=0.16)
# Total Precip from NLDAS
ax.plot(met_cs, color='black', alpha=0.8, linestyle='--', label='NLDAS Precip.')
# SWE
ax.plot(swe.index, np.array(_clm_swe)[:,:,0].mean(axis=0), color='C0', alpha=1.0, label='CLM SWE')
# Monthly Precip
ndays_month = [(_prcp_sum.index[j+1]-_prcp_sum.index[j]).days for j in range(len(_prcp_sum)-1)] + [30]
ax2 = ax.twinx()
ax2.bar(x=_prcp_sum.index, height=_prcp_sum, width=ndays_month, color='C1', alpha=1.0)
ax2.invert_yaxis()

ax.set_ylim(0,1100)
ax.set_yticks(ticks=np.arange(0,801,200), labels=np.arange(0,801,200))
ax.set_yticks(ticks=np.arange(0,801,50), minor=True)
ax.set_ylabel('Cumulative\nFlux (mm)')
ax.yaxis.set_label_coords(-0.12, 0.45)

ax2.set_ylim(400,0)
ax2.set_yticks(ticks=np.arange(0,201,50), labels=np.arange(0,201,50))
ax2.set_yticks(ticks=np.arange(0,201,25), minor=True)
ax2.set_ylabel('Precip.\n(mm/month)')
ax2.yaxis.set_label_coords(1.13, 0.72)

#ax.grid()
ax.margins(x=0.01)
ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10]))
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
#ax.set_ylim(0,1000)
#ax.yaxis.set_major_locator(ticker.MultipleLocator(200))
#ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
for z in np.arange(2017,2021).reshape(2,2):
    ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.04, color='red')
#fig.text(0.93, 0.50,  r'Q [mm/day]', color='black', va='center', rotation='vertical')
ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.26), handlelength=1.0, labelspacing=0.25, handletextpad=0.5, columnspacing=0.5, fontsize=13, frameon=False)
#plt.savefig('./figures/NLDAS_Precip.png', dpi=300)
plt.savefig(os.path.join(fdir, 'NLDAS_Precip.png'), dpi=300)
plt.show()












#-----------------------------------------------------------------------
#
# ET Age Dynamics
#
#-----------------------------------------------------------------------
land_surf = pd.read_csv('elevation_v4.sa', skiprows=1, header=None)
xx  = np.arange(0,559*1.5125,1.5125)
xx_ = np.arange(0,559)

fig, ax = plt.subplots(figsize=(2*3.17, 2))
ax.plot(xx, land_surf)
ax.grid()
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.minorticks_on()
ax.margins(x=0.01)
ax.set_xlabel('Distance (X)')
ax.set_ylabel('Elevation (m)')
fig.tight_layout()
plt.show()



et_age     = pd.DataFrame()
et_age_top = pd.DataFrame()
et_age_mid = pd.DataFrame()
et_age_bot = pd.DataFrame()

et_mass     = pd.DataFrame()
et_mass_top = pd.DataFrame()
et_mass_mid = pd.DataFrame()
et_mass_bot = pd.DataFrame()

et_xloc = pd.DataFrame()

for j in range(len(dirlist)):
    _et_age = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_ETAge.1721.pk'))
    for t in list(_et_age.keys())[1:]:
        _ett = _et_age[t]
        # --- hillslope-scale ---
        _wt    = _ett['mass'] / _ett['mass'].sum()
        et_age.loc[t,j]  = ((_ett['age'] * _wt).sum() / 8760)
        et_mass.loc[t,j] = _ett['mass'].mean()
        # --- top of hillslope ---
        msk  = np.where(_ett['x']<150, True, False)
        _wt  = _ett['mass'][msk] / _ett['mass'][msk].sum()
        et_age_top.loc[t,j]  = ((_ett['age'][msk] * _wt).sum() / 8760)
        et_mass_top.loc[t,j] = _ett['mass'][msk].mean()
        # --- middle of hillslope ---
        msk  = np.where((_ett['x']>150)&(_ett['x']<400), True, False)
        _wt  = _ett['mass'][msk] / _ett['mass'][msk].sum()
        et_age_mid.loc[t,j]  = ((_ett['age'][msk] * _wt).sum() / 8760)
        et_mass_mid.loc[t,j] = _ett['mass'][msk].mean()
        # --- bottom of hillslope ---
        msk  = np.where((_ett['x']>400)&(_ett['x']<600), True, False)
        _wt  = _ett['mass'][msk] / _ett['mass'][msk].sum()
        et_age_bot.loc[t,j]  = ((_ett['age'][msk] * _wt).sum() / 8760)
        et_mass_bot.loc[t,j] = _ett['mass'][msk].mean()

# cleanup
def build_wy_df(df):
    df.index = dates
    df['wy'] = set_wy(df)[0]
    df_ = df.groupby(by='wy').mean()
    return df_

et_age_wy  = build_wy_df(et_age)
et_mass_wy = build_wy_df(et_mass)
# --
et_age_top_wy  = build_wy_df(et_age_top)
et_mass_top_wy = build_wy_df(et_mass_top)
# --
et_age_mid_wy  = build_wy_df(et_age_mid)
et_mass_mid_wy = build_wy_df(et_mass_mid)
# --
et_age_bot_wy  = build_wy_df(et_age_bot)
et_mass_bot_wy = build_wy_df(et_mass_bot)                         


# -- Plot --
month_map  = {10:'O',11:'N',12:'D',1:'J',2:'F',3:'M',4:'A',5:'M',6:'J',7:'J',8:'A',9:'S'}
colors = plt.cm.twilight(np.linspace(0,1,len(month_map)+1))


fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10,5))
fig.subplots_adjust(top=0.96, bottom=0.15, left=0.15, right=0.8, hspace=0.18, wspace=0.32)
for j in range(len(dirlist)):
    _et_age  = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_ETAge.1721.pk'))
    _et_dict = pd.concat(_et_age)
    #_et_dict['age_wt'] = (_et_dict['age']*_et_dict['mass']/_et_dict['mass'].sum())
        
    _et_dict['wy']    = [date_map.loc[m,'wy'] for m in _et_dict.index.get_level_values(0)-1]
    _et_dict['year']  = [date_map.loc[m,'Date'].year for m in _et_dict.index.get_level_values(0)-1]
    _et_dict['month'] = [date_map.loc[m,'Date'].month for m in _et_dict.index.get_level_values(0)-1]

    # trim out lower part of domain?
    _et_dict = _et_dict[_et_dict['x']<540]
    # pick single water year?
    _et_dict = _et_dict[_et_dict['wy']==2019]
    

    r,c = j//3, j%3
    ax = axes[r,c]
    ax.scatter(_et_dict['x']*1.5125, _et_dict['age']/8760, marker='.', color='grey', alpha=0.5)#c=_et_dict['month'])
    for i,z in enumerate(_et_dict['month'].unique()):
        _msk = np.where(_et_dict['month'] == z, True, False)
        _wt  = _et_dict['mass'][_msk] / _et_dict['mass'][_msk].sum()
        ax.scatter((_et_dict['x'][_msk]*_wt).sum()*1.5125, (_et_dict['age'][_msk]*_wt).sum()/8760, color=colors[i], label=month_map[z])
        
    ax.minorticks_on()
    ax.grid()
    ax.text(0.05, 0.9, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    #ax.set_yscale('log')
fig.text(0.5, 0.04, 'ET Hillslope Position (m)', ha='center')
fig.text(0.02, 0.55, 'ET Mean Age (years)', va='center', rotation='vertical')
axes[0,2].legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), fontsize=12.5, handlelength=1, labelspacing=0.25)
plt.show()



# -- Timeseries plot -- 
fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(6.5, 7.5))
fig.subplots_adjust(hspace=0.14, top=0.98, bottom=0.08, right=0.95, left=0.15)
for j in range(len(dirlist)):
    _et_age  = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/ecoslim_ETAge.1721.pk'))
    _et_dict = pd.concat(_et_age)
    #_et_dict['age_wt'] = (_et_dict['age']*_et_dict['mass']/_et_dict['mass'].sum())
        
    _et_dict['wy']    = [date_map.loc[m,'wy'] for m in _et_dict.index.get_level_values(0)-1]
    _et_dict['date']  = [date_map.loc[m,'Date'] for m in _et_dict.index.get_level_values(0)-1]

    # trim out lower part of domain?
    _et_dict = _et_dict[_et_dict['x']<540]
    
    ax = axes[j]

    _age_flux_wt = []
    _date_map = date_map[np.where((date_map.wy==2017)|(date_map.wy==2018)|(date_map.wy==2019),True,False)]
    for i in _date_map.index:
        try:
            _age_flux_wt.append((_et_dict.loc[i+1,:]['age']*_et_dict.loc[i+1,:]['mass']/_et_dict.loc[i+1,:]['mass'].sum()).sum()/8760)
        except KeyError:
            _age_flux_wt.append(np.NaN)
    _age_flux_wt = np.array(_age_flux_wt)
   
    
    ax.plot(_date_map['Date'], _age_flux_wt)
        
    ax.grid()
    #ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[10,2,6]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.tick_params(axis='x', rotation=20, pad=0.1)
    #ax.set_ylabel(labs[j], labelpad=1.0)
    ax.text(0.015, 0.88, labs[j], horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.0))
    ax.margins(x=0.01)
    for z in np.arange(2017,2019).reshape(1,2):
        ax.axvspan(pd.to_datetime('{}-10-01'.format(z[0])), pd.to_datetime('{}-09-30'.format(z[1])), alpha=0.06, color='red')
fig.text(0.02, 0.55, 'ET Mean Age (years)', va='center', rotation='vertical')
#axes[0,2].legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), fontsize=12.5, handlelength=1, labelspacing=0.25)
plt.show()







#------------------------------------------------------------
#
# Water Tables Dyamics
#
#------------------------------------------------------------



land_surf = pd.read_csv('elevation_v4.sa', skiprows=1, header=None)
xx = np.arange(0,559*1.5125,1.5125)
xx_ = np.arange(0,559)




# Using the Parflow pftools scipt
_dir = 1

pf = pd.read_pickle(os.path.join(dirlist[_dir], 'parflow_out/pf_out_dict_0021.pk'))
wtd = pd.DataFrame.from_dict(pf['wtd'])
wtd.index = pf_2_dates('1999-10-01', '2021-09-30', '24H')
wtd = wtd[wtd.index>'2016-09-30']
wy_ = np.array(set_wy(wtd)[0])

#
# Plot Water Table Depths for single model and all years
#
#month_map  = {1:'O',2:'N',3:'D',4:'J',5:'F',6:'M',7:'A',8:'M',9:'J',10:'J',11:'A',12:'S'}
#months     = list(month_map.values())
month_map  = {10:'O',11:'N',12:'D',1:'J',2:'F',3:'M',4:'A',5:'M',6:'J',7:'J',8:'A',9:'S'}
colors = plt.cm.twilight(np.linspace(0,1,len(month_map)+1))


fig, axes = plt.subplots(nrows=len(np.unique(wy_)), ncols=1, figsize=(5, 6))
fig.subplots_adjust(top=0.94, bottom=0.12, left=0.2, right=0.78, hspace=0.1)
for y in range(len(np.unique(wy_))):
    ax = axes[y]
    wy = np.unique(wy_)[y]
    wt = wtd[np.where(wy_==wy, True, False)]
    # Plot all timesteps
    ax.fill_between(x=xx, y1=wt.min(axis=0), y2=wt.max(axis=0), color='grey', alpha=0.4)
    # Plot the first of the months
    #first_month = np.where(wt.index.day==1)[0]
    #for j in range(len(first_month)):
    #    jj = first_month[j]
    #    ax.plot(xx, wt.iloc[jj,:], color=colors[j], alpha=1.0, label='{}'.format(months[j]))
    # Plot average of each month
    for j in range(len(month_map)):
        mm = list(month_map.keys())[j]
        ax.plot(xx, wt[np.where(wt.index.month==mm,True,False)].mean(axis=0), color=colors[j], alpha=1.0, label='{}'.format(month_map[mm]))
    #ax.text(0.1, 0.85, wy_[y], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.minorticks_on()
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.set_ylim(wtd.to_numpy().min()-1, wtd.to_numpy().max()+2)
    ax.invert_yaxis()
    ax.margins(x=0.01)
    ax.set_ylabel(np.unique(wy_)[y])
    if y != len(axes)-1:
        ax.tick_params(axis='x', labelbottom=False)
    ax.grid(axis='y')
    ax.axvline(404*1.5125, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(494*1.5125, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(508*1.5125, color='grey', linestyle='--', alpha=0.5)
axes[0].set_title(labs[_dir], fontsize=13)
axes[0].legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), fontsize=12.5, handlelength=1, labelspacing=0.25)
fig.text(0.03, 0.5, 'Water Table Depth (mbls)', va='center', rotation='vertical')
axes[len(np.unique(wy_))-1].set_xlabel('Distance (m)')
plt.savefig(os.path.join(fdir, 'wt_temp_{}.png'.format(''.join(labs[_dir].split()))),dpi=300)
plt.show()



#
# Plot Water Table Depths for all models and single years
#
wy = 2019
fig, axes = plt.subplots(nrows=len(dirlist), ncols=1, figsize=(5, 6.5))
fig.subplots_adjust(top=0.94, bottom=0.12, left=0.22, right=0.78, hspace=0.14)
for j in range(len(dirlist)):
    pf = pd.read_pickle(os.path.join(dirlist[j], 'parflow_out/pf_out_dict_0021.pk'))
    wtd = pd.DataFrame.from_dict(pf['wtd'])
    wtd.index = pf_2_dates('1999-10-01', '2021-09-30', '24H')
    wtd = wtd[wtd.index>'2016-09-30']
    wy_ = np.array(set_wy(wtd)[0])

    ax = axes[j]
    wt = wtd[np.where(wy_==wy, True, False)]
    # Plot all timesteps
    ax.fill_between(x=xx, y1=wt.min(axis=0), y2=wt.max(axis=0), color='grey', alpha=0.5)
    # Plot the first of the months
    #first_month = np.where(wt.index.day==1)[0]
    #for i in range(len(first_month)):
    #    jj = first_month[i]
    #    ax.plot(xx, wt.iloc[jj,:], color=colors[i], alpha=1.0, label='{}'.format(months[i]))
    # Plot average of each month
    for jj in range(len(month_map)):
        mm = list(month_map.keys())[jj]
        ax.plot(xx, wt[np.where(wt.index.month==mm,True,False)].mean(axis=0), color=colors[jj], alpha=1.0, label='{}'.format(month_map[mm]))
    #ax.text(0.1, 0.85, wy_[y], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.minorticks_on()
    ax.set_ylim(wtd.to_numpy().min()-1, wtd.to_numpy().max()+2)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    if wtd.max().max() < 20.0:
        ax.set_ylim(-1, 20.0)
    else:
        ax.set_ylim(-1, 55.0)
    ax.invert_yaxis()
    ax.margins(x=0.01)
    ax.set_ylabel(labs[j])
    if j != len(dirlist)-1:
        ax.tick_params(axis='x', labelbottom=False)
    ax.grid(axis='y')
    ax.axvline(404*1.5125, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(494*1.5125, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(508*1.5125, color='grey', linestyle='--', alpha=0.5)
axes[0].set_title('WY{}'.format(wy), fontsize=13)
axes[0].legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), fontsize=12.5, handlelength=1, labelspacing=0.25)
fig.text(0.02, 0.5, 'Water Table Depth (mbls)', va='center', rotation='vertical')
axes[len(dirlist)-1].set_xlabel('Distance (m)')
plt.savefig(os.path.join(fdir,'wt_spatial.{}.png'.format(wy)), dpi=300)
plt.savefig(os.path.join(fdir,'wt_spatial.{}.svg'.format(wy)), format='svg')
plt.show()


























