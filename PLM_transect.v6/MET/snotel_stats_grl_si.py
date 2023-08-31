import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pdb

import matplotlib.dates as mdates
import matplotlib.ticker as ticker

plt.rcParams['font.size'] = 14



def set_wy(df):
    dates     = df.copy().index
    yrs       = dates.year
    yrs_      = np.unique(yrs)[1:]
    wy_inds_  = [np.where((dates > '{}-09-30'.format(i-1)) & (dates < '{}-10-01'.format(i)), True, False) for i in yrs_]
    wy_inds   = np.array([wy_inds_[i]*yrs_[i] for i in range(len(yrs_))]).sum(axis=0)
    first_yrs = [(wy_inds==i).argmax() for i in yrs_]
    return list(wy_inds), list(first_yrs)


def pf_2_dates(startdate, enddate, f):
    '''Assumes ParFlow outputs every 24 hours'''
    s = pd.to_datetime(startdate)
    e = pd.to_datetime(enddate)
    d_list = pd.date_range(start=s, end=e, freq=f)
    # Drop Leap years again
    d_list_ = d_list[~((d_list.month == 2) & (d_list.day == 29))]
    return d_list_



#---------------------------------
#
# Import and clean SNOTEL data
#
#--------------------------------- 
# Column Labels
#    Date,
#    Snow Water Equivalent (in) Start of Day Values,
#    Precipitation Accumulation (in) Start of Day Values,
#    Air Temperature Maximum (degF),
#    Air Temperature Minimum (degF),
#    Air Temperature Average (degF),
#    Precipitation Increment (in)


fname = './Butte_snotel.txt'

header = ['Date','SWE','P_accum','T_max','T_min','T_avg','P_inc']
dds    = pd.read_csv(fname, delimiter=',', skiprows=64, header=None, names=header, index_col='Date', parse_dates=True)

dds.loc[:,['SWE','P_accum','P_inc']] *= 25.4 # convert SWE from inches to mm


dds.loc[:,['T_max','T_min','T_avg']] = (dds.loc[:,['T_max','T_min','T_avg']]-32)*5/9 # temperature in C

# begin to interpolate any Nans
#dds = dds.resample('1D').interpolate('quadratic',limit_direction='both', pad=True)

dds['wy'] = set_wy(dds)[0] # add WY column

t1 = '2016-10-01'
t2 = '2021-09-30'
_dds = dds[t1:t2]


dds_2000 = dds[dds['wy']==2000]
first_month = np.array([np.where((dds_2000.index.month==m) & (dds_2000.index.day==1))[0][0] for m in [10,11,12,1,2,3,4,5,6,7,8,9]])
month_labs = ['O','N','D','J','F','M','A','M','J','J','A','S']


# Plot SWE curves
fig, ax = plt.subplots(figsize=(6,4))
for i,y in enumerate(dds['wy'].unique()):
    _dd = dds[dds['wy']==y]
    ax.plot(np.arange(len(_dd)), _dd['SWE'], color='grey', alpha=0.3)
    
for i,y in enumerate(_dds['wy'].unique()):
    _dd = _dds[_dds['wy']==y]
    ax.plot(np.arange(len(_dd)), _dd['SWE'], color='C{}'.format(i), alpha=1.0, label='WY{}'.format(y))

# --- Cleanup ---
ax.set_xticks(ticks=first_month, labels=month_labs)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.grid(alpha=0.2)
ax.margins(x=0.02)
ax.set_ylabel('SWE (mm)')
ax.legend(handlelength=1.0, labelspacing=0.25, handletextpad=0.5)
fig.tight_layout()
#plt.savefig('Butte_SWE_curves.png', dpi=300)
plt.show()



# Plot cumulative precip curves
fig, ax = plt.subplots(figsize=(6,4))
for i,y in enumerate(dds['wy'].unique()):
    _dd = dds[dds['wy']==y]
    ax.plot(np.arange(len(_dd)), _dd['P_accum'], color='grey', alpha=0.3)
    
for i,y in enumerate(_dds['wy'].unique()):
    _dd = _dds[_dds['wy']==y]
    ax.plot(np.arange(len(_dd)), _dd['P_accum'], color='C{}'.format(i), alpha=1.0, label='WY{}'.format(y))

# --- Cleanup ---
ax.set_xticks(ticks=first_month, labels=month_labs)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.grid(alpha=0.2)
ax.margins(x=0.02)
ax.set_ylabel('Precip. (mm/year)')
ax.legend(handlelength=1.0, labelspacing=0.25, handletextpad=0.5)
fig.tight_layout()
#plt.savefig('Butte_SWE_curves.png', dpi=300)
plt.show()



fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(6,6))
fig.subplots_adjust(top=0.98, bottom=0.1, left=0.15, right=0.96)
for i,y in enumerate(dds['wy'].unique()):
    _dd = dds[dds['wy']==y]
    ax[0].plot(np.arange(len(_dd)), _dd['SWE'], color='grey', alpha=0.25)
    ax[1].plot(np.arange(len(_dd)), _dd['P_accum'], color='grey', alpha=0.25)
    
for i,y in enumerate(_dds['wy'].unique()):
    _dd = _dds[_dds['wy']==y]
    ax[0].plot(np.arange(len(_dd)), _dd['SWE'], color='C{}'.format(i), linewidth=2.0, alpha=1.0, label='WY{}'.format(y))
    ax[1].plot(np.arange(len(_dd)), _dd['P_accum'], color='C{}'.format(i), linewidth=2.0, alpha=1.0, label='WY{}'.format(y))

# --- Cleanup ---
for i in [0,1]:
    ax[i].set_xticks(ticks=first_month, labels=month_labs)
    ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[i].grid(alpha=0.2)
    ax[i].margins(x=0.02)
ax[0].set_ylabel('SWE (mm)')
ax[1].set_ylabel('Precip. (mm/year)')
ax[1].set_xlabel('Month')
ax[0].text(0.035, 0.9, '(A)', fontweight='bold', horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.035, 0.9, '(B)', fontweight='bold', horizontalalignment='left', verticalalignment='center', transform=ax[1].transAxes)
ax[0].legend(handlelength=1.0, labelspacing=0.25, handletextpad=0.5)
plt.savefig('SI_Butte_SWE_curves.png', dpi=300)
plt.show()




