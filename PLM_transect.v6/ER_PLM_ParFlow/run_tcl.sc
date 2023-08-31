#!/bin/bash

# set directory names, these must match what is set in the tcl scripts
d_spin='wy_spinup'
d00_16='wy_2000_2016'
d17_21='wy_2017_2021'


#tcl_p='/home/nt126396/software/parflowTree/tcl/local/bin/tclsh8.6'
#tcl_p='/home/nt126396/Software/parflowTree/tcl/local/bin/tclsh8.6'
 


#
# Run the spinup
echo 'Running Spinup'
#$tcl_p $(pwd)/transect_spinup_trm.tcl
tclsh $(pwd)/transect_spinup_trm.tcl


#
# Sparcify 2000 to 2016 to 20 day outputs
echo 'Sparcifying...'
python3 sparse_vel_files.py


#
# Run 2000 to 2016
echo 'Running 2000-2016...'
mkdir $d00_16
cp $d_spin/clm.rst.00000.* $d00_16
cp $d_spin/clm_restart.tcl $d00_16
#cp $d_spin/$d_spin.out.press.03650.pfb $d80_16
#$tcl_p $(pwd)/transect_2000-2016_trm.tcl
tclsh $(pwd)/transect_2000-2016_trm.tcl


#
# Run 2017 to 2021
echo 'Running 2017-2021...'
mkdir $d17_21
cp $d00_16/clm.rst.00000.* $d17_21
cp $d00_16/clm_restart.tcl $d17_21
#cp $d80_16/$d80_16.out.press.13505.pfb $d17_21
#$tcl_p $(pwd)/transect_2017-2021_trm.tcl
tclsh $(pwd)/transect_2017-2021_trm.tcl

