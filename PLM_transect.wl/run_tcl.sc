#!/bin/bash

# set directory names, these must match what is set in the tcl scripts
d_spin='wy_spinup'
d17_21='wy_2017_2021'


#tcl_p='/home/nt126396/software/parflowTree/tcl/local/bin/tclsh8.6'
#tcl_p='/home/nt126396/Software/parflowTree/tcl/local/bin/tclsh8.6'
 
#
# Run the spinup
echo 'Running Spinup'
tclsh $(pwd)/transect_spinup_trm.tcl


#
# Run 2017 to 2021
echo 'Running 2017-2021...'
mkdir $d17_21
cp $d_spin/clm.rst.00000.* $d17_21
cp $d_spin/clm_restart.tcl $d17_21
#$tcl_p $(pwd)/transect_2017-2021_trm.tcl
tclsh $(pwd)/transect_2017-2021_trm.tcl

