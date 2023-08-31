#!/bin/bash

# Note: the sed function below works differently on MacOS versus Linux. This script is designed for linux on big computeres


source activate parflow

wdir=$(pwd)

ss=0   # starting model number
ee=15  # ending model number

ncores=12

outlog=runlog_${ss}_${ee}.txt

echo $(date) >> $outlog
echo "$ncores cores" >> $outlog

for i in $(seq -f "%05g" $ss $ee)
do
  start=`date +%s`  

  rdir="Run"$i
  mkdir $rdir
  echo "Running $i of $ee"

  # Start copying files into new dir
  cp "./Training_Set/subsurface_$i.pfb" $rdir
  cp "./Training_Set/subsurface_$i.mat" $rdir

  cp "./pull_parflow_wl_v6.py" "./et_to_pickle.py" $rdir
  cp "run_tcl.sc" "transect_spinup.tcl" "transect_2012-2021.tcl" $rdir
  cp "slope_x_v4.pfb" "slope_y_v4.pfb" $rdir
  cp "wy_spinupA.out.press.00730.pfb" $rdir  

  # Jump into rundir
  cd $rdir

  # Update ParFlow tcl with correct run number 
  sed -i "s/subsurface_00000/subsurface_$i/g" transect_spinup.tcl 
  sed -i "s/subsurface_00000/subsurface_$i/g" transect_2012-2021.tcl 
  
  sed -i "s/Process.Topology.P 24/Process.Topology.P $ncores/g" transect_spinup.tcl 
  sed -i "s/Process.Topology.P 24/Process.Topology.P $ncores/g" transect_2012-2021.tcl 

  # Run ParFlow
  ./run_tcl.sc

  # Pull waterlevel outputs
  python pull_parflow_wl_v6.py  
  python et_to_pickle.py

  # Runtimes
  end=`date +%s`
  runtime=$((end-start))
  echo "  Runtime=$runtime s" 

  # Back to parent dir
  cd $wdir 
  
  # Save time runtimes
  echo "Run $i = $runtime s" >> $outlog
done
