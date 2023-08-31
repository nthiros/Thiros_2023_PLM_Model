#!/bin/bash

# Note: the sed function below works differently on MacOS versus Linux. This script is designed for linux on big computeres


#conda activate parflow

wdir=$(pwd)

ss=64   # starting model number
ee=127  # ending model number


for i in $(seq -f "%05g" $ss $ee)
do
  
  echo "Running $i of $ee"
  rdir="Run"$i

  cp "./pull_parflow_wl_v6.py" "./et_to_pickle.py" "parflow_post_proc.py" $rdir

  # Jump into rundir
  cd $rdir

  # Pull waterlevel outputs
  python pull_parflow_wl_v6.py  
  python et_to_pickle.py
  python parflow_post_proc.py

  # Back to parent dir
  cd $wdir 
done
