#!/bin/bash

cp ../ER_PLM_ParFlow/pull_parflow_wl_v6.py .
cp ../ER_PLM_ParFlow/read_vtk.110222.py .
cp ../ER_PLM_ParFlow/et_to_pickle.py .
cp ../ER_PLM_ParFlow/parflow_post_proc.py .


python pull_parflow_wl_v6.py
python read_vtk.110222.py
python et_to_pickle.py
#python sat_press_to_pickle.py
python parflow_post_proc.py
