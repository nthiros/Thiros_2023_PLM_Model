#!/bin/bash

start_time="$(date -u +%s)"

./run_tcl.sc
./run_eco.sc

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Runtime: $elapsed seconds"
