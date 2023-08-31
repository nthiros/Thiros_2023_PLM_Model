#!/bin/bash

export OMP_NUM_THREADS=24

#---
# Set file files
# and directories

# spinup -- longtime run
dspn1='ecoslim_spinup_lng'
mkdir $dspn1
cp slimin.txt.spinup_lng $dspn1/slimin.txt
dspn1_=$(awk 'NR==1{print $1}' $dspn1/slimin.txt)


# spinup
dspn='ecoslim_spinup'
mkdir $dspn
cp slimin.txt.spinup $dspn/slimin.txt
dspn_=$(awk 'NR==1{print $1}' $dspn/slimin.txt)


# 2000-2016
d00_16='./ecoslim_2000_2016'
mkdir $d00_16
cp slimin.txt.2000_2016 $d00_16/slimin.txt
d00_16_=$(awk 'NR==1{print $1}' $d00_16/slimin.txt)


# 2017-2021
d17_21='./ecoslim_2017_2021'
mkdir $d17_21
cp slimin.txt.2017_2021 $d17_21/slimin.txt
d17_21_=$(awk 'NR==1{print $1}' $d17_21/slimin.txt)


f1='_exited_particles.bin'
f2='_particle_restart.bin'


echo 'Ecoslim'

#---
# spinup1 -- longtime
#---
cd $dspn1
echo 'Working on long spinup...'
cp ../../elevation_v4.pfb .
EcoSLIM.exe
cd ..


#---
# spinup
#---
wdir_o=$dspn1
fnam_0=$dspn1_

wdir_c=$dspn
fnam_1=$dspn_

cd $wdir_c
cp ../../elevation_v4.pfb .
cp ../$wdir_o/*.bin .
mv $fnam_0$f1 $fnam_1$f1
mv $fnam_0$f2 $fnam_1$f2
echo 'Working on spinup...'
EcoSLIM.exe
cd ..


#---
# 2000-2016
#---
wdir_o=$dspn
fnam_0=$dspn_

wdir_c=$d00_16
fnam_1=$d00_16_

cd $wdir_c
cp ../../elevation_v4.pfb .
cp ../$wdir_o/*.bin .
mv $fnam_0$f1 $fnam_1$f1
mv $fnam_0$f2 $fnam_1$f2
echo 'Working on 2000-2016...'
EcoSLIM.exe
cd ..


#---
# 2017-2021
#---
wdir_o=$d00_16
fnam_0=$d00_16_

wdir_c=$d17_21
fnam_1=$d17_21_

cd $wdir_c
cp ../../elevation_v4.pfb .
cp ../$wdir_o/*.bin .
mv $fnam_0$f1 $fnam_1$f1
mv $fnam_0$f2 $fnam_1$f2
echo 'Working on 2017-2021...'
EcoSLIM.exe
cd ..

