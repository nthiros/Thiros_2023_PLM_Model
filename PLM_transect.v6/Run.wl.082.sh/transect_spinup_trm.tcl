# This is a spinup
# Forcing with long-term average singal of 1979 to 2021 MET data (original NLDAS files hourly transience)
# Total simulation time is 10 years
# Using 1hr timesteps

#---------------------------------------------------------
# Setup Rundir and Restarts
#---------------------------------------------------------
set runname wy_spinup
file  mkdir "./$runname"
cd "./$runname"

# This is for CLM -- time zero is istep=1, restart anything greater than 1
set istep 1

# This will reset istep variable
if { $istep > 1 } then {
  source clm_restart.tcl
}

#---------------------------------------------------------
# Import the ParFlow TCL package
#---------------------------------------------------------
set tcl_precision 17

lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*

#---------------------------------------------------------
# File input version number
#---------------------------------------------------------
pfset FileVersion 4

#---------------------------------------------------------
# Process Topology --
#---------------------------------------------------------
pfset Process.Topology.P        24
pfset Process.Topology.Q        1
pfset Process.Topology.R        1


set met_fname "met.avg.10yr.1hr.txt"
set ind_fname "subsurface_00082.pfb"
set mat_fname "subsurface_00082.mat"
set ip        "../wy_spinupA.out.press.00730.pfb"
if { $istep == 1 } then {
        file copy -force "../../slope_x_v4.pfb" .
	file copy -force "../../slope_y_v4.pfb" .
	file copy -force "../$ind_fname" .
	file copy -force "../../MET/$met_fname" .
	file copy -force "../$mat_fname" .
        file copy -force "../../clm_inputs/drv_clmin_shrub.dat" "./drv_clmin.dat"
        file copy -force "../../clm_inputs/drv_vegp.dat"  .
        file copy -force "../../clm_inputs/drv_vegm_shrub.dat"  .
}

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX               559
pfset ComputationalGrid.NY               1
pfset ComputationalGrid.NZ               32

pfset ComputationalGrid.DX		 1.5125
pfset ComputationalGrid.DY		 1.5125
pfset ComputationalGrid.DZ		 10.0

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
pfset GeomInput.Names                 	"domain_input indi_input"

pfset GeomInput.domain_input.GeomName  		domain
pfset GeomInput.domain_input.InputType  	Box

#---------------------------------------------------------
# Domain Geometry 
#---------------------------------------------------------
pfset Geom.domain.Lower.X                        0.0
pfset Geom.domain.Lower.Y                        0.0
pfset Geom.domain.Lower.Z                        0.0
 
pfset Geom.domain.Upper.X                        845.4875
pfset Geom.domain.Upper.Y                        1.5125
pfset Geom.domain.Upper.Z                        320.0
pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"



#-----------------------------------------------------------------------------
# Simulation Time
#-----------------------------------------------------------------------------
# Units are hours
# 10 year simulation
set runlength [expr (10*8760)]

# Start counter=0 for spinup, or $istep if restarting in same directory
set startcount 0

# This will extend previous run to total of runlength hours. 
# ie. runlength is end time, not added time, thus end of previous run must be < runlength
set starttime 0.0
set stoptime $runlength

puts -----
puts startcount:$startcount
puts istep:$istep
puts starttime:$starttime
puts endtime:$stoptime
puts -----

pfset TimingInfo.BaseUnit        1.0
pfset TimingInfo.StartCount      $startcount
pfset TimingInfo.StartTime       $starttime
pfset TimingInfo.StopTime        $stoptime
pfset TimingInfo.DumpInterval    24
pfset TimeStep.Type              Constant
pfset TimeStep.Value             1.0



### START .mat file
source [pwd]/$mat_fname

#---------------------------------------------------------
# variable dz assignments
#---------------------------------------------------------

#---------------------------------------------------------
# Indicator Geometry Input
#---------------------------------------------------------

#-----------------------------------------------------------------------------
# Permeability 
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Porosity 
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------
### END .mat file


#-----------------------------------------------------------------------------
# Permeability Tensors
#-----------------------------------------------------------------------------
pfset Perm.TensorType                   TensorByGeom
pfset Geom.Perm.TensorByGeom.Names      "domain"
pfset Geom.domain.Perm.TensorValX       1.0d0
pfset Geom.domain.Perm.TensorValY       1.0d0
pfset Geom.domain.Perm.TensorValZ       1.0d0


#-----------------------------------------------------------------------------
# Specific Storag
#-----------------------------------------------------------------------------
pfset SpecificStorage.Type			Constant
pfset SpecificStorage.GeomNames			"domain"
pfset Geom.domain.SpecificStorage.Value 	1.0e-5


#-----------------------------------------------------------------------------
# Phases  -- Setting to 1.0 allows for the calculation of K instead of k
#-----------------------------------------------------------------------------
pfset Phase.Names				"water"
pfset Phase.water.Density.Type	      	 	 Constant
pfset Phase.water.Density.Value	       		 1.0

pfset Phase.water.Viscosity.Type		Constant
pfset Phase.water.Viscosity.Value		1.0


#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names		""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------
pfset Geom.Retardation.GeomNames	""

pfset Gravity				1.0


#-----------------------------------------------------------------------------
#Set CLM Restart
#-----------------------------------------------------------------------------
if { $istep > 1 } then {
 exec cp ../clm_inputs/drv_clmin.restart.dat ./drv_clmin.dat

set  nproc [expr [pfget Process.Topology.P]*[pfget Process.Topology.Q]*[pfget Process.Topology.R]]

for { set i 0 } { $i < $nproc } { incr i 1 } {
      set fname_rst [format "clm.rst.%05d.$i" [expr $istep]]
      exec cp clm.rst.00000.$i $fname_rst
    }
}

#-----------------------------------------------------------------------------
# Domain -- Report/Upload the domain problem that is defined above
#-----------------------------------------------------------------------------
pfset Domain.GeomName			"domain"

#-----------------------------------------------------------------------------
# Mobility -- Mobility between phases -- Only one phase in this problem
#-----------------------------------------------------------------------------
pfset Phase.water.Mobility.Type		Constant
pfset Phase.water.Mobility.Value	1.0

#-----------------------------------------------------------------------------
# Wells
#----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names 				"constant"
pfset Cycle.constant.Names			"alltime"
pfset Cycle.constant.alltime.Length    		  1
pfset Cycle.constant.Repeat           		  -1
 
#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-upper.BCPressure.Type                  DirEquilRefPatch
pfset Patch.x-upper.BCPressure.Cycle                 "constant"
pfset Patch.x-upper.BCPressure.RefGeom                domain
pfset Patch.x-upper.BCPressure.RefPatch               z-upper
pfset Patch.x-upper.BCPressure.alltime.Value         0.443843365

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"  
pfset Patch.z-lower.BCPressure.alltime.Value	       0.0 

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

pfset Patch.z-upper.BCPressure.Type                   OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle                  "constant"
pfset Patch.z-upper.BCPressure.alltime.Value          0.0

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------
pfset TopoSlopesX.Type 		"PFBFile"
pfset TopoSlopesX.GeomNames 	"domain"
pfset TopoSlopesX.FileName	slope_x_v4.pfb

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------
pfset TopoSlopesY.Type 		"PFBFile"
pfset TopoSlopesY.GeomNames	"domain"
pfset TopoSlopesY.FileName  	slope_y_v4.pfb

#---------------------------------------------------------
# Mannings Coefficient
#---------------------------------------------------------
pfset Mannings.Type		 "Constant"
pfset Mannings.GeomNames         "domain"
pfset Mannings.Geom.domain.Value  5.52e-6


#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------
pfset PhaseSources.water.Type				Constant
pfset PhaseSources.water.GeomNames			domain
pfset PhaseSources.water.Geom.domain.Value		0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
pfset KnownSolution                                    NoKnownSolution

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------
pfset Solver						Richards
pfset Solver.MaxIter					250000

pfset Solver.TerrainFollowingGrid			True

pfset Solver.Nonlinear.MaxIter				400
pfset Solver.Nonlinear.ResidualTol			1e-5
pfset Solver.Nonlinear.EtaValue				0.001

pfset Solver.PrintSubsurf				False
pfset Solver.Drop					1E-15
pfset Solver.AbsTol					1E-10
pfset Solver.MaxConvergenceFailures			7

pfset Solver.Nonlinear.UseJacobian			True 
pfset Solver.Nonlinear.StepTol				1e-19
pfset Solver.Nonlinear.AbsTol				1e-5

pfset Solver.Linear.MaxIter				15000
pfset Solver.Nonlinear.Globalization			LineSearch
pfset Solver.Linear.KrylovDimension			90
pfset Solver.Linear.MaxRestarts				6

pfset Solver.Linear.Preconditioner			PFMG
pfset Solver.Linear.Preconditioner.PCMatrixType		FullJacobian

#Using kinematic wave
#pfset OverlandFlowDiffusive				1

#------------------
# CLM
#------------------
# Output setup
pfset Solver.LSM                                      CLM
pfset Solver.CLM.CLMFileDir                           "clm_outputs"
pfset Solver.CLM.Print1dOut                           False
pfset Solver.BinaryOutDir                             False
pfset Solver.CLM.CLMDumpInterval                      24

# Evap and veg stress functions/parameters
pfset Solver.CLM.EvapBeta                             Linear
pfset Solver.CLM.VegWaterStress                       Pressure 
pfset Solver.CLM.ResSat                               0.08
pfset Solver.CLM.WiltingPoint                         -150.0
pfset Solver.CLM.FieldCapacity                        0.0

# LW WAS:
#pfset Solver.CLM.VegWaterStress                       Saturation
#pfset Solver.CLM.ResSat                               0.1
#pfset Solver.CLM.WiltingPoint                         0.12
#pfset Solver.CLM.FieldCapacity                        0.98

# Met forcing and timestep setup
pfset Solver.CLM.MetForcing                           1D
pfset Solver.CLM.MetFileName			      $met_fname 
pfset Solver.CLM.MetFilePath			      "./"

pfset Solver.CLM.IstepStart                           $istep

# Irrigation setup
#pfset Solver.CLM.IrrigationTypes                      none

# Writing output:
pfset Solver.PrintSubsurfData                         True
pfset Solver.PrintPressure                            True
pfset Solver.PrintSaturation                          True

pfset Solver.WriteCLMBinary                           False
pfset Solver.WriteSiloCLM			      False
pfset Solver.PrintCLM				      True
pfset Solver.CLM.ReuseCount				1
pfset Solver.CLM.WriteLogs				True
pfset Solver.CLM.WriteLastRST				True
pfset Solver.CLM.DailyRST				True
pfset Solver.CLM.SingleFile				True
pfset Solver.CLM.RootZoneNZ                           4
pfset Solver.CLM.SoiLayer                             4

pfset Solver.WriteSiloSubsurfData		      False
pfset Solver.WriteSiloPressure			      False
pfset Solver.WriteSiloVelocities                      False
pfset Solver.WriteSiloSaturation                      False

# SLIM keys required
pfset Solver.PrintEvapTrans                           True
pfset Solver.PrintVelocities                          True

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
if {$istep == 1 } {
  pfset ICPressure.Type                                   PFBFile
  pfset ICPressure.GeomNames                              domain

  pfset Geom.domain.ICPressure.FileName                   $ip
  pfdist						  $ip
  }

if {$istep > 1 } {
  set outnum [expr (int($istep*1/24))]

  pfset ICPressure.Type                                   PFBFile
  pfset ICPressure.GeomNames                              domain
  #set ip [format $runname.out.press.%05d.pfb $startcount]
  set ip [format $runname.out.press.%05d.pfb $outnum]
  pfset Geom.domain.ICPressure.FileName                   $ip
  pfdist                                                  $ip
}


#--------------------------------------------------------
# Distribute Inputs
#--------------------------------------------------------
pfset ComputationalGrid.NZ              1
pfdist slope_x_v4.pfb
pfdist slope_y_v4.pfb
pfset ComputationalGrid.NZ              32

pfdist $ind_fname

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfrun      $runname
pfwritedb  $runname
pfundist   $runname
pfundist   $ind_fname
pfundist   "../$ip"

cd "../"
