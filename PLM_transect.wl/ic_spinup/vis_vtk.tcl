#-----------------------------------------------------------------------------
# Properly visualize terrain following grid
# and variable dz_layers
#-----------------------------------------------------------------------------

### Import the ParFlow TCL package
lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

pfset     FileVersion    4	

### set the parflow run dir
set runname       "wy_ic"
set outtime       "00000"


### load in parflow outputs
# dem
set dem [pfload -pfb elevation.pfb]
# permeability field
set perm [pfload -pfb "$runname/$runname.out.perm_z.pfb"]
# a saturation field
set sat [pfload -pfb "$runname/$runname.out.satur.$outtime.pfb"]
# a pressure field
set prs [pfload -pfb "$runname/$runname.out.press.$outtime.pfb"]


### set the variable dz layer thicknesses (not the multipliers
set dzlist "32\
            10.0 10.0 10.0 10.0 10.0 8.0 8.0 6.0 6.0\
            4.0 4.0 2.0 2.0 2.0 1.0 1.0 1.0 0.5\
            0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\
            0.5 0.5 0.5 0.25 0.25"


### Generate new vtk that can be visualized in Paraview
pfvtksave $perm -vtk "tfg.out.Perm.vtk" -var "perm" -dem $dem -tfg $dzlist
pfvtksave $sat -vtk "tfg.out.Sat.$outtime.vtk" -var "sat" -dem $dem -tfg $dzlist
pfvtksave $prs -vtk "tfg.out.Press.$outtime.vtk" -var "press" -dem $dem -tfg $dzlist
