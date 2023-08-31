lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

pfset     FileVersion    4


# example usage w/ command line
# "tclsh file_conversion.sa2pfb.tcl elevation.sa"

# Command line argument needs sa file name
set infile $argv

# Change to .pfb file
set sp [split $infile .]
set base [lindex $sp 0]
set sa ".sa"
set safile $base$sa

puts "working on $infile --> $safile"

#Converting sa to pfb and silo
set ind [pfload -pfb $infile]
pfsave $ind -sa $safile

