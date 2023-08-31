lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

pfset     FileVersion    4


# example usage w/ command line
# "tclsh file_conversion.sa2pfb.tcl elevation.sa"

# Command line argument needs sa file name
set safile $argv

# Change to .pfb file
set sp [split $safile .]
set base [lindex $sp 0]
set pfb ".pfb"
set pfbfile $base$pfb

puts "working on $safile --> $pfbfile"

#Converting sa to pfb and silo
set ind [pfload -sa $safile]
pfsave $ind -pfb $pfbfile

