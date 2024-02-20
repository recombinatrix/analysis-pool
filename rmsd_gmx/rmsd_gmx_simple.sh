#!/bin/bash

rmsd() {

# runs rmsd calculations using gromacs tools
# expecys that your input files have a consistent naming scheme
# requires gro, xtc and ndx files 
# works faster with stripped / dried trajectories

# expected folders and files

# ./calcs/                      - this is where the RMSD will be saved, as an .xvg

for sys in ${systems[@]}; do

    printf "4\n4\n"   | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_BB.xvg 

    done
}



# save a list of systems by their common name

systems=( "hOCT1_MF1B_p3" "hOCT1_MF1A_p2" "hOCT1_MF1A_p3" "hOCT1_MF1B_p6" )  # example names used in my analysis

# run rmsd on everything in your list

rmsd
