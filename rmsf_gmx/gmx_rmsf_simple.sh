#!/bin/bash


# runs rmsf calculations for a protein using gromacs tools
# expecys that your input files have a consistent naming scheme
# requires gro, xtc and ndx files 

    # expected folders and files

    # ./calcs/                      - this is where the RMSF will be saved, as an .xvg
    # ./pdb/                        - this is where structures labed with rmsf in the bfactor column will be saved 

# first define the function

rmsf() {

for sys in ${systems[@]}; do
    # calc protein and ligand rmsf for all systems in the $systems variable
    # requires $systems to be a envirnmental variable as an array.  I do this because passing arrays to bash functions is a pain

    # ./calcs               - this is where the RMSF will be saved, as an .xvg
    # ./pdb                 - this is where structures labed with rmsf in the bfactor column will be saved 

    # entire protein RMSF relative to self
    printf "1\n"   | gmx rmsf -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_RMSF.xvg -oq ${sys}_ALL_RMSF_bf.pdb -res  


    mv *.pdb pdb/ # for some reason generating the file directly in the subfolder is failing, move it into a subdomain called pdb
    done

}

# define which systems will be run, using their standard name

systems=( "hOCT1_DTZ_p2-1" "hOCT1_DTZ_p2-2" "hOCT1-ASP357ASPH_DTZ_p2-1"  ) # example system names

# run rmsf

rmsf
