#!/bin/bash


# runs rmsf calculations for a protein using gromacs tools
# expecys that your input files have a consistent naming scheme
# requires gro, xtc and ndx files 

    # expected folders and files

    # ./calcs/                      - this is where the RMSF will be saved, as an .xvg
    # ./calcs/concat.py             - combines all the .xvg files into a single csv, extracting metadata from the filename 
    # ./calcs_subdomain/            - this is where the RMSF for subdomains will be saved as an .xvg
    # ./calcs_subdomain/concat.py   - combines all the .xvg files into a single csv, extracting metadata from the filename 
    # ./pdb/                        - this is where structures labed with rmsf in the bfactor column will be saved 

# first define the function

rmsf() {

for sys in ${systems[@]}; do
    # calc protein and ligand rmsf for all systems in the $systems variable
    # requires $systems to be a envirnmental variable as an array.  I do this because passing arrays to bash functions is a pain

    # ./calcs               - this is where the RMSF will be saved, as an .xvg
    # ./calcs_subdomain     - this is where the RMSF for subdomains will be saved as an .xvg
    # ./pdb                 - this is where structures labed with rmsf in the bfactor column will be saved 

    # entire protein RMSF relative to self
    printf "1\n"   | gmx rmsf -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_RMSF.xvg -oq ${sys}_ALL_RMSF_bf.pdb -res  

    # ======== SUDOMAIN SECTION START ========    

    # do it again with two subdomains fit to themselves
    # this is useful if you have a protein with multiple domains moving independantly, where relative motion perturbs the fit
    # for example, the hOCT1 transmembrane domain is quite static, but the extracellular domain was extremely flexible
    # the ECD was large enough I wanted to exclude relative domain/domain motion from the fit & calculation
    # in this example, we're using index group 20 (hOCT1 TMD, all atom) and group 18 (hOCT1 ECD, all atom)

    printf "20\n"  | gmx rmsf -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs_subdomain/${sys}_TMD_RMSF.xvg -oq ${sys}_TMD_RMSF_bf.pdb -res  
    printf "18\n"  | gmx rmsf -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs_subdomain/${sys}_ECD_RMSF.xvg -oq ${sys}_ECD_RMSF_bf.pdb -res  

    cat ${sys}_TMD_RMSF_bf.pdb ${sys}_ECD_RMSF_bf.pdb > ${sys}_TMDfit_and_ECDfit_bf.pdb # stitch the self fit TMD.pdb and the self fit ECD.pdb together into one file

    # ======== SUDOMAIN SECTION END ========    

    mv *.pdb pdb/ # for some reason generating the file directly in the subfolder is failing, move it into a subdomain called pdb
    done

}

# define which systems will be run, using their standard name

systems=( "hOCT1_DTZ_p2-1" "hOCT1_DTZ_p2-2" "hOCT1-ASP357ASPH_DTZ_p2-1"  )

# run rmsf

rmsf

# clean up data by combining all individual calcs into a single csv

cd calcs
./concat.py  
cd ../calcs_subdomain
./concat.py  
