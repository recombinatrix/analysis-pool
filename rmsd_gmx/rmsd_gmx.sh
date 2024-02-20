#!/bin/bash

rmsd() {

# runs rmsd calculations using gromacs tools
# expecys that your input files have a consistent naming scheme
# requires gro, xtc and ndx files 
# works faster with stripped / dried trajectories

# my memory is that gromacs won't let you easilly fit to one thing, and calculate RMSD of a second thing
# when you want to calculate ligand RMSD, you need prepare a preexisting trajectory that has been fit to the protein backbone
# my memory is that this allows you to calculate ligand RMSD fit to the protein using nofit
# however looking at this script, I seem not to have done that.  Perhaps that only applied to ligand clustering?


# expected folders and files

# ./calcs/                      - this is where the RMSD will be saved, as an .xvg
# ./calcs/concat.py             - combines all the .xvg files into a single csv, extracting metadata from the filename 

    # example index groups used in this script:

    #  4: backbone
    # 13: ligand (in ligand bearing systems)
    # 14: N-bundle all atom
    # 15: N-bundle backbone
    # 16: C-bundle all atom
    # 17: C-bundle backbone
    # 18: ECD all atom
    # 19: ECD backbone
    # 20: TMD all atom
    # 21: TMD backbone
    # 22: TMD backbone plus ligand
    # 23: TMD all atom plus ligand
    # 24: Entire protein plus ligand


for sys in ${systems[@]}; do
    # calc protein and ligand rmsd for all systems
    # requires $systems to be a envirnmental variable as an array.  I do this because passing arrays to bash functions is a pain

    # for all rmsd calcs, structure is fit to backbone of either entire protein or of a subdomain

    # entire protein RMSD fit to self
    printf "4\n4\n"   | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_fit_prot_calc_prot_BB.xvg 
    
    # TMD RMSD fit to self
    printf "21\n21\n" | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_fit_TMD_calc_TMD_BB.xvg 

    # ECD RMSD fit to TMD backbone
    printf "21\n19\n" | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_fit_TMD_calc_ECD_BB.xvg

    # ECD rmsd fit to ECD backbone

    printf "19\n19\n" | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_fit_ECD_calc_ECD_BB.xvg


    if [ $sys != "hOCT1_apo_rebuild" ]; then # calc ligand rmsd relative to TMD backbone only in ligand bearing systems
    # doing it this way required me to padd the index group files for apo systems with empty lines, so every index file had the same items at the same number
    # this is a pain but it worked

        if [ $sys != "hOCT1-ASP357ASPH_apo_D357prot_stripped" ]; then 

        # ligand RMSD relative to starting frame, fit to backbone
        printf "21\n13\n" | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../clean/stripped/${sys}_stripped.gro -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_fit_TMD-start_calc_LIGAND_ALL.xvg

        # ligand RMSD relative to cluster structure, fit to backbone

        printf "21\n13\n" | gmx rms -f ../clean/stripped/${sys}_stripped.xtc -s ../cluster/pdb/${sys}_TMDBBlig_clusters.pdb -n ../clean/stripped/${sys}_stripped.ndx -xvg none -o calcs/${sys}_fit_TMD-clust-1_calc_LIGAND_ALL.xvg

        # cluster structure pdbs contain more than one structure, gromacs will only look at the first structure.  If I want to look at any secondary clusters I will need to do those calculations manually
        fi
    fi
    done
}



systems=( "hOCT1_MF1B_p3" "hOCT1_MF1A_p2" "hOCT1_MF1A_p3" "hOCT1_MF1B_p6" "hOCT1_MF1B_p6-clust-switch" "hOCT1-ASP357ASPH_MF1B_p3" "hOCT1-ASP357ASPH_MF1A_p2" "hOCT1-ASP357ASPH_MF1A_p3" "hOCT1-ASP357ASPH_MF1B_p6" )
rmsd

# systems=( "hOCT1-ASP357ASPH_FNRR_p3" "hOCT1-ASP357ASPH_FNRR_p7" "hOCT1-ASP357ASPH_FNSS_p3" "hOCT1-ASP357ASPH_FNSS_p4" "hOCT1_FNRR_p3" "hOCT1_FNRR_p7" "hOCT1_FNSS_p3" "hOCT1_FNSS_p4" )
# systems=("hOCT1_MF1B_p6" )
# rmsd

# systems=( "hOCT1_THIA_p1" "hOCT1_THIA_p3-2" "hOCT1_THII_p1" "hOCT1_THII_p2" "hOCT1-ASP357ASPH_THIA_p1" "hOCT1-ASP357ASPH_THIA_p3-2" "hOCT1-ASP357ASPH_THII_p1" "hOCT1-ASP357ASPH_THII_p2")
# systems=( "hOCT1-ASP357ASPH_THIA_p3-2" )

# rmsd

# systems=( "hOCT1_PPNS_p4" "hOCT1_PPNR_p2-3" "hOCT1_PPNS_pCRYO" "hOCT1-ASP357ASPH_PPNS_pCRYO" )
# rmsd

# systems=( "hOCT1_apo_rebuild" "hOCT1-ASP357ASPH_apo_D357prot"  )
# rmsd

# systems=( "hOCT1_DTZ_p2-1" "hOCT1_DTZ_p2-2" "hOCT1-ASP357ASPH_DTZ_p2-1" 
# systems=("hOCT1-ASP357ASPH_DTZ_p2-2"  )
# rmsd

# systems=("hOCT1_FNSS_p4r1ex"  )
systems=("hOCT1_FNRR_p7r3ex"  )
rmsd

# combine all rmsd data into a single csv

cd calcs
./concat.py
cd ..

# # testing:

# # gmx rms -f ../clean/stripped/hOCT1_FNRR_p7_stripped_fit.xtc -n ../clean/stripped/hOCT1_FNRR_p7_stripped.ndx -s ../cluster/pdb/hOCT1_FNRR_p7_TMDBBlig_clusters.pdb -xvg none -o clustertest.xvg
# # looks like it ran a comparison against only the first cluster in the record


# gmx rms -f ../clean/stripped/hOCT1_FNRR_p7_stripped.xtc -n ../clean/stripped/hOCT1_FNRR_p7_stripped.ndx -s ../cluster/pdb/hOCT1_FNRR_p7_TMDBBlig_clusters.pdb -xvg none -o clustertest.xvg
