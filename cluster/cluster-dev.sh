#!/bin/bash


# taken directly from prod analysis
# this preobably needs to be checked and refactored, and commented

# trajectory clustering for all systems
# The goal is to identify the cluster structures of the TMD, the ligand aligned to the TMD backbone, and the entire TMD/ligand complex
# had to rerun these as my original clustering only saved the coodinates of the ligand and TMD, which meant I couldn't use them for RMSD calcs.

# index groups

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

clust_trj() {
        for sys in ${systems[@]}; do

            # fit stripped trj to TMD BB so I can look at ligand RMSD in isolation
        printf "21\n0\n" | gmx trjconv -f ../clean/stripped/${sys}_stripped.xtc -n ../clean/stripped/${sys}_stripped.ndx -s ../clean/stripped/${sys}_stripped.gro -o ../clean/stripped/${sys}_stripped_fit.xtc -fit rot+trans 2>&1 | tee trjlog/${sys}_trjlog.txt
        
        done
}

clusters() {

    for sys in ${systems[@]}; do

    # if I want to do them, entire protein clusters would look like this
    # at this stage I'm not going to bother as I expect ECD movement will dominate these


    # printf "4\n1\n" | gmx cluster -f ../clean/stripped/${sys}_stripped_fit.xtc -n ../clean/stripped/${sys}_stripped.ndx -s ../clean/stripped/${sys}_stripped.gro -g logs/${sys}_prot.log -cl pdb/${sys}_prot_clusters.pdb -cutoff 0.2 -dist dist/${sys}_prot_rmsd-dist.xvg -o xpm/${sys}_prot_rmsd-clust.xpm -method gromos

    # TMD clusters

    printf "20\n20\n" | gmx cluster -f ../clean/stripped/${sys}_stripped_fit.xtc -n ../clean/stripped/${sys}_stripped.ndx -s ../clean/stripped/${sys}_stripped.gro -g logs/${sys}_TMD.log -cl pdb/${sys}_TMD_clusters.pdb -cutoff 0.2 -dist dist/${sys}_TMD_rmsd-dist.xvg -o xpm/${sys}_TMD_rmsd-clust.xpm -method gromos

        if [ $sys != "hOCT1_apo_rebuild" ]; then 
            if [ $sys != "hOCT1-ASP357ASPH_apo_D357prot" ]; then 

                # calculate ligand aligned to TMD backbone
                # the TMD backbone has already been fit to itself, hence use -nofit to look at changes to ligand RMSD within the cavity
                # this should give clusters based on ligand rmsd, showing the distribution of ligand poses
        
                printf "13\n0\n" | gmx cluster -f ../clean/stripped/${sys}_stripped_fit.xtc -nofit -n ../clean/stripped/${sys}_stripped.ndx -s ../clean/stripped/${sys}_stripped.gro -g logs/${sys}_TMDBBlig.log -cl pdb/${sys}_TMDBBlig_clusters.pdb -dist dist/${sys}_TMDBBlig_rmsd-dist.xvg -o xpm/${sys}_TMDBBlig_rmsd-clust.xpm -method gromos -cutoff 0.2 -clid clid/${sys}_TMDBBlig_membership.xvg

                # calculate TMD all atom + ligand aligned to TMD backbone
                # the ligand RMSD dominates as the TMD backbone has already been fit to itself, but side chain fluctuation can also play a role.  This might be useful for getting constraints for refining the cryoEM model, but it may be very noisy.  A future step might be just looking at residues from the core 

               printf "23\n0\n" | gmx cluster -f ../clean/stripped/${sys}_stripped_fit.xtc -nofit -n ../clean/stripped/${sys}_stripped.ndx -s ../clean/stripped/${sys}_stripped.gro -g logs/${sys}_TMDlig.log -cl pdb/${sys}_TMDlig_clusters.pdb -cutoff 0.2 -dist dist/${sys}_TMDlig_rmsd-dist.xvg -o xpm/${sys}_TMDlig_rmsd-clust.xpm -method gromos 

            fi
        fi

    done

}

systems=( "hOCT1_MF1B_p3" "hOCT1_MF1A_p2" "hOCT1_MF1A_p3" "hOCT1_MF1B_p6" )

clust_trj
clusters

