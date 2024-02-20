#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
import pandas as pd

from MDAnalysis.analysis import distances

from tqdm import tqdm  # progress bar

top = 'path/to/topology.gro' # path to .gro, .tpr or .pdb
trj = 'path/to/trajectory.xtc' # path to .xtc or .trr
name = 'protein_ligand_contacts'  # the name of this calculation.  will be used to generate filenames
prot_str = 'protein' # selection string for protein
lig_str = 'resname Org Org2'  # selection string for ligand, in this case residues named Org and Org2

# the two selection strings should not overlap

# variables used in calc, you can set these as you require

cutoff = 4 # cutoff distance in A
highfreq = 0.8 # what % of frames constitutes a high frequency contact 

# load the universe

u = mda.Universe(top,trj)

# define a function to calculate contacts
# by default, returns c_df, a dataframe listing all contact events per frame, and count, a dataframe showing the per-frame frequency of contact by reidue pairs

def contacts(u,prot_str,lig_str,
        cutoff=4,  # distances in A <= cuttoff constitute a contact
        savedist=False,  # set to true if you want to save all pairwise ligand/residue distances calculated for every frame.  With large selections this is very slow.  this produces third output, the per frame minimum distance between each residue pair 
        savecount=True, # do you want to save the count as a csv, named 
        savecrude=True, # do you want to save the raw pairwise per-frame contacts as a csv?  This will be very large.  If savedist is true, this will also save the pairwise distances between all residue paris per frame.  That file will be gigantic.   
        skip0=False,  # set to true if you want to exclude the starting frame from the calculations
        highfreq=0.8,  # percentage, generates a dataframe showing only high frequency interactions.  can be set to false to skip that step.  if savecount is true, saves this as a csv
        name='contacts',  # what is the name of the calculation?  this is used in filenames
        rounding=2, # how many digits are distances rounded to?  Only really relevant for savedist 
        ):

    ligand = u.select_atoms(lig_str)
    protein = u.select_atoms(prot_str)

    out = "_".join(name.split(" "))

    #check if ligand and protein selections are overlapping, fail if they are

    if not set(protein.atoms).isdisjoint(set(ligand.atoms)):
        print('ligand and protein selections are overlapping.  Exiting without performing calculations') 
        print(f'protein selection: {prot_str}') 
        print(f'ligand selection: {lig_str}') 
        return (None, None)

    # I don't like using ligand and protein, I want terms for "thing contacting & thing being contacted"
    # ligand and protein will do for now

    pdict = {} # dictionary linking resids to resnames
    # this may not work correctly if you have multiple residues with the same resid
    # you probably should not be doing that that anyway 

    for res in protein.residues:
        pdict[int(res.resid)]=str(res.resname)
    for lig in ligand.residues:
        pdict[int(lig.resid)]=str(lig.resname)

    c_touch = [] # timestep, prot resid,ligand resid, only for entries where dist <= cutoff
    c_dist = [] # 

    trjlen = len(u.trajectory) - int(skip0) # ignore initial frame if skip0 is true

    for ts in tqdm(u.trajectory[int(skip0):len(u.trajectory)], desc="Calculating contacts across trajectory"):  
        # right now this starts at frame 1 to ignore initial coordinates; 
        # I think in the future I should start processing trajectories to dump the initial frame as a gro, 
        # and not include it in the xtc, otherwise it gets triple counted in concatenation

        for res in protein.residues:
            for lig in ligand.residues:
                dist = round(distances.distance_array(res.atoms.positions, lig.atoms.positions).min(),rounding)  
                if dist <= cutoff: 
                    c_touch.append([ts.frame,res.resid,lig.resid])

                # this next lines are used to generate a list of pairwise distances between the ligand and each residue, for every frame 
                # with large selections (eg an entire protein), this is very slow and produces a gigantic dataset
                # needs optimising, only use it if you need it, or know what you are doing
                
                if savedist:
                    here = [ts.frame,res.resid,lig.resid,dist]             
                    c_dist.append(here)
                    np.array(c_dist)

    c_df = pd.DataFrame(c_touch,columns=['frame','prot resid','lig resid'])
    if savecrude: c_df.to_csv(f"{out}_crude_contacts.csv")


    count = c_df.groupby(['prot resid','lig resid']).count()
    count['frac'] = count['frame']/trjlen
    count['prot resname'] = [pdict[i[0]] for i in count.index]
    count['lig resname'] = [pdict[i[1]] for i in count.index]
    count.set_index(['prot resname','lig resname'], append=True,inplace=True)
    count = count.reorder_levels(['lig resid','lig resname','prot resid','prot resname',])
    count.drop('frame',axis=1,inplace=True)
    if savecount: count.reset_index().to_csv(f"{out}.csv",index=False)
    if highfreq: 
        print(count.loc[count.frac >= highfreq])
        if savecount: count.loc[count.frac >= highfreq].reset_index().to_csv(f"{out}_highfreq.csv",index=False)

    if savedist:
        d_df = pd.DataFrame(c_dist,columns=['frame','prot resid','lig resid','min dist'])
        if savecrude: d_df.to_csv(f"{out}_crude_distances.csv")
        return (c_df,count,d_df)

    else: return(c_df,count)

# now, run the contacts

c_df, count = contacts(u,prot_str=prot_str,lig_str=lig_str,
        cutoff=cutoff,
        savedist=False,
        savecount=True,
        savecrude=True,
        skip0=True,
        name=name,
        highfreq=highfreq)

# done


