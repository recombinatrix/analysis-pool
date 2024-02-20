#!/usr/bin/env python                                                       
# coding: utf-8 

import numpy as np                                                              
import pandas as pd                                                             


import argparse
import os

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms


site_dict = {
    "None" : "None",
    "S1" : "S1",
    "upper" : "LAS", 
    "lower" : "LAS",
    "flex"  : "LAS",
    "phenyl"  : "LAS",
    "antiphenyl" : "LAS",
}

i_dict = {
    # dictionary of selection strings for inhibitors
    "E73" : "resname EN73" ,
    "E25R" : "resname E25R" ,
    "E25S" : "resname E25S" ,

}


# this is a mess of complexity demon nonsense and needs refactoring
# it's very old and I know I've made a simpler one

parser = argparse.ArgumentParser()
parser.add_argument("protein", type=str, help="protein name")
parser.add_argument("membrane", type=str, help="POPC_CHOL")
parser.add_argument("inhibitor", type=str, help="None E25R E25S E73")
parser.add_argument("pose", type=str, help="None upper lower flex phenyl antiphenyl S1")

args = parser.parse_args()

sys_str = "_".join([args.inhibitor, args.pose, args.protein , args.membrane ,""])
print(sys_str)

flist = os.listdir("clean/")
syslist = [f for f in flist if sys_str in f]
syslist = [f for f in syslist if not "offset" in f]

f0 = ["clean/" + f for f in syslist if "gro" in f]
ftrj = ["clean/" + f for f in syslist if "xtc" in f]
rep_d = {}
for i in sorted([ [x[0], x[1]] for x in zip (f0,ftrj)], key=lambda x:x[0]): # the sorting puts the reps in ascending order
    rep_d[i[0].split("_")[5]] = i 

print(rep_d)

def load_sim_rep(rep_d, rep):
    u = mda.Universe(rep_d[i][0], rep_d[i][1])
    return(u)

def load_sim_cat(rep_d):
    trjlist = []
    for i in rep_d.keys():
        trjlist.append(rep_d[i][1])
        coord = (rep_d[i][0])
    u = mda.Universe( coord, trjlist )
    return(u)


def rmsd_aligner(rep_d, force=False):
    fname = ( "rmsd/trj/" + sys_str + 'aligned_self.xtc')
    u = load_sim_cat(rep_d)
    if (os.path.exists(fname) == False) or (force == True) :
        print('run alignment: ' + fname)
        aligner = align.AlignTraj(u, u, select='protein', filename=fname, verbose=True)
        aligner.run()
        ('alignment complete')
    else:
        print(fname + " already exists, use force=True to force calculation of a new alignment")
    return(fname)

def run_rmsd_rep(rep_d,force=False):
    trj_aln = rmsd_aligner(rep_d, force=False)
    for i in rep_d.keys():
        coord = (rep_d[i][0])
    u_aligned, ref = mda.Universe(coord, trj_aln) , mda.Universe(coord)
    u_aligned.trajectory[-1]  # set mobile trajectory to last frame
    ref.trajectory[0]  # set reference trajectory to first frame
    u_protein = u_aligned.select_atoms('protein')
    ref_protein = ref.select_atoms('protein')
    print("running protein rmsd for " + sys_str)
    protein_RMSD = rms.RMSD(u_aligned,  # universe to align PROTEIN 
             u_aligned,  # reference universe or atomgroup
             select='backbone',  # group to superimpose and calculate RMSD
             groupselections=['protein'],  # groups for RMSD
             ref_frame=0, verbose=True)   # frame index of the reference
    protein_RMSD.run()

    df = pd.DataFrame(protein_RMSD.rmsd,
                    columns=['Frame', 'Time (ps)',
                            'backbone','protein rmsd'])
    if args.inhibitor in i_dict:
        print("running inhibitor rmsd for system: " + sys_str)
        i_selstr = i_dict[args.inhibitor]
        u_i = u_aligned.select_atoms( i_selstr )
        ref_i = ref.select_atoms( i_selstr )     
        i_RMSD  = rms.RMSD(u_aligned,  # universe to align SUBSTRATE                    
             u_aligned,  # reference universe or atomgroup                  
             select= 'backbone or ' + i_selstr,  # group to superimpose and calculate RMSD  
             groupselections=[ i_selstr ],  # groups for RMSD                
             ref_frame=0, verbose=True)  # frame index of the reference                   
        i_RMSD.run()                                                          
                                                                            
        i_df = pd.DataFrame(i_RMSD.rmsd,                              
                  columns=['Frame', 'Time (ps)',                            
                           'ref','inhibitor rmsd'])          
        df['inhibitor rmsd'] = i_df['inhibitor rmsd']
    else: print("Did not run RMSD for inhibitor\"" + args.inhibitor + "\", not defined in the inhibitors dictionary")

    df["protein"] = args.protein
    df["membrane"] = args.membrane
    df["inhibitor"] = args.inhibitor
    df["pose"] = args.pose
    df["site"] = site_dict[args.pose]
    if args.pose == "phenyl":
        df["replicate"] = np.repeat([1,2,3,4,5], 501)
    elif (args.pose == "S1") and (args.inhibitor == "E25R"):
        df["replicate"] = np.repeat([2,3], 501)
    else:
        df["replicate"] = np.repeat([1,2,3], 501)


    return(df)

def summarise_rmsd(df):
    if args.pose == "phenyl":
        reps = [1,2,3,4,5,"all"]
    elif (args.pose == "S1") and (args.inhibitor == "E25R"):
        reps = [2,3,'all']
    else:
        reps = [1,2,3,"all"]

    for rep in reps:
        if rep != "all":
            r_df = df.loc[df["replicate"] == rep]
            fname = "rmsd/summary/" + sys_str + "r" + str(rep) + "_rmsd_summary.csv"
        else:
            r_df = df
            fname = "rmsd/summary/" + sys_str + "cat1-3_rmsd_summary.csv"
        summary_df = pd.DataFrame()
        summary_df["protein"] = [args.protein]
        summary_df["membrane"] = [args.membrane]
        summary_df["inhibitor"] = [args.inhibitor]
        summary_df["pose"] = [args.pose]
        summary_df["site"] = [site_dict[args.pose]]
        summary_df['protein rmsd mean'] = [r_df['protein rmsd'].mean()]
        summary_df['protein rmsd std'] = [r_df['protein rmsd'].std()]
        if not args.inhibitor in i_dict:
            summary_df['inhibitor rmsd mean'] = [r_df['inhibitor rmsd'].mean()]
            summary_df['inhibitor rmsd std'] = [r_df['inhibitor rmsd'].std()]
        print(summary_df)
        summary_df.to_csv(fname, index=False)
    return

df = run_rmsd_rep(rep_d,force=False)

data_fname = ("rmsd/data/" + sys_str + "rmsd.csv")
df.to_csv(data_fname, index=False)


summarise_rmsd(df)



