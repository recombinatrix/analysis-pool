#!/usr/bin/env python                                                       
# coding: utf-8 

import numpy as np                                                              
import pandas as pd                                                             

# get file list

model_d = {
    'hOCT1' : 'D357-',
    'hOCT1-ASP357ASPH' : 'D357H'
}

import os
flist = os.listdir()
cflist = [f for f in flist if "xvg" in f and "#" not in f]

# load data file by file

rep = np.repeat([1,2,3],501)
r_df = pd.DataFrame(columns=["Timestep","RMSD"])
for i in cflist:
    # parse filename

    prot=model_d[i.split("_")[0]]
    lig=i.split("_")[1]
    pose=i.split("_")[2]
    fit=i.split("_")[4]
    calc=i.split("_")[6]


    df = pd.read_csv(i,names=["Timestep","RMSD"],dtype=float,delim_whitespace=True)
    df["model"]=prot
    df["ligand"]=lig
    df["pose"]=pose
    df["calc for"]=calc
    df["fit to"]=fit
    df["Timestep"]=df.loc[:,"Timestep"] / 1000 
    df["replicate"]=rep.astype(int)

    r_df = pd.concat((r_df,df))

# data is now ready for processing
r_df.to_csv("../RMSD_data.csv", index=False)

