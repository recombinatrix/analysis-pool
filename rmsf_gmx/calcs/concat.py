#!/usr/bin/env python                                                       
# coding: utf-8 

import matplotlib.pyplot as plt                                                 
import numpy as np                                                              
import pandas as pd                                                             
import scipy as sp
import seaborn as sns

# get file list

import os
flist = os.listdir()
# this was for putting apo first; not doing that anymore because who cares
# aflist = [f for f in flist if "xvg" in f and "apo" in f and "#" not in f]
# rflist = aflist + [f for f in flist if "xvg" in f and "#" not in f and "apo" not in f] # put apo first so I can calculate difference from apo

rflist = [f for f in flist if "xvg" in f and "#" not in f] # put apo first so I can calculate difference from apo
mdict = {'hOCT1-ASP357ASPH' : 'D357H', 'hOCT1' : 'D357-'}


# load data file by file
r_df = pd.DataFrame(columns=["resid","RMSF"])

for i in rflist:
    # parse filename

    model=mdict[i.split("_")[0]]
    lig=i.split("_")[1]
    pose=i.split("_")[2]



    df = pd.read_csv(i,names=["resid","RMSF"],dtype=float,delim_whitespace=True)
    df['resid']=df['resid'].astype(int)
    df["model"]=model
    df["ligand"]=lig
    df["pose"]=pose
    r_df = pd.concat((r_df,df))
    
# data is now ready for processing
r_df['ligand'] = r_df['ligand'].str.swapcase() # swap case so lowercase a comes before capitals 
r_df = r_df.sort_values(by=['ligand','pose','model','resid'])
r_df['ligand'] = r_df['ligand'].str.swapcase() # swap case to originalback

r_df.to_csv("../RMSF_data.csv", index=False)

