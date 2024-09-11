import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import skimage
import sklearn
import spateo as st
import dynamo as dyn
import cv2
import anndata as ad
import pandas as pd
from scipy import stats
import os


import scanpy as sc
from collections import Counter
import anndata as ad
import loompy as lp
import numpy as np

import scanpy as sc


#####set pathway
os.chdir("/cellshape/")
rawdir = "/h5ad_RAw/"
csvdir = "/meta/"

#####let us creat processed cycle code to finish cellsegmentation
for key in [i [:-8] for i in os.listdir(csvdir)]:
    print(key)
    key1 = key[:-2]
    print(key1)
    df = pd.read_csv(csvdir+key+"meta.csv",index_col=0)
    data.obs.index=data.obs.index.astype("int64")
    data.obs=pd.merge(data.obs,df,right_index=True,left_index=True,how="left")
    data.obs=data.obs.fillna("others")
    data.obs.index=data.obs.index.astype("object")
    data1=data[data.obs['first_type'].isin(['EC','Fibroblast','Immune','Macrophage','Modulated_SMC','SMC'])]
    data1.obs.index=data1.obs.index.astype("int64")
    
    
    ###Create RGB array
    red = np.zeros_like(data1.uns["seg"])
    green = np.zeros_like(data1.uns["seg"])
    blue = np.zeros_like(data1.uns["seg"])
    ###Extract the corresponding number for the cell type
    EC_list=data1.obs_names[data1.obs['first_type']=="EC"]
    Fibroblast_list=data1.obs_names[data1.obs['first_type']=="Fibroblast"]
    Immune_list=data1.obs_names[data1.obs['first_type']=="Immune"]
    Macrophage_list=data1.obs_names[data1.obs['first_type']=="Macrophage"]
    Modulated_SMC_list=data1.obs_names[data1.obs['first_type']=="Modulated_SMC"]
    SMC_list=data1.obs_names[data1.obs['first_type']=="SMC"]
    
    
    ###Assign RGB colors to each cell
    ###EC 
    red[np.isin(data1.uns["seg"], EC_list)] = 35
    green[np.isin(data1.uns["seg"], EC_list)] = 139
    blue[np.isin(data1.uns["seg"], EC_list)] = 69
    
    ###Fib
    red[np.isin(data1.uns["seg"], Fibroblast_list)] = 139
    green[np.isin(data1.uns["seg"], Fibroblast_list)] = 74
    blue[np.isin(data1.uns["seg"], Fibroblast_list)] = 75
    
    ###modulated_SMC
    red[np.isin(data1.uns["seg"], Modulated_SMC_list)] = 161
    green[np.isin(data1.uns["seg"], Modulated_SMC_list)] = 207
    blue[np.isin(data1.uns["seg"], Modulated_SMC_list)] = 250
    
    ###immune
    red[np.isin(data1.uns["seg"], Immune_list)] = 206
    green[np.isin(data1.uns["seg"], Immune_list)] = 18
    blue[np.isin(data1.uns["seg"], Immune_list)] = 86
    
    ###Macrophage
    red[np.isin(data1.uns["seg"], Macrophage_list)] = 223
    green[np.isin(data1.uns["seg"], Macrophage_list)] = 101
    blue[np.isin(data1.uns["seg"], Macrophage_list)] = 176
    
    #####SMC
    red[np.isin(data1.uns["seg"], SMC_list)] = 31
    green[np.isin(data1.uns["seg"], SMC_list)] = 120
    blue[np.isin(data1.uns["seg"], SMC_list)] = 180
    
    RGB = np.dstack((blue,green,red))
    
    cv2.imwrite(key+'_celltype.tif',RGB)
