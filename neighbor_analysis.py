## Weihua Guo, Ph.D.
## 12/14/2023, Sat., 09:25 PM

import os
os.system('clear')

import glob
import numpy as np
import pandas as pd
import datetime as dt
import itertools as its
import scipy.spatial.distance as ssd

download_dir = "..."# NOTE: Please change the ... into the directory you were saving the data
download_dir = "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/ihc"# NOTE: Please change the ... into the directory you were saving the data

panel_dir = download_dir+"/mDC_Th_Panel"
in_dir = panel_dir+"/annotated_detections"
fr_dir = panel_dir+"/fixed_radius"
nn_dir = panel_dir+"/nearest_neighbors"

fid='*_extracted_detection.csv'
all_detect_files=glob.glob(in_dir+'/'+fid)
frs=[10,20,30,40,50]

fr_flag=False
nn_flag=True

print("Start to extract the spatial information...\n")
for idf in all_detect_files:
    ipid=idf.split('/')[-1].split('_ex')[0]
    print("Processing "+ipid)
    print("\tReading the QuPath output...")
    tst=dt.datetime.now()
    tmp_df=pd.read_csv(idf, index_col=0)
    print(tmp_df)
    print(tmp_df.columns.tolist())
    if nn_flag:
        nn_df = tmp_df[['Image', 'Name', 'Class', 'Parent', 'ROI', 'Centroid X µm', 'Centroid Y µm', 'phenotype']]
        for iph in nn_df['phenotype'].unique():
            nn_df[iph]=0
        for ir in range(tmp_df.shape[0]):
            tmp_x=tmp_df['Centroid X µm'].iloc[ir]
            tmp_y=tmp_df['Centroid Y µm'].iloc[ir]
            for iph in nn_df['phenotype'].unique():
                ph_df=tmp_df.loc[tmp_df['phenotype']==iph,:]
                ph_df['dist']=ssd.cdist(XA=[[tmp_x, tmp_y]], XB=ph_df.loc[:,ph_df.columns.str.contains('Centroid')])[0]
                ph_df = ph_df.loc[ph_df['dist'] > 0.0,:]
                min_idx = ph_df[['dist']].idxmin().tolist()
                nn_df[iph].iloc[ir]=ph_df[['dist']].loc[min_idx[0]]
        nn_df.to_csv(nn_dir+"/"+ipid+"_nearest_neighbor_distances.csv")

    if fr_flag:
        for ifr in frs:
            fr_df = tmp_df[['Image', 'Name', 'Class', 'Parent', 'ROI', 'Centroid X µm', 'Centroid Y µm', 'phenotype']]
            fr_df['fixed_r'] = ifr
            for iph in fr_df['phenotype'].unique():
                fr_df[iph]=0
            print(fr_df)
            for ir in range(tmp_df.shape[0]):
                print(ir)
                tmp_x=tmp_df['Centroid X µm'].iloc[ir]
                tmp_y=tmp_df['Centroid Y µm'].iloc[ir]

                x_mask=(tmp_df['Centroid X µm']>=tmp_x-ifr) & (tmp_df['Centroid X µm']<=tmp_x+ifr) & (tmp_df['Centroid X µm']!=tmp_x)
                y_mask=(tmp_df['Centroid Y µm']>=tmp_y-ifr) & (tmp_df['Centroid Y µm']<=tmp_y+ifr) & (tmp_df['Centroid Y µm']!=tmp_y)
                
                tmp_nghb_df = tmp_df.loc[x_mask&y_mask]
                if tmp_nghb_df.shape[0]>0:
                    for inghb in range(tmp_nghb_df.shape[0]):
                        fr_df[tmp_nghb_df['phenotype'].iloc[inghb]].iloc[ir] += 1
            fr_df.to_csv(fr_dir+"/"+ipid+"_fr"+str(ifr)+"_cell_count.csv")
