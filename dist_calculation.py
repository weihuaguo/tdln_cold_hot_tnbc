## Calculate the distances between the center cells and other cells of interests
## Calculate the average distances per sample

## Weihua Guo, Ph.D.
## 05/08/2021 3:47 PM

import os
os.system('clear')
import warnings
warnings.filterwarnings('ignore')

import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime as dt
from scipy.spatial import distance as distc
from sklearn.cluster import DBSCAN
from kneed import KneeLocator
from sklearn.neighbors import NearestNeighbors
import sklearn.utils
from sklearn.preprocessing import StandardScaler

data_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/extract_spatial_info"
res_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/dist_comp"
dist_type = 'euclidean'

fid='*_extracted_spatial_info.csv'
targets=['mDC']
neighbors=['CD8']
out_name='center_'+targets[0]+'_ngb_'+neighbors[0]
dist_cate = [25, 50, 100]
split_step=5000

all_spt_files=glob.glob(data_dir+'/'+fid)
# print(all_spt_files)

for ispt in all_spt_files:
    ipid=ispt.split('/')[-1].split('_e')[0]
    print("Processing "+ipid)
    print("\tReading spatial information...")
    tst=dt.now()
    tmp_df=pd.read_csv(ispt, index_col=0)
    target_df=tmp_df.loc[tmp_df['phenotype']==targets[0],:] # TODO: Multiple targets
    neighb_df=tmp_df.loc[tmp_df['phenotype'].str.contains(neighbors[0]),:]
    print("\t\t"+str(dt.now()-tst))

    '''
    print("\tCalculating distances...")
    tst=dt.now()
    tmp_dist=distc.cdist(target_df[['Centroid X µm', 'Centroid Y µm']], neighb_df[['Centroid X µm', 'Centroid Y µm']], metric=dist_type)
    tmp_dist = pd.DataFrame(tmp_dist, index=["mDC_"+str(x) for x in target_df.index.values.tolist()], 
            columns = [neighb_df.iloc[i,-1]+'_'+str(neighb_df.index.values.tolist()[i]) for i in range(0,neighb_df.shape[0])])
    tmp_dist.reset_index(inplace=True)
    print("\t\t"+str(dt.now()-tst))

    print("\tSaving distance dataframe...")
    tst=dt.now()
    tmp_dist.to_pickle(res_dir+"/"+ipid+"_dist_"+out_name+".pkl")
    print("\t\t"+str(dt.now()-tst))
    '''

    ## NOTE: Split due to the limited RAM!!!
    ## Every 5000 targets one time
    end_point=np.ceil(target_df.shape[0]/split_step)
#    print(end_point)
#    print(target_df.shape)
    sp=0
    ssst=dt.now()
    for ied in range(1, int(end_point+1)):
        ep=split_step*ied
        if ep>target_df.shape[0]:
            ep=target_df.shape[0]
        print("\t\tSubset: "+str(ep))

        sub_target=target_df.iloc[sp:ep]
        print("\tCalculating distances...")
        tst=dt.now()
        tmp_dist=distc.cdist(sub_target[['Centroid X µm', 'Centroid Y µm']], neighb_df[['Centroid X µm', 'Centroid Y µm']], metric=dist_type)
        tmp_dist = pd.DataFrame(tmp_dist, index=["mDC_"+str(x) for x in sub_target.index.values.tolist()], 
                columns = [neighb_df.iloc[i,-1]+'_'+str(neighb_df.index.values.tolist()[i]) for i in range(0,neighb_df.shape[0])])
        tmp_dist.reset_index(inplace=True)
        print("\t\t"+str(dt.now()-tst))


#        sub_dist=tmp_dist.iloc[sp:ep]
#        print(sub_dist)
        print("\tMelting the distance dataframe...")
        tst=dt.now()
        melt_dist = pd.melt(tmp_dist, 
                value_vars=[neighb_df.iloc[i,-1]+'_'+str(neighb_df.index.values.tolist()[i]) for i in range(0,neighb_df.shape[0])], 
                id_vars = ['index'])
        print("\t\t"+str(dt.now()-tst))

        melt_dist=melt_dist.loc[melt_dist['value']<=max(dist_cate),:]
        icts=0
        dddst=dt.now()
        for ids in dist_cate:
            print("\t\t\tRadius: "+str(ids))
            tmp_melt=melt_dist.loc[melt_dist['value']<=ids,:]
            ngb_split=tmp_melt['variable'].str.split("_", n=2, expand=True)
            tmp_melt[['ngb', 'ingb']]=tmp_melt['variable'].str.split("_", expand=True)
            tmp_melt['radius']='Within '+str(ids)+' microns'
            tmp_res=pd.DataFrame(tmp_melt.groupby(['index','ngb'])['value'].describe().reset_index())
            tmp_res['radius']='Within '+str(ids)+' microns'
            if icts==0:
                merge_res=tmp_res
                merge_melt=tmp_melt
            else:
                merge_res=pd.concat([merge_res, tmp_res])
                merge_melt=pd.concat([merge_melt, tmp_melt])
            icts+=1
            print("\t\t\t\t"+str(dt.now()-dddst))
        merge_res['pid']=ipid
        if sp==0:
            final_res=merge_res
            final_dist=tmp_dist
            final_melt=merge_melt
        else:
            final_res=pd.concat([final_res,merge_res])
            final_dist=pd.concat([final_dist,tmp_dist])
            final_melt=pd.concat([final_melt, merge_melt])
        print(tmp_dist.shape)
        sp=ep
    print("\tSaving distance dataframe...")
    tst=dt.now()
    final_dist.to_pickle(res_dir+"/"+ipid+"_dist_"+out_name+".pkl")
    print("\t\t"+str(dt.now()-tst))

    print("\t\tTotal time cost for summarization: "+str(dt.now()-ssst))
    final_res.to_csv(res_dir+"/"+ipid+"_neighbor_summary_"+out_name+".csv")
    print(final_res.shape)
    print(target_df.shape)
    print(neighb_df.shape)
'''
    tst=dt.now()
    melt_dist = pd.melt(tmp_dist, 
            value_vars=[neighb_df.iloc[i,-1]+'_'+str(neighb_df.index.values.tolist()[i]) for i in range(0,neighb_df.shape[0])], 
            id_vars = ['index'])
    print("\t\t"+str(dt.now()-tst))

    print("\tSaving melted dataframe...")
    tst=dt.now()
#    melt_dist.to_pickle(res_dir+"/"+ipid+"_melt_dist_"+out_name+".pkl")
    print("\t\t"+str(dt.now()-tst))

    melt_dist=melt_dist.loc[melt_dist['value']<=max(dist_cate),:]
    icts=0
    for ids in dist_cate:
        print("\t\t\tRadius: "+str(ids))
        tmp_melt=melt_dist.loc[melt_dist['value']<=ids,:]
        ngb_split=tmp_melt['variable'].str.split("_", n=2, expand=True)
        tmp_melt[['ngb', 'ingb']]=tmp_melt['variable'].str.split("_", expand=True)
        tmp_res=pd.DataFrame(tmp_melt.groupby(['index','ngb'])['value'].describe().reset_index())
        tmp_res['radius']='Within '+str(ids)+' microns'
        if icts==0:
            merge_res=tmp_res
        else:
            merge_res=pd.concat([merge_res, tmp_res])
        icts+=1
    merge_res['pid']=ipid
    merge_res.to_csv(res_dir+"/"+ipid+"_neighbor_summary_"+out_name+".csv")
    print(tmp_dist.shape)
    print(target_df.shape)
    print(neighb_df.shape)
'''
