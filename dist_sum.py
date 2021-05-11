## Summarize the cell-wise distances into sample-wise distances
## Output files for statistical test

## Weihua Guo, Ph.D
## 05/09/2021

def mean_return(x):
    result = {
            'cell_number': x['count'].sum(),
            'dist_mean': x['sum'].sum()/x['count'].sum()
            }
    return pd.Series(result).round(0)

import os
os.system('clear')
import warnings
warnings.filterwarnings('ignore')

import glob
import pandas as pd
import numpy as np
from datetime import datetime as dt

data_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/extract_spatial_info"
input_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/dist_comp"
output_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/dist_sum"

dist_type = 'euclidean'

fid='*_extracted_spatial_info.csv'
targets=['mDC']
neighbors=['CD8']
out_name='center_'+targets[0]+'_ngb_'+neighbors[0]
dist_cate = [25, 50, 100, 200]

case_info=pd.read_excel("/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/patient_annotation_v2.xlsx", sheet_name='Sheet1')
# print(case_info)

all_spt_files=glob.glob(data_dir+'/'+fid)
all_dist_files=glob.glob(input_dir+"/*_dist_"+out_name+".pkl")
# print(all_dist_files)

for idf in all_dist_files:
    ipid=idf.split('/')[-1].split('_dist')[0]
    print("Processing "+ipid)
    print("\tReading distances...")
    tmp_df=pd.read_pickle(idf)
    print(tmp_df.shape)
    exit()
    tmp_df['sum']=tmp_df['count']*tmp_df['mean']
    tmp_res=tmp_df.groupby(['ngb','radius'],as_index=False).agg({'count':'sum', 'sum': 'sum'})
    tmp_res['pid_mean']=tmp_res['sum']/tmp_res['count']
    print(tmp_df)
    print(tmp_res)
    tmp_res=tmp_df.groupby(['ngb','radius']).apply(mean_return)
    print(tmp_res)

    exit()

