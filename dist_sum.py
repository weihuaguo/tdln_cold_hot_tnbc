## Summarize the cell-wise distances into sample-wise distances
## Output files for statistical test

## Weihua Guo, Ph.D
## 05/09/2021

def mean_return(x):
    result = {
            'cell_number': x['count'].sum(),
            'dist_mean': x['sum'].sum()/x['count'].sum(),
            'dist_std': x['mean'].std(),
            }
    return pd.Series(result)
def final_mean_return(x):
    result = {
            'cell_number': x['cell_number'].sum(),
            'dist_mean': (x['cell_number']*x['dist_mean']).sum()/x['cell_number'].sum(),
            'dist_std': x['dist_std'].mean(),
            }
    return pd.Series(result)

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
print(case_info)

all_spt_files=glob.glob(data_dir+'/'+fid)
all_dist_files=glob.glob(input_dir+"/*_dist_"+out_name+".pkl")
all_dist_files=glob.glob(input_dir+"/*_neighbor_summary_"+out_name+".csv")

# print(all_dist_files)
ic=0
for idf in all_dist_files:
    ipid=idf.split('/')[-1].split('_neigh')[0]
    print("Processing "+ipid)
    print("\tReading distances...")
    tmp_df=pd.read_csv(idf, index_col=0)
    print(tmp_df.shape)
    tmp_df['sum']=tmp_df['count']*tmp_df['mean']
    tmp_res=tmp_df.groupby(['ngb','radius'],as_index=False).apply(mean_return)
#    print(tmp_res)

    tmpt_res=tmp_df.groupby(['radius'],as_index=False).apply(mean_return)
    tmpt_res['ngb']='Total CD8'
#    print(tmpt_res)

    tmp_res=pd.concat([tmp_res,tmpt_res])
    tmp_res['pid']=ipid
#    print(tmp_res)
    tmp_res['cohort']=case_info.loc[case_info['pid']==ipid,'cohort'].values.tolist()[0]
    tmp_res['Relapsed']=case_info.loc[case_info['pid']==ipid,'Relapsed'].values.tolist()[0]
    tmp_res['RFST']=case_info.loc[case_info['pid']==ipid,'RFST'].values.tolist()[0]
    tmp_res['Deceased']=case_info.loc[case_info['pid']==ipid,'Deceased'].values.tolist()[0]
    tmp_res['OST']=case_info.loc[case_info['pid']==ipid,'OST'].values.tolist()[0]

    if ic==0:
        merge_res=tmp_res
    else:
        merge_res=pd.concat([merge_res, tmp_res])
    ic+=1
merge_res.to_csv(output_dir+'/all_avg_distance_perSample_perDC_results.csv')

final_res=merge_res.groupby(['ngb','radius','cohort'],as_index=False).apply(final_mean_return)
print(final_res)
final_res.to_csv(output_dir+'/all_avg_distance_perCohort_perSample_perDC_results.csv')
