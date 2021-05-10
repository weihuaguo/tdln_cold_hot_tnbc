## Summarize the cell-wise distances into sample-wise distances
## Output files for statistical test

## Weihua Guo, Ph.D
## 05/09/2021

import os
os.system('clear')
import warnings
warnings.filterwarnings('ignore')

import glob
import pandas as pd
import numpy as np
from datetime import datetime as dt

input_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/dist_comp"
output_dir="/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/dist_sum"

dist_type = 'euclidean'

fid='*_extracted_spatial_info.csv'
targets=['mDC']
neighbors=['Th']
out_name='center_'+targets[0]+'_ngb_'+neighbors[0]
dist_cate = [25, 50, 100, 200]

case_info=pd.read_excel("/home/weihua/mnts/data_plee/Group/weihua/dc_lamp_ln/post_qupath/patient_annotation_v2.xlsx")
print(case_info)
