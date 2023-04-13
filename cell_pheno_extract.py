## Read the QuPath detections (cell-wise annotation)
## Convert the QuPath cell phenotyping names into cell names of interests
## Extract the spatial information of the cells of interests

## Weihua Guo, Ph.D.
## 5/8/2021, Sat., 2:27 PM

import os
os.system('clear')

import glob
import numpy as np
import pandas as pd
import datetime as dt
import itertools as its
import scipy.spatial.distance as ssd

download_dir = "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir = download_dir+"/input_file"
out_dir = download_dir+"/Fig3/annotated_detection"

detect_dir=in_dir+"/qupath_detection" # Detection result directory
fid='*_qupath_cell_classifier_results.txt'

all_detect_files=glob.glob(detect_dir+'/'+fid)
cell_pheno_chart=pd.read_excel(in_dir+"/Input3_panel_cell_type_annot_v2.xlsx", index_col=0, engine='openpyxl')

print("Prepare for cell phenotype merging...\n")
cell_pheno_dict={}
all_classes_oi=[]
for ict, ir in cell_pheno_chart.iterrows():
    print(ict)
    y_markers=cell_pheno_chart.columns[ir=='Y'].tolist()
    na_markers=cell_pheno_chart.columns[ir.isna()].tolist()
    all_na_comb=[]
    for i in range(len(na_markers)+1):
        combl=list(its.combinations(na_markers,i))
        all_na_comb+=combl
    all_classes=[]
    for inc in all_na_comb:
        tmp_markers=list(its.permutations(y_markers+list(inc)))
        tmp_class=[",".join(list(i))+" positive" for i in tmp_markers] # NOTE: Check the QuPath output format
#        print(tmp_class)
        all_classes+=tmp_class
    cell_pheno_dict[ict]=all_classes
    if len(list(set(all_classes)&set(all_classes_oi)))!=0:
        warning("Overlapped cell phenotyping!!!")
    all_classes_oi+=all_classes
# print(cell_pheno_dict)

print("Start to extract the spatial information...\n")
for idf in all_detect_files:
    ipid=idf.split('/')[-1].split('_q')[0]
    print("Processing "+ipid)
    print("\tReading the QuPath output...")
    tst=dt.datetime.now()
    tmp_df=pd.read_table(idf, sep = '\t')
#    print(tmp_df['Class'].value_counts())
    tmp_df['pid']=ipid
    tmp_df['phenotype']='Others'
    tmp_df.loc[tmp_df['Class'].isna(),'phenotype']='Unknown'
    print("\t::Reading cost "+str(dt.datetime.now()-tst))
    # Cell phenotype merging
    print("\tMerging the phenotypes...")
    tst=dt.datetime.now()
    for key in cell_pheno_dict:
        tmp_mask=tmp_df['Class'].isin(cell_pheno_dict[key])
        tmp_df.loc[tmp_mask, 'phenotype']=key
    print("\t::Merging cell phenotypes cost "+str(dt.datetime.now()-tst))

#    print(tmp_df['phenotype'].value_counts())
    cln_df=tmp_df.loc[~tmp_df['Class'].isna(),:]
    cln_df=cln_df.loc[cln_df['phenotype']!='Others',:]
    cln_df.to_csv(out_dir+"/"+ipid+"_extracted_detection.csv")
    cln_spt_df=cln_df[['Image', 'Class', 'Centroid X µm', 'Centroid Y µm', 'pid', 'phenotype']]
    cln_spt_df.to_csv(out_dir+"/"+ipid+"_extracted_spatial_info.csv") 
