## First attempt for spatial clustering
## Weihua Guo, Ph.D.
## 05/27/2020

def dbscan_cust_func(xydf, slide_name, plot_pf, wsl_vis=False, sample_weight="NONE", min_pts=5):
    st = dt.now()
    res_dict = {}
    res_dict['slide_name'] = slide_name
    res_dict['sample_weight'] = sample_weight
    res_dict['plot_prefix'] = plot_pf
    print("Customized DBSCAN with optimal EPS selection algorithm...")
    print("\tSave directory: "+plot_pf)
    print("\tSetting sample weights: "+sample_weight)
    xydf['sample_weight'] = 1.0
    pheno_counts = xydf['phenotype'].value_counts()

    print("\tScaling...")
    wsl_xy = StandardScaler().fit_transform(xydf[['x','y']])
    print("\tPlotting k-distance graph for optimal eps...")
    neigh = NearestNeighbors(n_neighbors=min_pts)
    nbrs = neigh.fit(wsl_xy)
    distances, indices = nbrs.kneighbors(wsl_xy)
    dist_df = pd.DataFrame(distances, columns=[str(x) for x in range(0,min_pts)], index=xydf['phenotype'])
    dist_df = dist_df.sort_values(by=[str(min_pts-1)])
    dist_df = dist_df.reset_index()
    dist_df = dist_df.reset_index()

    kn = KneeLocator(x=dist_df['index'].tolist(), y=dist_df[str(min_pts-1)].tolist(), S=5.0,
            curve='convex', direction='increasing')
    print("\tOptimal eps is "+str(dist_df.loc[kn.knee, str(min_pts-1)]))
    res_dict['optimal_eps'] = dist_df.loc[kn.knee, str(min_pts-1)]

    melt_dist_df = pd.melt(dist_df, id_vars=['index', 'phenotype'], 
            value_vars=[str(x) for x in range(1,min_pts)])
    melt_dist_df['variable'] = melt_dist_df['variable'].astype("str")
    sns.scatterplot(data=melt_dist_df, x='index', y='value', hue='variable', 
            style='phenotype', linewidth=0, s=0.2)
    plt.vlines(kn.knee, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
    plt.hlines(dist_df.loc[kn.knee, str(min_pts-1)], plt.xlim()[0], plt.xlim()[1], linestyles='-.')
    plt.savefig(plot_pf+sample_weight+'_kdist_graph.png', dpi=300)
    plt.clf()
    plt.close()

    print("\tClustering...")
    db = DBSCAN(eps=dist_df.loc[kn.knee, str(min_pts-1)], 
            min_samples=min_pts, n_jobs=12).fit(wsl_xy, sample_weight=xydf['sample_weight'])
    labels = db.labels_
    res_df = xydf
    res_df = res_df.assign(dbscan_xy = labels)
    res_df['dbscan_clst'] = 'C' + res_df['dbscan_xy'].astype(str)
    res_df.to_csv(plot_pf+sample_weight+'_results.csv')

    print("\tAnalyzing cluster results...")
    res_dict['outlier_cell_num'] = res_df['dbscan_clst'].str.contains("C-").sum()

    res_dict['total_cell_num'] = xydf.shape[0]
    res_dict['outlier_rel_cell_num'] = res_dict['outlier_cell_num']/res_dict['total_cell_num']
    clst_cell_counts = res_df.groupby(['dbscan_clst']).size().reset_index(name='cell_num')
    res_dict["largest_clst_cell_num"] = clst_cell_counts['cell_num'].max()
    res_dict["largest_clst"] = clst_cell_counts.loc[clst_cell_counts['cell_num'].argmax(), 'dbscan_clst']
    out_clst = clst_cell_counts['dbscan_clst'][clst_cell_counts['dbscan_clst'].str.contains('-')].tolist()
    res_dict['total_clst_num'] = res_df['dbscan_clst'].value_counts().shape[0]-len(out_clst)


    clst_counts = res_df.groupby(['phenotype', 'dbscan_clst']).size().reset_index(name='cell_num')
    pivot_clst_cts = clst_counts.pivot(index='dbscan_clst', columns='phenotype', values='cell_num')
    if res_dict['total_clst_num'] == 0:
        figwidth = 3
    else:
        figwidth = res_dict['total_clst_num']
    pivot_clst_cts.plot.bar(stacked=True, logy=True, figsize=(0.15*figwidth,3))
    plt.savefig(plot_pf+sample_weight+'_stack_bar_pheno_clst.png', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    clst_counts.to_csv(plot_pf+sample_weight+'_raw_cluster_counts.csv')
    clst_counts = clst_counts.loc[~clst_counts['dbscan_clst'].str.contains('-')]
    clst_counts['r2pheno'] = clst_counts['cell_num'].div(clst_counts.groupby('phenotype')['cell_num'].transform('sum'))
    clst_counts['r2clst'] = clst_counts['cell_num'].div(clst_counts.groupby('dbscan_clst')['cell_num'].transform('sum'))

    clst_counts['rel_tag'] = 'mixed'
    only_mask = (clst_counts['r2clst'] == 1.00)&(clst_counts['phenotype'] == 'Th1')
    clst_counts.loc[only_mask, 'rel_tag'] = 'th1_only'
    res_dict['th1_only_cluster'] = only_mask.sum()
    only_mask = (clst_counts['r2clst'] == 1.00)&(clst_counts['phenotype'] == 'Th2')
    clst_counts.loc[only_mask, 'rel_tag'] = 'th2_only'
    res_dict['th2_only_cluster'] = only_mask.sum()
    only_mask = (clst_counts['r2clst'] == 1.00)&(clst_counts['phenotype'] == 'DC')
    clst_counts.loc[only_mask, 'rel_tag'] = 'dc_only'
    res_dict['dc_only_cluster'] = only_mask.sum()

    for iclst in clst_counts['dbscan_clst'].unique().tolist():
        tmp_mask = clst_counts['dbscan_clst']==iclst
        tmp_cts = clst_counts.loc[tmp_mask]
        tmp_only_mask = tmp_cts['rel_tag'].str.contains('only')
        if tmp_cts.shape[0] == 2:
            if any(tmp_cts['phenotype'].str.contains("DC")):
                if any(tmp_cts['phenotype'].str.contains("Th1")):
                    clst_counts.loc[tmp_mask, 'rel_tag'] = 'dc_th1_only'
                else:
                    clst_counts.loc[tmp_mask, 'rel_tag'] = 'dc_th2_only'
            else:
                clst_counts.loc[tmp_mask, 'rel_tag'] = 'th1_th2_only'
        elif tmp_cts.shape[0] == 3:
            th1_val = tmp_cts.loc[tmp_cts['phenotype'] == "Th1", 'r2pheno'].tolist()[0]
            th2_val = tmp_cts.loc[tmp_cts['phenotype'] == "Th2", 'r2pheno'].tolist()[0]

            if th1_val > th2_val:
                clst_counts.loc[tmp_mask, 'rel_tag'] = 'th1_domed'
            elif th1_val > th2_val:
                clst_counts.loc[tmp_mask, 'rel_tag'] = 'equal_ths'
            else:
                clst_counts.loc[tmp_mask, 'rel_tag'] = 'th2_domed'
    clst_counts.to_csv(plot_pf+sample_weight+'_cluster_counts.csv')
    ndc_counts = clst_counts.drop_duplicates(subset='dbscan_clst')

    res_dict['dc_th1_only_clst'] = (ndc_counts['rel_tag'] == 'dc_th1_only').sum()
    res_dict['dc_th2_only_clst'] = (ndc_counts['rel_tag'] == 'dc_th2_only').sum()
    res_dict['th1_th2_only_clst'] = (ndc_counts['rel_tag'] == 'th1_th2_only').sum()
    res_dict['th1_domed_clst'] = (ndc_counts['rel_tag'] == 'th1_domed').sum()
    res_dict['th2_domed_clst'] = (ndc_counts['rel_tag'] == 'th2_domed').sum()
    res_dict['eq_ths_clst'] = (ndc_counts['rel_tag'] == 'equal_ths').sum()
    res_dict['th1_sum_clst'] = res_dict['dc_th1_only_clst'] + res_dict['th1_domed_clst']
    res_dict['th2_sum_clst'] = res_dict['dc_th2_only_clst'] + res_dict['th2_domed_clst']

    if wsl_vis:
        clst_num = res_df['dbscan_clst'].value_counts().shape[0]
        if clst_num <= 20:
            npcol = 2
        elif clst_num <=40:
            npcol = 4
        elif clst_num <= 60:
            npcol = 6
        elif clst_num <= 80:
            npcol = 7
        elif clst_num <= 100:
            npcol = 8
        else:
            npcol = 10
        print("\tPlotting for phenotype...")
        wsl_pheno_plot = sns.scatterplot(x="x", y="y", hue="phenotype", data=res_df,
                s=0.15, alpha=0.5, linewidth=0, hue_order=['Th1','DC','Th2'])
        figure = wsl_pheno_plot.get_figure()
        figure.savefig(plot_pf+'phenotype_wsl.png', dpi=1200)
        plt.clf()
        plt.close()

        print("\tPlotting for DBSCAN...")
        sns.set_palette("colorblind")
        sns.scatterplot(x="x", y="y", hue="dbscan_clst", data=res_df, s=0.15, alpha=0.5, linewidth=0)
        plt.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', ncol=npcol)
        plt.savefig(plot_pf+sample_weight+'_cluster_wsl.png', dpi=1200, bbox_inches='tight')
        plt.clf()
        plt.close()

        print("\tPlotting for DBSCAN plus phenotype...")
        sns.scatterplot(x="x", y="y", hue="dbscan_clst", data=res_df,
                style="phenotype", alpha=0.5, s=0.15, linewidth=0)
        plt.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', ncol=npcol)
        plt.savefig(plot_pf+sample_weight+'_phenotype_dbscan_wsl.png', dpi=1200, bbox_inches='tight')
        plt.clf()
        plt.close()

        print("\tPlotting for DBSCAN plus phenotype (Facet)...")
        g = sns.FacetGrid(col="dbscan_clst", hue="phenotype", data=res_df, 
                col_wrap=9, 
                sharex=True, sharey=True)
        g.map(sns.scatterplot, 'x','y', linewidth=0, alpha=0.2, s=2.0)
        g.add_legend()
        plt.savefig(plot_pf+sample_weight+'_phenotype_dbscan_wsl_split.png', 
                dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()
    else:
        print("\tNo plot will be generated...")
    print("\t"+str(dt.now()))
    print("Time cost: "+str(dt.now()-st))
    print("+-"*45)
#    print(res_dict)
    return res_dict

import os
os.system('clear')
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import product
from config import config
from datetime import datetime as dt
from scipy.spatial import distance as distc
from sklearn.cluster import DBSCAN
from kneed import KneeLocator
from sklearn.neighbors import NearestNeighbors
import sklearn.utils
from sklearn.preprocessing import StandardScaler

data_dir = "/home/weihua/mnts/group_plee/Weihua/NanoString_Hot_Cold/dc_lamp_th12/qupath_processed/distant_analysis_simple_v4"
res_dir = "/home/weihua/mnts/group_plee/Weihua/NanoString_Hot_Cold/dc_lamp_th12/qupath_processed/cluster_test_v4/"
min_pts = 5

allTxt = glob.glob(data_dir+'/*.txt')
merge_res_dict = {}
for itxt in allTxt:
    txt_name = itxt.split('/')[-1].split('.t')[0]
    res_folder = "_".join(txt_name.split('_')[0:-5]+['dbscan'])
    exp_dir = os.path.join(res_dir, res_folder)
    if not os.path.exists(exp_dir):
        os.mkdir(exp_dir)
    tmp_pf = exp_dir+'/'+res_folder+'_'
    sn = '_'.join(tmp_pf.split('/')[-1].split('_')[4:-2])
#    print(tmp_pf)
    tmp_xy = pd.read_csv(itxt, sep = "\t")
    tmp_res_dict = dbscan_cust_func(xydf=tmp_xy, slide_name = sn, plot_pf=tmp_pf)
    merge_res_dict[tmp_res_dict['slide_name']+'_all_none'] = tmp_res_dict

    merge_res_df = pd.DataFrame.from_dict(merge_res_dict)
    print(merge_res_df)
    print(":) "*20+"\n\n\n")
merge_res_df = pd.DataFrame.from_dict(merge_res_dict)
merge_res_df.to_csv(res_dir+"final_merged_results_06012020.csv")
