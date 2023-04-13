# Weihua Guo, Ph.D.,
# For LN cold vs Hot clustering metrics
# 06/03/2020
def clst_metric(df, df_name, clst_nums = range(2,13)):
    print("Screening for "+df_name+" dataframe...")
    score_dict = {}
    for i in clst_nums:
        tmp_index = df_name+'_C'+str(i)
        score_dict[tmp_index] = {}
        km = KMeans(n_clusters=i, random_state=0).fit(df)
        preds = km.predict(df)
        
        print("Score for number of cluster(s) {}: {}".format(i,km.score(df)))
        score_dict[tmp_index]['km_scores'] = -km.score(df)
        
        silhouette = silhouette_score(df, preds)
        score_dict[tmp_index]['silhouette_score'] = silhouette
        print("Silhouette score for number of cluster(s) {}: {}".format(i,silhouette))
        
        db = davies_bouldin_score(df, preds)
        score_dict[tmp_index]['davies_bouldin_score'] = db
        print("Davies Bouldin score for number of cluster(s) {}: {}".format(i,db))

        print("-+"*45)
    score_df = pd.DataFrame.from_dict(score_dict).T
    score_df = score_df.reset_index()
    split_info = score_df['index'].str.split('_', n=2, expand=True)
    split_info.columns = ['cohort', 'clst_num']
    res_df = pd.concat([score_df, split_info], sort=False, axis=1)
    return res_df

import os
os.system('clear')

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score,v_measure_score

download_dir = "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir = download_dir+"/Fig2"
out_dir = download_dir+"/Fig2"


work_dir = "/home/weihua/mnts/group_plee/Weihua/nanostring_results/outline_and_lists/weihua_results/"
plot_dir = "/home/weihua/mnts/group_plee/Weihua/nanostring_results/outline_and_lists/plots/"

low_corr_file = "_corr_coef_low_pearson_immune_only_xself_.csv"
high_corr_file = "_corr_coef_high_pearson_immune_only_xself_.csv"

low_df = pd.read_csv(in_dir+low_corr_file, sep = ",", index_col=0)
cold_score_df = clst_metric(low_df, df_name="cold")

high_df = pd.read_csv(in_dir+high_corr_file, sep = ",", index_col=0)
hot_score_df = clst_metric(high_df, df_name="hot")

clst_score_df = pd.concat([cold_score_df, hot_score_df], axis=0, sort=False)
clst_score_df.to_csv(out_dir+"cluster_score_results_06032020.csv")
