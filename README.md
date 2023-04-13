# tdln_cold_hot_tnbc
Scripts for comparing tumor draining lymph nodes with cold and hot primary triple negative breast tumors

1. limma_ns_de.R
Gene differential expression analysis (SF1 + SD1)
2. corr_analysis.R
Calculate the correlation coefficients for each two genes within cold and hot TDLN cohorts. (SD3)
3. heatmap_plotter.R
Plot the correlation matrices with different clustering orders (Fig 2A).
4. cluster_metric.py
Calculate silhouttee score for clustering results
5. leave_one_out_corr_analysis.R
Calculate the PCCs with each one sample exclusion
6. leave_one_out_merge_analysis.R
Concatenate the above LOO results and implement student's t-test between cold and hot TDLNs
7. leave_one_out_vis.R
Plot volcano plot for LOO results (SF2)
