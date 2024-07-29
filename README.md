# tdln_cold_hot_tnbc
Scripts for comparing tumor draining lymph nodes with cold and hot primary triple negative breast tumors
RT: Run time based on a i9-9900KF CPU + 64GB RAM computer with Ubuntu LTS 20.04 OS.
Anaconda based virtual environment was used to run all the following scripts. Detailed package versions are available in Methods section of the manuscript.

1. limma_ns_de.R (RT ~ 3 mins)
Gene differential expression analysis (SF1 + SD1)
2. corr_analysis.R (RT ~ 15 mins)
Calculate the correlation coefficients for each two genes within cold and hot TDLN cohorts. (SD3)
3. heatmap_plotter.R (RT ~ 2 mins) Note: Must run after corr_analysis.R
Plot the correlation matrices with different clustering orders (Fig 2A).
4. cluster_metric.py (RT ~ 2 mins) Note: Must run after corr_analysis.R
Calculate silhouttee score for clustering results
5. leave_one_out_corr_analysis.R (RT ~ 120 mins)
Calculate the PCCs with each one sample exclusion
6. leave_one_out_merge_analysis.R (RT ~ 10 mins) Note: Must run after leave_one_out_corr_analysis.R
Concatenate the above LOO results and implement student's t-test between cold and hot TDLNs
7. leave_one_out_vis.R (RT ~ 10 mins) Note: Must run after leave_one_out_merge_analysis.R
Plot volcano plot for LOO results (SF2)
8. cell_pheno_extract.py (RT ~ 60 mins)
Annotate the cell phenotypes and extract the usable spatial information
9. annotation_examination_mdc_panel.R
Generate the ridge plot and heatmaps of marker expression for each annotated cell type (mDC-Th panel)
10. annotation_examination_il4_panel.R
Generate the ridge plot and heatmaps of marker expression for each annotated cell type (MC-IL4 panel)
11. auc_cold_hot.R
Draw AUC curves for mast cell validation from H-DAB staining
12. manual_calculate_cell_type_score.R
Calulcate cell type score
13. t_activate_prolifer_check.R
Examine the T cell activation/profiling/differential marker expression
14. gse88715_de.R
Differential expression on GSE88715 dataset
15. gse88715_timer_comp.R
Compare TIMER 2.0 results between cold and hot cohort
16. spatial_cluster.py
DBSCAN for each slide with k-distance graph
