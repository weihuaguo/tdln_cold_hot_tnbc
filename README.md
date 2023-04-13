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
9. dist_calculation.py (RT ~ 3 days) Note: Must run after cell_pheno_extract.py
Calculat the average distance between mDC and other neighbor cells for each mDC in each slide
10. dist_sum.py (RT ~ 30 mins) Note: Must run after dist_calculation.R
Summarize the average distances and output the data for student's t-test
