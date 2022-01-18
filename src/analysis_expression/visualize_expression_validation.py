# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses scatter plots to visualize the gene expression-validation of DTox interpretation results on Tox21 datasets, comparing the compound differential expression proportion compounds among significant DTox paths vs among background DTox paths 


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests


## 0. Input arguments
result_file	= 'data/compound_target_probability_tox21_interpret_expression/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_expression_validation_all_result.tsv'	# name of input compound differential expression proportion result file 
pv_file		= 'data/compound_target_probability_tox21_interpret_expression/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_expression_validation_result_compare_summary.tsv'	# name of input proportion comparison p-value result file  
plot_file	= 'plot/compound_target_probability_tox21_interpret_expression/compound_target_fingerprint_maccs_probability_tox21_interpret_expression_validation_result_compare_'	# folder name of output plot file 
query_fdr	= 0.05	# number that indicates query pathway FDR threshold for extracting validation results of interest  
query_rule	= 'gamma-epsilon_0.001_0.1'	# string that indicates query LPR propagation rule along with parameter values for extracting validation results of interest
query_dt	= [('10uM', '6h'), ('10uM', '24h'), ('1.11uM', '24h')]	# list that indicates query dose and time combinations for extracting validation results of interest
query_assay	= ['tox21-aromatase-p1', 'tox21-mitotox-p1', 'tox21-pxr-p1', 'tox21-rt-viability-hepg2-p1']	# list that indicates query Tox21 asssay datasets for extracting validation results of interest 

## 1. Process input validation result files   
# read in differential expression proportion of compounds as data frame, extract rows of interest based on specified LRP propagation rule and FDR threshold  
result_df = pd.read_csv(result_file, sep = '\t', header = 0)
result_df = result_df[result_df.rule == query_rule]
result_df = result_df[result_df.fdr_threshold == query_fdr]
# read in proportion comparison among significant DTox paths vs among background DTox paths as data frame, extract rows of interest based on specified LRP propagation rule and FDR threshold    
pv_df = pd.read_csv(pv_file, sep = '\t', header = 0) 
pv_df = pv_df[pv_df.rule == query_rule]
pv_df = pv_df[pv_df.fdr_threshold == query_fdr]
# iterate by dose and time combination, extract rows of interest from proportion comparison data frame based on each dose and time values  
qd_pv_df_list = []
for qd in query_dt:
	qd_pv_df = pv_df[(pv_df.dose == qd[0]) & (pv_df.time == qd[1])]
	qd_pv_df_list.append(qd_pv_df)	
pv_df = pd.concat(qd_pv_df_list)

## 2. Specify size and parameters of scatter plot 
fig, gs_ax = plt.subplots(len(query_dt), len(query_assay),  figsize = (6*len(query_assay), 6*len(query_dt)))
plt.rc('font', size = 35)
xy_lim = [[0, 0.3], [0, 0.1], [0, 0.3], [0, 0.25]]
xy_line = [[0, 0.3], [0, 0.1], [0, 0.3], [0, 0.25]]
xy_tick = [[0, 0.1, 0.2, 0.3], [0, 0.05, 0.1], [0, 0.1, 0.2, 0.3], [0, 0.1, 0.2]]

## 3. Make scatter plots 
# iterate by query Tox21 assay dataset, make scatter plots showing the comparison of differential expression proportion of compounds from each assay      
for lqa in range(0, len(query_assay)):
	# select rows relevant to current Tox21 assay from compound differential expression proportion data frame   
	lqa_assay = query_assay[lqa]
	lqa_result_df = result_df[result_df.assay == lqa_assay]
	# obtain the t test P-value for the comparison for current Tox21 assay, use FDR to perform multiple testing correction
	lqa_pv = pv_df[pv_df.assay == lqa_assay].t_pv.values	
	lqa_reject, lqa_fdr, _, _ = multipletests(lqa_pv, method = 'fdr_bh')
	# iterate by query dose-time combination, make scatter plot showing the comparison of differential expression proportion of compounds from each combination 
	for lqd in range(0, len(query_dt)):
		# obtain the differential expression proportion among significant DTox paths and among background DTox paths for the current dose-time combination  
		lqd_dt = query_dt[lqd]
		lqd_result_df1 = lqa_result_df[lqa_result_df.dose == lqd_dt[0]]
		lqd_result_df = lqd_result_df1[lqd_result_df1.time == lqd_dt[1]]
		lqd_x = lqd_result_df.background_sig_ratio.values
		lqd_y = lqd_result_df.interpret_sig_ratio.values
		# specify color for each proportion point: blue for significant > background, gray otherwise  
		lqd_point_color = []
		for lx in range(0, len(lqd_x)):
			if lqd_y[lx] > lqd_x[lx]:
				lqd_point_color.append('blue')
			else:
				lqd_point_color.append('gray')
		# make scatter plot showing the comparison between differential expression proportion of compounds among significant DTox paths vs among background DTox paths 
		gs_ax[lqd][lqa].scatter(lqd_x, lqd_y, s = 15, c = lqd_point_color)
		# add dashed diagonal line  
		lqd_y_max = np.ceil(np.max(lqd_y)/0.05) * 0.05
		lqd_xy_line = [0, lqd_y_max]
		gs_ax[lqd][lqa].plot(lqd_xy_line, lqd_xy_line, '--k')
		# specify label and range for axes  
		lqd_xy_lim = [-lqd_y_max/20, lqd_y_max*1.05]
		gs_ax[lqd][lqa].set_xlim(lqd_xy_lim)
		gs_ax[lqd][lqa].set_xlabel('Expected DE paths %', fontsize = 25)
		gs_ax[lqd][lqa].set_ylim(lqd_xy_lim)
		gs_ax[lqd][lqa].set_ylabel('Observed DE paths %', fontsize = 25)
		# add text strings that specifies FDR value (red if FDR < 0.05, black otherwise)  
		lqd_fdr_char = 'FDR = ' + str(np.round(lqa_fdr[lqd], 2))
		if lqa_reject[lqd]:
			lqa_col = 'red'
		else:
			lqa_col = 'black'
		gs_ax[lqd][lqa].text(0.8, 0.05, lqd_fdr_char, horizontalalignment = 'center', verticalalignment = 'center', transform = gs_ax[lqd][lqa].transAxes, size = 22, c = lqa_col)
		# specify ticks for axes 
		if lqd_y_max < 0.16:
			lqd_ticks = np.arange(0, 1 + lqd_y_max/0.05) * 0.05
		else:
			lqd_ticks = np.arange(0, 1 + np.floor(lqd_y_max/0.1)) * 0.1
			if lqd_y_max/0.05 % 2 == 1:
				lqd_ticks = np.insert(lqd_ticks, len(lqd_ticks), lqd_y_max)	
		gs_ax[lqd][lqa].set_xticks(lqd_ticks)
		gs_ax[lqd][lqa].set_yticks(lqd_ticks)
		gs_ax[lqd][lqa].tick_params(axis = 'both', which = 'major', labelsize = 22)
		# remove axes on the top and right   
		gs_ax[lqd][lqa].spines['right'].set_visible(False)
		gs_ax[lqd][lqa].spines['top'].set_visible(False)

## 4. Save boxplot
plt.tight_layout()
plt.savefig(plot_file + query_rule + '.pdf')
plt.close()

