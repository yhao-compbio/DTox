# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions for visualizing DTox model results 


## Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## This function uses barplot to visualize performance comparison across different machine learning method implementations on Tox21 datasets 
def visualize_dataset_performane_comparison(file_df, ds_annot_df, plot_files):
	## 0. Input arguments 
		# file_df: data frame that contains performance file names of different method implementations  
		# ds_annot_df: data frame taht cotains name annotations of Tox21 datasets 
		# plot_files: folder name of output visualization results  
	
	## 1. Read in performance data frame of each method implementation on Tox21 datasets  
	# iterate by method implementation 
	perf_df_list = []
	for i in range(0, file_df.shape[0]):
		# read in performance metric data frame of current method implementation  
		i_perf_df = pd.read_csv(file_df.iloc[i,:].perf_file, sep = '\t', header = 0)
		perf_df_list.append(i_perf_df)	

	## 2. Visualize performance comparison across different method implementations by barplot
	# specify metrics to be compared: AUROC and balanced accuracy  
	compare_measures = ['testing_auc', 'testing_bac']
	compare_measure_names = ['Area under ROC curve', 'Balanced accuracy']
	compare_measure_lims = [(0.5, 0.92), (0.5, 0.8)]
	# iterate by Tox21 dataset  
	for j in range(0, i_perf_df.shape[0]): 
		# obtain full name of current Tox21 assay dataset  
		j_ds = i_perf_df.dataset_name.values[j]
		j_ds_row = ds_annot_df.loc[ds_annot_df.protocol_name == j_ds].iloc[0,:]	
		j_ds_annot = j_ds_row.assay_target.split(' (')[0]
		# add cell line name to dataset name if it is a viability assay (to differentiate two viability assays)
		if j_ds_annot == 'Cell viability':
			j_ds_annot = j_ds_annot + ' (' + j_ds_row.cell_line + ')'
		# specify figure and font size 
		f = plt.figure(figsize = (2*len(compare_measures), 4))
		plt.rc('font', size = 18)
		plt.rc('axes', titlesize = 18)
		plt.rc('axes', labelsize = 18)
		plt.rc('xtick', labelsize = 14)
		plt.rc('ytick', labelsize = 15)
		plt.rc('legend', fontsize = 14)
		# iterate by performance metric of interest 
		for lcm in range(0, len(compare_measures)):
			# extract performance metric values and 95% CIs of different method implementations
			cm = compare_measures[lcm]
			cm_value = []
			cm_error = [] 
			for pdl in perf_df_list:
				cm_value.append(pdl[cm].values[j])
				cm_error.append(pdl[cm + '_ci'].values[j])
			x_pos = np.arange(0, len(cm_value))
			# make barplot showing the extracted performance metric values and 95% CIs 
			ax = f.add_subplot(int('1' + str(len(compare_measures)) + str(lcm+1)))
			ax.bar(x_pos, cm_value, yerr = cm_error, align = 'center', ecolor = 'black', capsize = 2)
			# remove axes on the top, bottom, and right
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			# set y axis (performance metric) range and label
			ax.set_ylim(compare_measure_lims[lcm])
			ax.set_ylabel(compare_measure_names[lcm])
			# set x axis (different methods) ticks 
			ax.set_xticks(x_posperformance metric)
			ax.set_xticklabels(file_df.method.values, rotation = 90)
		# add dataset name as title to the barplot   
		f.suptitle(j_ds_annot, size = 20) 
		# save barplot 
		j_out = plot_files + '_method_performance_comparison_' + j_ds + '.pdf'
		plt.tight_layout()
		plt.savefig(j_out)
		plt.close()
	
	return 1
