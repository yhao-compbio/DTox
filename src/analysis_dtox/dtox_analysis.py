# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions used in DTox model result anaysis


## Module
import sys
import numpy as np
import pandas as pd


## The function finds optmial hyperparameter setting of supervised learning model based on input performance metric values   
def find_optimal_hyperparameter_setting(perf_df, metric_name, output_folder):
	## 0. input arguments 
		# perf_df: data frame that contains performance metric values of all hyperparameter settings  
		# metric_name: name of metric to be used for optmial hyperparameter identification   
		# output_folder: folder name of output hyperparameter setting result file  

	## 1. 
	# when metric of interest involves loss function, optimal setting resides in minimal metric value   
	if 'loss' in metric_name:
		a_order = True 
	# for other metrics, optmial setting resides in maximal metric value
	else:	
		a_order = False 
	# obtain names of all datasets, iterate by dataset  
	datasets = perf_df['dataset_name'].unique()
	optimal_perf_list = []
	for ds in datasets:
		# select performance metric values of all hyperparameter settings on the dataset   
		ds_perf_df = perf_df[perf_df['dataset_name'] == ds]
		# rank all performance metric values, identify optimal hyperparameter setting (min or max value, depend on metric) 
		ds_perf_df = ds_perf_df.sort_values(metric_name, axis = 0, ascending = a_order)	
		optimal_perf_list.append(ds_perf_df.iloc[0, :])
	# aggregate identified hyperparameter settings of all datasets in data frame form  
	optimal_perf_df = pd.concat(optimal_perf_list, axis = 1).T	
	# write optimal hyperparameter setting results to output file  
	output_file = output_folder + '_optimal_performance_summary_by_' + metric_name + '.tsv'
	optimal_perf_df.to_csv(output_file, sep = '\t', index = False) 

	return optimal_perf_df	
	
