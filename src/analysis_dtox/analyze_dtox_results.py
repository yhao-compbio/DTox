# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script identifies optimal hyperparameter setting of machine learning method implementation, then compares and visualizes model performance across different method implementations  


## Modules
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'src/analysis_dtox/')
import dtox_analysis
import dtox_plot


## Main function
def main(argv):
	## 0. Input arguments: 
		# argv 1: string that indicates type of analysis task: 'simple' for analyzing simple machine learning model results; 'dtox' for analyzing DTox model results   
		# argv 2: for 'simple' task: input file that contains model performance metric values of all hyperparameter settings; for 'dtox' task: input file that contains model performance file names of different method implementations (will compare basic DTox implementation with other method implementations in the analysis) 
		# argv 3: name of metric to be used for optimal hyperparameter identification (e.g. 'training_root_loss' for DTox models, 'training_log_loss' for simple machine learning models) 
		# argv 4: folder name of output analysis results
		# argv 5: folder name of output visualization results 
	
	## 1. Perform analysis according to the specified task  
	# analyze simple machine learning model results 
	if argv[1] == 'simple':
		# read in performance metric values of all hyperparameter settings as data frame, find optimal hyperparameter setting of each dataset  
		all_perf_df = pd.read_csv(argv[2], sep = '\t', header = 0)
		op_perf_df = dtox_analysis.find_optimal_hyperparameter_setting(all_perf_df, argv[3], argv[4])
	# analyze DTox model results 
	if argv[1] == 'dtox': 
		# read in performance file names of method implementations, obtain performance file name for DTox basic implementation
		perf_file_df = pd.read_csv(argv[2], sep = '\t', header = 0)
		imple_file = perf_file_df.loc[perf_file_df.method == 'VNN'].iloc[0,:].perf_file
		# read in performance metric values of all hyperparameter settings in DTox basic implementation as data frame, find optimal hyperparameter setting of each dataset 
		imple_perf_df = pd.read_csv(imple_file, sep = '\t', header = 0)
		op_perf_df = dtox_analysis.find_optimal_hyperparameter_setting(imple_perf_df, argv[3], argv[4])
		# use barplot to visualize performance comparison across different method implementations on Tox21 datasets  
		dataset_annot_file = 'https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_assay_info.tsv'
		dataset_annot_df = pd.read_csv(dataset_annot_file, sep = '\t', header = 0)		
		method_compare = dtox_plot.visualize_dataset_performane_comparison(perf_file_df, dataset_annot_df, argv[5])
	
	return 1
 

## Call main function
main(sys.argv) 
