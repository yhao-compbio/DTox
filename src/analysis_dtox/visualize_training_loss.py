# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses line charts to visualize evolution of training/testing loss over epoches during DTox learning process 


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## 0. Input arguments 
optimal_file		= 'data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv'	# name of optimal DTox model performance (on Tox21 datasets) file  
tox21_annot_file        = '/home/yunhao1/project/tox_data/downloads/tox21/tox21_assay_info.tsv'	# name of Tox21 assay annotation file 
loss_folder		= 'data/compound_target_probability_tox21_implementation/'	# name of DTox model training/testing loss file folder   
out_folder		= 'plot/compound_target_probability_tox21_implementation/training_loss/compound_target_fingerprint_maccs_probability'	# name of output line chart folder  

## 1. Visualize evolution of training/testing loss over epoches during DTox learning process by line chart  
# read in optimal DTox model performance as data frame 
optimal_df = pd.read_csv(optimal_file, sep = '\t', header = 0)
# read in Tox21 assay annotation info as data frame   
tox21_annot_df = pd.read_csv(tox21_annot_file, sep = '\t', header = 0)
# iterate by Tox21 dataset  
for ods in range(0, optimal_df.shape[0]):
	# obtain name of current Tox21 dataset   
	ods_dataset = optimal_df.dataset_name[ods]
	ods_target = tox21_annot_df[tox21_annot_df.protocol_name == ods_dataset].assay_target.values[0].split(' (')[0]
	# add cell line name to dataset name if it is a viability assay (to differentiate two viability assays)
	if ods_target == 'Cell viability':
		ods_cell = tox21_annot_df[tox21_annot_df.protocol_name == ods_dataset].cell_line.values[0]
		ods_target = ods_target + ' (' + ods_cell + ')'
	# obtain training/testing loss file of current Tox21 dataset, read in training/testing loss info as data frame  
	ods_loss_file = loss_folder + ods_dataset.split('tox21-')[1] + '/compound_target_fingerprint_maccs_probability_' + ods_dataset + '_whole_data.tsv_' + optimal_df.hyperparameter_setting[ods] + '_loss.tsv' 
	ods_loss_df = pd.read_csv(ods_loss_file, sep = '\t', header = 0)	
	# obtain the training and testing loss at each epoch of learning process   
	train_loss = ods_loss_df.training_total_loss.values
	test_loss = ods_loss_df.testing_total_loss.values
	# obtain the epoch when optimal model was reached  
	optimal_point = ods_loss_df.shape[0] - 20
	# compute appropriate y-axis range for plotting  
	min_loss = np.min(np.concatenate([train_loss, test_loss]))
	max_loss = np.max(np.concatenate([train_loss, test_loss]))	
	if min_loss - np.floor(min_loss) > 0.5: 
		plt_min = np.floor(min_loss) + 0.5	
	else:
		plt_min = np.floor(min_loss)
	if plt_min == 0:
		plt_min = np.floor(min_loss * 10)/10
	if np.ceil(max_loss) - max_loss > 0.5:
		plt_max = np.ceil(max_loss) - 0.5
	else:
		plt_max = np.ceil(max_loss)
	# specify figure and font size 
	plt.figure(figsize = (6, 6))
	plt.rc('font', size = 25)
	plt.rc('axes', titlesize = 25)
	plt.rc('axes', labelsize = 25)
	plt.rc('xtick', labelsize = 22)
	plt.rc('ytick', labelsize = 22)
	plt.rc('legend', fontsize = 22)
	# make line chart showing the evolution of training loss over epoches  
	plt.plot(ods_loss_df.epoch.values, ods_loss_df.training_total_loss.values, '-', color = 'tab:blue', label = 'training')
	# make line chart showing the evolution of testing loss over epoches 
	plt.plot(ods_loss_df.epoch.values, ods_loss_df.testing_total_loss.values, '-', color = 'tab:orange', label = 'testing')
	# add red dashed line showing the epoch when optimal model was reached   
	plt.plot([optimal_point, optimal_point], [plt_min, plt_max], '--r', label = 'optimal')
	# set x-axis label and range 
	plt.xlabel('Epoch')
	plt.xlim([-5, 205])
	# set y-axis label and range 
	plt.ylabel('Loss function')
	plt.ylim([plt_min, plt_max])
	# add title (Tox21 dataset name) and legend   
	plt.title(ods_target)
	if ods == 0:
		plt.legend(loc = 'upper left', frameon = False)
	# save line chart	
	plt.tight_layout()
	plt.savefig(out_folder + '_' + ods_dataset + '_training_loss_plot.pdf')
	plt.close()

