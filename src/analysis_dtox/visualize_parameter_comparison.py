# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses barplot to visualize comparison of DTox and MLP model statistics across Tox21 datasets  


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments
parameter_file	= 'data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_model_parameter_summary.tsv'	# name of input model parameter stat file 
plot_file	= 'plot/compound_target_probability_tox21_implementation/parameter_comparison/compound_target_probability_tox21_implementation_optimal_model_parameter_summary_barplot.pdf'	# name of output barplot file 
 
## 1. Read in parameter stat of DTox and match MLP model as data frame  
parameter_df = pd.read_csv(parameter_file, sep = '\t', header = 0) 
# convert ratio to percentage
parameter_df.dtox_training_parameters_ratio = parameter_df.dtox_training_parameters_ratio * 100
parameter_df.dtox_mlp_parameters_ratio = parameter_df.dtox_mlp_parameters_ratio * 100

## 2. Specify plotting parameters 
# figure and font size 
fig, gs_ax = plt.subplots(1, 5, figsize = (24, 6), sharey = True)
plt.rc('font', size = 20)
# data columns to be plotted, column names, and x-axis ranges
plot_columns = ['n_training', 'n_pathway_module', 'n_dtox_parameters', 'dtox_training_parameters_ratio', 'dtox_mlp_parameters_ratio']
plot_xlabels = ['#Training samples', '#Pathway modules', '#VNN parameters', 'VNN samples/parameters %', 'VNN/MLP parameters %']
plot_xlims = [[0, 6500], [0, 610], [0, 65000], [0, 32], [0, 6]]

## 3. Visualize comparison of specified data columns across Tox21 datasets by barplot 
# iterate by data column
for i in range(0, len(plot_columns)):
	# make barplot showing the comparison of current data column across Tox21 datasets   
	sns.barplot(x = plot_columns[i], y = 'assay_name', data = parameter_df, color = 'blue', ax = gs_ax[i])
	# adjust width of each bar 
	new_h = 0.6
	for ap in gs_ax[i].patches:
		ap_y = ap.get_y()
		ap_h = ap.get_height()
		ap_c = ap_y + ap_h/2
		ap.set_y(ap_c - new_h/2)
		ap.set_height(new_h)
	# add red dashed line showing the mean value across all Tox21 datasets 
	col_mean = parameter_df[plot_columns[i]].mean()
	gs_ax[i].plot([col_mean, col_mean], [-0.5, parameter_df.shape[0]-0.5], '--r')
	# set x-axis range and label 
	gs_ax[i].set_xlim(plot_xlims[i])
	gs_ax[i].set_xlabel(plot_xlabels[i], fontsize = 20)
	# set y-axis range and label
	gs_ax[i].set_ylim([parameter_df.shape[0]-0.5, -0.5])
	gs_ax[i].set_ylabel(None)
	gs_ax[i].tick_params(axis = 'both', which = 'major', labelsize = 18)

## 4. Save barplot
sns.despine()
plt.tight_layout()
plt.savefig(plot_file)
plt.close()

