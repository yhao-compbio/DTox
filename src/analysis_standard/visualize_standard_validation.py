# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses density plot and barplot to visualize the standard pattern-validation of DTox interpretation results on Tox21 datasets, comparing the observed and expected proportion of validated compounds   


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm


## 0. Input arguments
sample_file	= 'data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_standard_validation_bg_sample.tsv'	# name of input sampled expected proportion of validated compounds file  
compare_file	= 'data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_standard_validation_result_compare_summary.tsv'	# name of input observed/expected proportion stats file  
assay_file	= 'https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_assay_reactome_map.tsv'	# name of Tox21 assay to Reactome pathway-receptor file 
plot_file	= 'plot/compound_target_probability_tox21_interpret_standard/compound_target_fingerprint_maccs_probability_tox21_interpret_standard_validation_result_compare_'	# folder name of output plot files 
query_rule	= 'gamma-epsilon_0.001_0.1'	# string that indicates query LPR propagation rule along with parameter values for extracting validation results of interest

## 1. Process input DTox interpretation-standard validation files  
# read in sampled expected proportion of validated compounds as data frame, select rows of the data frame relevant to the specified LRP propagation rule  
sample_df = pd.read_csv(sample_file, sep = '\t', header = 0)
sample_df = sample_df[sample_df.rule == query_rule]
# read in observed/expected proportion stats of validated compounds as data frame, select rows of the data frame relevant to the specified LRP propagation rule 
compare_df = pd.read_csv(compare_file, sep = '\t', header = 0)
compare_df = compare_df[compare_df.rule == query_rule]

## 2. Visualize standard pattern-validation results of Tox21 assay datasets by density plots and barplots 
# read in Tox 21 assay to standard Reactome pathway-receptor pattern map as data frame  
assay_df = pd.read_csv(assay_file, sep = '\t', header = 0) 
# iterate by Tox21 assay, visualize standard pattern-validation results of each assay   
all_assays = compare_df.assay.values
for aa in all_assays:
	# specify figure and font size of density plot and barplot  
	fig = plt.figure(figsize = (6, 3))
	plt.rc('font', size = 20)
	plt.rc('axes', titlesize = 20)
	plt.rc('axes', labelsize = 20)
	plt.rc('xtick', labelsize = 20)
	plt.rc('ytick', labelsize = 20)
	plt.rc('legend', fontsize = 20)
	gs = fig.add_gridspec(2, 1, hspace = 0, wspace = 0)
	gs_ax = gs.subplots(sharex = True)
	# select sampled expected proportion of validated compounds relevant to the current assay, make density plot showing the distribution of sampled expected proportions  
	aa_sample_df = sample_df[sample_df.assay == aa]
	sns.distplot(aa_sample_df.bg_sample.values, hist = True, fit = norm, kde = False, bins = 20, color = 'darkblue', ax = gs_ax[0])	
	# set figure title (Tox21 assay name) and x-axis range (0-1) 
	gs_ax[0].set_xlim(0, 1)
	gs_ax[0].axis('off')
	aa_target_name = assay_df[assay_df.assay == aa].assay_target.values[0]
	gs_ax[0].set_title(aa_target_name, loc = 'right', color = 'red')
	# select observed/expected proportion stat of validated compounds relevant to the current assay, make barplot showing the comparion between observed and expected proportion 
	aa_compare_df = compare_df[compare_df.assay == aa]
	bar_values1 = [aa_compare_df.standard_observed_ratio.values[0], aa_compare_df.standard_bg_med.values[0]]
	bar_values2 = 1 - np.array(bar_values1)
	bar_pos = ['Observed', 'Expected']
	gs_ax[1].barh(bar_pos, bar_values1, height = 0.6, color = 'tab:purple', edgecolor = 'black', linewidth = 1.2)
	gs_ax[1].barh(bar_pos, bar_values2, height = 0.6, left = bar_values1, color = 'white', edgecolor = 'black', linewidth = 1.2)
	# add error bar showing the 95% confidence interval of expected proportion  
	bar_err_x = [aa_compare_df.standard_bg_95ci_lower.values[0], aa_compare_df.standard_bg_95ci_upper.values[0]]
	bar_err_y = ['Expected', 'Expected'] 
	gs_ax[1].plot(bar_err_x, bar_err_y, linewidth = 1.5, color = 'k')	
	# set x-axis range (0-1)
	gs_ax[1].set_xlim(0, 1)
	gs_ax[1].set_xlabel('Validated compounds %')
	gs_ax[1].spines['top'].set_visible(False)
	gs_ax[1].spines['right'].set_visible(False)
	# save barplot
	plt.tight_layout()
	plt.savefig(plot_file + query_rule + '_' + aa + '.pdf')
	plt.close()
