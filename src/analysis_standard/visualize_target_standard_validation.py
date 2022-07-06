# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022 
## This script uses line charts to visualize the validation performance comparison among DTox, LIME, and Read-across regarding the interpretation task of connecting active compounds to their respective target receptor in four nuclear receptor assays


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm


## 0. Input arguments
compare_file	= 'data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_target_standard_validation_result_compare_summary.tsv'	
assay_file	= 'https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_assay_reactome_map.tsv'	# name of Tox21 assay to Reactome pathway-receptor file 
plot_file	= 'plot/compound_target_probability_tox21_interpret_standard/compound_target_fingerprint_maccs_probability_tox21_interpret_target_standard_validation_result_compare'	# folder name of output plot files 

## 1. Read in validation performance comparison file 
compare_df = pd.read_csv(compare_file, sep = '\t', header = 0)

## 2. Visualize validation performance comparison by line charts  
# read in Tox 21 assay to standard Reactome pathway-receptor pattern map as data frame  
assay_df = pd.read_csv(assay_file, sep = '\t', header = 0) 
all_assays = assay_df.assay.values
# specify figure and font size of density plot and barplot  
fig = plt.figure(figsize = (3*len(all_assays), 3))
plt.rc('font', size = 13.5)
plt.rc('axes', titlesize = 13.5)
plt.rc('axes', labelsize = 13.5)
plt.rc('xtick', labelsize = 13.5)
plt.rc('ytick', labelsize = 13.5)
plt.rc('legend', fontsize = 11)
# iterate by assays of interest  
for laa in range(0, len(all_assays)):
	# validation performance comparison metrics of current assay 
	aa = all_assays[laa]
	aa_compare_df = compare_df[compare_df.assay == aa]
	# creat subplot  
	fig_n = int(str(1) + str(len(all_assays)) + str(laa + 1))
	ax = fig.add_subplot(fig_n)
	# make line chart plotting proportion of validated compounds ~ TC (similarity) threshold of Read-across 
	acd_x = aa_compare_df.tc_threshold.values 
	ax.plot(acd_x, aa_compare_df.standard_observed_ratio_drugbank.values*100, marker = '.', label = 'RAx (DrugBank)', color = 'tab:orange')
	ax.plot(acd_x, aa_compare_df.standard_observed_ratio_comptoxai.values*100, marker = '.', label = 'RAx (ComptoxAI)', color = 'tab:pink')
	ax.plot(acd_x, aa_compare_df.standard_observed_ratio_lime_strict.values*100, marker = '.', label = 'LIME (strict)', color = 'tab:blue')
	ax.plot(acd_x, aa_compare_df.standard_observed_ratio_lime_loose.values*100, marker = '.', label = 'LIME (lax)', color = 'tab:cyan')
	ax.plot(acd_x, aa_compare_df.standard_observed_ratio_dtox.values*100, marker = '.', label = 'DTox', color = 'tab:red')
	# add axis labels and ticks 
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('TC threshold for RAx')
	ax.set_ylabel('Validated compounds %')
	ax.set_ylim(-5, 100)
	ax.set_yticks(np.arange(0, 101, 25))
	# add main title  
	aa_tar_name = assay_df[assay_df.assay == aa].assay_target.values[0]
	aa_title = aa_tar_name.split('(')[1].split(')')[0]
	ax.set_title(aa_title, color = 'red')
	# add figure legend specifying method labels  
	if laa == 0:
		ax.legend(frameon = False, loc = 'upper right', bbox_to_anchor = (1, 1.05))
# save line chart 
plt.tight_layout()
plt.savefig(plot_file + '_line_chart.pdf')
plt.close()
