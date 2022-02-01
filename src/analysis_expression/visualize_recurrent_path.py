# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses barplot to visualize the frequency of recurrent differentially expressed DTox paths from model interpretation results on Tox21 dataset of interest  


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments
recurrent_file	= 'data/compound_target_probability_tox21_interpret_expression/validation_summary/gamma-epsilon_0.001_0.1_expression_validation_sig_path_recurrent.tsv'	# name of input recurrent differentially expressed DTox path file 
plot_file	= 'plot/compound_target_probability_tox21_interpret_expression/gamma-epsilon_0.001_0.1_expression_validation_sig_path_recurrent'	# name of output barplot file
query_assay	= 'tox21-aromatase-p1'	# string that indicates Tox21 assay of interest   
recurrent_cut	= 5	# number that indicates occurrence threshold for a differentially expressed DTox path to be included in barplot 

## 1. Process input recurrent differentially expressed DTox path file 
# read in recurrent differentially expressed DTox path info as data frame  
recurrent_df = pd.read_csv(recurrent_file, sep = '\t', header = 0) 
# select DTox path rows of interest based on specified Tox21 assay and occurrence threshold  
query_recurrent_df = recurrent_df[recurrent_df.assay == query_assay]
query_recurrent_df = query_recurrent_df[query_recurrent_df.n_compounds >= recurrent_cut]
# iterate by selected row, obtain the name of Reactome pathways at the end of each DTox path (lowest Reacomte hierarchy) 
query_recurrent_lowest_path = []
for qrds in range(0, query_recurrent_df.shape[0]):
	qrds_path_s = query_recurrent_df.recurrent_path.values[qrds].split(';')
	qrds_lp = qrds_path_s[-1]
	query_recurrent_lowest_path.append(qrds_lp)
query_recurrent_df['lowest_path'] = query_recurrent_lowest_path

## 2. Specify figure and font size of barplot 
plt.figure(figsize = (15, 5.5))
plt.rc('font', size = 35)
plt.rc('axes', titlesize = 35)
plt.rc('axes', labelsize = 35)
plt.rc('xtick', labelsize = 30)
plt.rc('ytick', labelsize = 30)
plt.rc('legend', fontsize = 30)

## 3. Make barplot
# iterate by DTox path to be plotted, assign plotting color to bar of each DTox path 
path_color = []
for qrlp in query_recurrent_lowest_path:
	# red if the path is related to TP53, lightblue otherwise
	if ('TP53' in qrlp) | ('TFAP2' in qrlp):
		path_color.append('salmon')
	else:
		path_color.append('lightblue')
# make barplot showing the occurrence frequency of selected differentially expressed DTox paths, which are named after the end pathways along each DTox path  
ax = sns.barplot(x = 'n_compounds', y = 'lowest_path', data = query_recurrent_df, palette = path_color)
# adjust width of each bar 
new_h = 0.5
for ap in ax.patches:
	ap_y = ap.get_y()
	ap_h = ap.get_height()
	ap_c = ap_y + ap_h/2
	ap.set_y(ap_c - new_h/2)
	ap.set_height(new_h)
# set figure title and axis labels  
ax.set(title = 'Recurrent DE pathways in aromatase assay')
ax.set(xlabel = 'Number of compounds')
ax.set(ylabel = None)

## 4. Save barplot 
sns.despine(left = True, right = True)
plt.savefig(plot_file + '_' + query_assay + '_' + str(recurrent_cut) + '.pdf', bbox_inches = 'tight')
plt.close()
