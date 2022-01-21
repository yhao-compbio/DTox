# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses boxplot and barplot to visualize the comparison of DTox HEK293 viability prediction results between positive and negative compounds of drug-induced kidney injury (DIKI) phenotypes


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments 
score_file	= 'data/compound_target_probability_tox21_prediction/rt-viability-hek293-p1/diki_phenotype_tox21-rt-viability-hek293-p1_dsstox_compound_model_pred_score.tsv'	# name of input DSSTox compound HEK293 viability prediction result file 
compare_file	= 'data/compound_target_probability_tox21_prediction/rt-viability-hek293-p1/diki_phenotype_tox21-rt-viability-hek293-p1_dsstox_compound_model_pred_score_compare.tsv'	# name of input DIKI phenotype comparison file
plot_file	= 'plot/compound_target_probability_tox21_prediction/diki_phenotype_tox21-rt-viability-hek293-p1_dsstox_compound_model_pred_score_compare'	# folder name of output plof files 

## 1. Process input comparison and prediction data 
# read in predicted cytotoxicity probability comparison result between DIKI + vs - compounds as data frame 
compare_df = pd.read_csv(compare_file, sep = '\t', header = 0)
# sort all DIKI phenotypes by HEK293 viability-association odds ratio (OR) in decreasing order, identify DIKI phenotypes with OR > 1  
compare_df1 = compare_df[compare_df.viability_odds_ratio > 1]
ae_order1 = compare_df1.adverse_event.values
# read in predicted cytotoxicity probability of DSSTox compounds as data frame, select rows relevant to identified DIKI phenotypes 
score_df = pd.read_csv(score_file, sep = '\t', header = 0)
score_df1 = score_df[score_df.adverse_event.isin(ae_order1)]

## 2. Visualize comparison of predicted cytotoxicity probability within DIKI phenotypes by boxplot 
# specify figure and font size of boxplot plot 
plt.figure(figsize = (10, 7))
plt.rc('font', size = 25)
plt.rc('axes', titlesize = 25)
plt.rc('axes', labelsize = 25)
plt.rc('xtick', labelsize = 25)
plt.rc('ytick', labelsize = 25)
plt.rc('legend', fontsize = 25)
# make boxplot showing the comparison of predicted cytotoxicity probability between positive and negative compounds of each DIKI phenotype 
ax = sns.boxplot(x = 'predicted_outcome', y = 'adverse_event', hue = 'toxicity_label', data = score_df1, order = ae_order1, hue_order = ['+', '-'], showfliers = False, palette = 'Set3', width = 0.7, linewidth = 1.5)
# add asterisk to show DIKI phenotypes with significant comparison P-value  
sig_x = []
sig_y = []
for cds in range(0, compare_df1.shape[0]):
	cdpcpv = compare_df1.prediction_compare_greater_pv.values[cds]
	if cdpcpv < 0.05:
		sig_x.append(-0.05)
		sig_y.append(cds + 0.4)
plt.plot(sig_x, sig_y, marker = '*', color = 'r', linestyle = 'None', markersize = 15, clip_on = False)
# set range, label, and ticks of axes   
ax.set_xlim([-0.02, 1])
ax.set(xlabel = 'Predicted cytotoxicity score')
ax.set_ylim([compare_df1.shape[0] - 0.5, -0.5])
ax.set(ylabel = None)
plt.yticks(np.arange(len(ae_order1)), ae_order1)
plt.legend(loc = 'upper right', frameon = False, bbox_to_anchor=(1.25, 0.3))
# save boxplot
plt.tight_layout()
sns.despine()
plt.savefig(plot_file + '_boxplot_greater.pdf')
plt.close()

## 3. Visualize HEK293 viability-association odds ratio of DIKI phenotypes by barplot 
# specify figure and font size of boxplot plot 
plt.figure(figsize = (6, 7))
plt.rc('font', size = 25)
plt.rc('axes', titlesize = 25)
plt.rc('axes', labelsize = 25)
plt.rc('xtick', labelsize = 25)
plt.rc('ytick', labelsize = 25)
plt.rc('legend', fontsize = 25)
# make barplot showing the HEK293 viability-association odds ratio of DIKI phenotypes   
ax = sns.barplot(x = 'viability_odds_ratio', y = 'adverse_event', data = compare_df1, color = 'salmon')
# adjust width of each bar
new_h = 0.5
for ap in ax.patches:
	ap_y = ap.get_y()
	ap_h = ap.get_height()
	ap_c = ap_y + ap_h/2
	ap.set_y(ap_c - new_h/2)
	ap.set_height(new_h)
# add error bars to show the 95% confidence interval of OR for each DIKI phenotype 
for cd in range(0, compare_df1.shape[0]):
	cd_x1 = compare_df1.viability_odds_ratio_95ci_lower.values[cd]
	cd_x2 = compare_df1.viability_odds_ratio_95ci_upper.values[cd]
	ax.plot([cd_x1, cd_x2], [cd, cd], color = 'black')
# add red dashed line to show x coordinate at OR = 1  
ax.plot([1, 1], [-0.5, compare_df1.shape[0] - 0.5], '--r')
# set range, label, and ticks of axes  
ax.axes.get_yaxis().set_visible(False)	
ax.set_xlim([-0.02, 3])
ax.set(xlabel = 'Odds ratio w. cytotoxicity')
ax.set_ylim([compare_df1.shape[0] - 0.5, -0.5])
ax.set(ylabel = None)
ax.set(yticklabels = [])
plt.gca().invert_xaxis()
# save boxplot
plt.tight_layout()
sns.despine(left = True, right = False)
plt.savefig(plot_file + '_barplot_greater.pdf')
plt.close()
