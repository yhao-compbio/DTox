# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script uses boxplot to visualize the distributions of predicted HEK293 cytotoxicity scores of compounds in EPA/DrugBank lists


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments 
score_file	= 'data/compound_target_probability_tox21_prediction/rt-viability-hek293-p1/dsstox_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hek293-p1_whole_data.tsv_rt_14_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred_by_list.tsv'
el_compare_file	= 'data/compound_target_probability_tox21_prediction/rt-viability-hek293-p1/dsstox_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hek293-p1_whole_data.tsv_rt_14_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred_epa_list_compare_result_summary.tsv'
dl_compare_file	= 'data/compound_target_probability_tox21_prediction/rt-viability-hek293-p1/dsstox_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hek293-p1_whole_data.tsv_rt_14_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred_drugbank_list_compare_result_summary.tsv'
plot_file	= 'plot/compound_target_probability_tox21_prediction/epa_list_tox21-rt-viability-hek293-p1_dsstox_compound_model_pred_score_compare'	# folder name of output plof files 
list_name_dict	= {'tbutylphenols': 'butylphenols', 'wikiantifungals': 'antifungals', 'npinsect': 'insecticides', 'iarc1': 'carcinogens'}

## 1. Process input comparison and prediction data
# read in predicted HEK293 cytotoxicity scores of compounds in EPA/DrugBank lists 
score_df = pd.read_csv(score_file, sep = '\t', header = 0) 
sd_list_names = score_df.list_name.unique()
# read in the result comparing each EPA list vs pos/neg controls  
el_compare_df = pd.read_csv(el_compare_file, sep = '\t', header = 0) 
el_compare_df = el_compare_df[el_compare_df.list_acronym.isin(sd_list_names)]
# read in the result comparing each DrugBank list vs pos/neg controls 
dl_compare_df = pd.read_csv(dl_compare_file, sep = '\t', header = 0)
dl_compare_df = dl_compare_df[dl_compare_df.list_acronym.isin(sd_list_names)]
# sort the lists by median of predicted HEK293 cytotoxicity scores
ecd_order = el_compare_df.sort_values("compound_pred_score_median", ascending = False).list_acronym.values
dcd_order = dl_compare_df.sort_values("compound_pred_score_median", ascending = False).list_acronym.values
all_order = np.concatenate([['HEK293 +'], ecd_order, dcd_order, ['HEK293 -']]) 
# specify plotting colors  
ecd_col = np.repeat('mistyrose', len(ecd_order))
dcd_col = np.repeat('lightblue', len(dcd_order))
all_col = np.concatenate([['red'], ecd_col, dcd_col, ['blue']])
ecd_jc = np.repeat('darkred', len(ecd_order))
dcd_jc = np.repeat('darkblue', len(dcd_order))
all_jc = np.concatenate([['red'], ecd_jc, dcd_jc, ['blue']])
# specify jitter styles 
score_df1 = score_df.copy()
score_df1.loc[score_df1.list_name == 'HEK293 +', 'predicted_cytotoxicity_score'] = np.nan
score_df1.loc[score_df1.list_name == 'HEK293 -', 'predicted_cytotoxicity_score'] = np.nan
# specify list name
all_name = []
for ao in all_order:
	if ao in list_name_dict.keys():
		all_name.append(list_name_dict[ao])
	else:
		all_name.append(ao)

## 2. Visualize distributions of predicted HEK293 cytotoxicity scores of compounds in EPA/DrugBank lists   
# specify figure and font size of boxplot plot 
plt.figure(figsize = (16, 8))
plt.rc('font', size = 25)
plt.rc('axes', titlesize = 25)
plt.rc('axes', labelsize = 25)
plt.rc('xtick', labelsize = 25)
plt.rc('ytick', labelsize = 25)
plt.rc('legend', fontsize = 25)
# make boxplot showing the distrubtion of predicted HEK293 cytotoxicity scores of compounds in EPA/DrugBank lists 
ax = sns.boxplot(x = 'list_name', y = 'predicted_cytotoxicity_score', data = score_df, showfliers = False, palette = all_col, order = all_order, width = 0.6, linewidth = 1.5)
ax = sns.stripplot(x = 'list_name', y = 'predicted_cytotoxicity_score', data = score_df1, palette = all_jc, order = all_order, jitter = 0.25, size = 2.5)
# add red asterisk under EPA list with no significant difference from pos control, blue asterisk under DrugBank list with no significant difference from neg control
sig_y = -0.05
for lao in range(0, len(all_order)):
	ao = all_order[lao]
	if ao in el_compare_df.list_acronym.values:
		ao_pv = el_compare_df[el_compare_df.list_acronym.isin([ao])].pos_compare_wilcox_pv.values
		if ao_pv > 0.05:
			plt.plot(lao, sig_y, marker = '*', c = 'r', linestyle = 'None', markersize = 15, clip_on = False)
	elif ao in dl_compare_df.list_acronym.values:
		ao_pv = dl_compare_df[dl_compare_df.list_acronym.isin([ao])].neg_compare_wilcox_pv.values
		if ao_pv > 0.05:
			plt.plot(lao, sig_y, marker = '*', c = 'b', linestyle = 'None', markersize = 15, clip_on = False)
	else:
		continue
# set axis label
ax.set_xticklabels(all_name, rotation = 90)
ax.set(ylabel = 'Predicted cytotoxicity score')
ax.set(xlabel = None)
# save boxplot
plt.tight_layout()
sns.despine()
plt.savefig(plot_file + '_boxplot.pdf')
plt.close()

