# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script uses barplot to visualize the frequency of target proteins and lowest level pathways among the cytotoxic compounds not linked to the viability-related pathways by DTox  


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


## 0. Input arguments
new_path_file	= 'data/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_compound_new_pathways_summary.tsv'
new_target_file	= 'data/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_compound_new_targets_summary.tsv'
N_plot		= 30
plot_folder	= 'plot/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_compound'
path_dict = {'Class C/3 (Metabotropic glutamate/pheromone receptors)': 'Metabotropic glutamate/pheromone receptors', 'Orexin and neuropeptides FF and QRFP bind to their respective receptors': 'Orexin and neuropeptides FF and QRFP bind to receptors', 'Activated PKN1 stimulates transcription of AR (androgen receptor) regulated genes KLK2 and KLK3': 'Activated PKN1 stimulates transcription of KLK2 and KLK3', 'Formyl peptide receptors bind formyl peptides and many other ligands': 'Formyl peptide receptors bind formyl peptides and others', 'GRB2:SOS provides linkage to MAPK signaling for Integrins ': 'GRB2:SOS provides linkage to MAPK for Integrins', 'Inactivation, recovery and regulation of the phototransduction cascade': 'Inactivation, recovery and regulation of phototransduction'}

## 1. Use barplot to visualize the frequency of lowest level pathways among the cytotoxic compounds not linked to the viability-related pathways by DTox
# read in the data frame that contains frequency of lowest level pathways among such cytotoxic compounds
new_path_df = pd.read_csv(new_path_file, sep = '\t', header = 0)
new_path_df = new_path_df.iloc[0:N_plot, :]
new_path_cates = new_path_df.general_pathway_name.unique()
# specifiy color palette for plotting, one color for general pathway category 
all_colors = sns.color_palette("tab10")
# iterate by lowest level pathways of interest  
new_path_colors = []
new_path_values = []
new_path_names = []
for i in range(0, N_plot):
	# specify the bar color based on the general pathway category of current lowest level pathway 
	i_cate = new_path_df.general_pathway_name[i]
	i_id = np.where(new_path_cates == i_cate)[0][0]
	new_path_colors.append(all_colors[i_id])
	# create text containing prevalence percentage 
	i_value = str(new_path_df.path_compound_count[i]) + '(' + str(int(round(new_path_df.path_compound_ratio[i], 3) * 100)) + '%)'
	new_path_values.append(i_value)
	# abbreviate name of lowest level pathway if it is too long 
	i_path = new_path_df.lowest_level_pathway_name[i]
	if i_path in path_dict.keys():
		new_path_names.append(path_dict[i_path])
	else:
		new_path_names.append(i_path)
# specify figure and font size of barplot 
plt.figure(figsize = (15, 20))
plt.rc('font', size = 25)
plt.rc('axes', titlesize = 30)
plt.rc('axes', labelsize = 30)
plt.rc('xtick', labelsize = 30)
plt.rc('ytick', labelsize = 30)
plt.rc('legend', fontsize = 25)
# make barplot showing the prevalence percentage of lowest level pathways among such cytotoxic compounds 
ax = sns.barplot(x = 'path_compound_count', y = 'lowest_level_pathway_name', data = new_path_df, palette = new_path_colors)
# adjust width of each bar 
new_h = 0.5
for ap in ax.patches:
	ap_y = ap.get_y()
	ap_h = ap.get_height()
	ap_c = ap_y + ap_h/2
	ap.set_y(ap_c - new_h/2)
	ap.set_height(new_h)
# add prevalence percentage as text along the bar
for i in range(0, N_plot):
	i_x = new_path_df.path_compound_count[i]
	ax.text(i_x, i+0.2, new_path_values[i], color = 'black', ha = 'left')
# set figure title and axis labels  
ax.set(title = 'Prevalent target proteins')
ax.set(xlabel = 'Number of HepG2-cytotoxic compounds')
ax.set(ylabel = None)
ax.set_yticklabels(new_path_names)
# add figure legend specifying the color labels (general pathway categories)  
legend_elements = []
for lpnc in range(0, len(new_path_cates)):
	legend_elements.append(Patch(facecolor = all_colors[lpnc], edgecolor = all_colors[lpnc], label = new_path_cates[lpnc]))
ax.legend(handles = legend_elements, loc = 'upper right', bbox_to_anchor = (1.1, 1.4), frameon = False)
# save barplot 
sns.despine(left = True, right = True)
plt.savefig(plot_folder + '_new_pathways_barplot.pdf', bbox_inches = 'tight')
plt.close()

## 2. Use barplot to visualize the frequency of target proteins among the cytotoxic compounds not linked to the viability-related pathways by DTox 
# read in the data frame that contains frequency of target proteins among such cytotoxic compounds 
new_tar_df = pd.read_csv(new_target_file, sep = '\t', header = 0)
new_tar_df = new_tar_df.iloc[0:N_plot, :]
# iterate by target proteins of interest
new_tar_colors = []
new_tar_values = []
for i in range(0, N_plot):
	# specify the bar color based on the general pathway category of current target protein 
	i_cate = new_tar_df.general_pathway_name[i]
	i_id = np.where(new_path_cates == i_cate)[0][0]
	new_tar_colors.append(all_colors[i_id])
	# create text containing prevalence percentage
	i_value = str(new_tar_df.path_compound_count[i]) + '(' + str(int(round(new_tar_df.path_compound_ratio[i], 3) * 100)) + '%)'
	new_tar_values.append(i_value)
# specify figure and font size of barplot 
plt.figure(figsize = (14, 20))
plt.rc('font', size = 25)
plt.rc('axes', titlesize = 30)
plt.rc('axes', labelsize = 30)
plt.rc('xtick', labelsize = 30)
plt.rc('ytick', labelsize = 30)
plt.rc('legend', fontsize = 25)
# make barplot showing the prevalence percentage of target proteins among such cytotoxic compounds 
ax = sns.barplot(x = 'path_compound_count', y = 'target_protein_name', data = new_tar_df, palette = new_tar_colors)
# adjust width of each bar 
new_h = 0.5
for ap in ax.patches:
	ap_y = ap.get_y()
	ap_h = ap.get_height()
	ap_c = ap_y + ap_h/2
	ap.set_y(ap_c - new_h/2)
	ap.set_height(new_h)
# add prevalence percentage as text along the bar 
for i in range(0, N_plot):
	i_x = new_tar_df.path_compound_count[i]
	ax.text(i_x, i+0.2, new_tar_values[i], color = 'black', ha = 'left')
# set figure title and axis labels  
ax.set(title = 'Prevalent target proteins')
ax.set(xlabel = 'Number of HepG2-cytotoxic compounds')
ax.set(ylabel = None)
# save barplot 
sns.despine(left = True, right = True)
plt.savefig(plot_folder + '_new_targets_barplot.pdf', bbox_inches = 'tight')
plt.close()

