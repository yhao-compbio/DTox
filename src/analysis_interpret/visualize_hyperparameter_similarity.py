# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses heatmap to visualize the similarity of significant DTox paths under different hyperparameter settings of layer-wise relevance propagation rule on Tox21 datasets 


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments
sim_data_file	= 'data/compound_target_probability_tox21_interpret_analysis/compound_target_fingerprint_maccs_probability_gamma-epsilon_path_similarity.tsv'	# name of input similarity index file 
sim_plot_file	= 'plot/compound_target_probability_tox21_interpret_analysis/compound_target_fingerprint_maccs_probability_gamma-epsilon_path_similarity.pdf'	# name of output heatmap file 

## 1. Process input matrix that contains Jaccard Index of significant DTox paths between hyperparameter setting pairs 
# read in symmetric matrix that contains Jaccard Index of significant DTox paths between hyperparameter setting pairs as data frame  
sim_data_df = pd.read_csv(sim_data_file, sep = '\t', header = 0)
# remove the redundant part of upper triangle of symmetric matrix 
sim_mask = np.triu(np.ones_like(sim_data_df, dtype = np.bool))
# iterate by hyperparameter setting (values of γ and ε), name each column and row by values of γ and ε 
sim_labels = []
for sc in sim_data_df.columns:
	sc_s = sc.split('_')
	sc_l = 'γ = ' + sc_s[1] + '\nε = ' + sc_s[2]
	sim_labels.append(sc_l)

## 2. Specify figure and font size of heatmap  
fig = plt.figure(figsize = (11, 9))
plt.rc('font', size = 18)
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)

## 3. Make heatmap showing the similarity of significant DTox paths between hyperparameter setting pairs 
ax = sns.heatmap(sim_data_df,
	vmin = 0, 
	vmax = 1,
	cmap = 'rocket_r', 
	annot = True,
	mask = sim_mask,
	xticklabels = sim_labels,
	yticklabels = sim_labels,
	cbar_kws = {'label': 'Jaccard index'}
)
ax.set_title('Similarity of identified VNN paths', fontdict = {'fontsize': 25});

## 4. Save boxplot
fig.tight_layout()
plt.savefig(sim_plot_file)
plt.close()


