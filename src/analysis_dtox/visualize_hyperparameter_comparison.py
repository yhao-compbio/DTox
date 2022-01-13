# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script uses heatmap and upsetplot to visualize normalized model performance of Tox21 datasets across root pathway settings 


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import upsetplot
from upsetplot import from_memberships
from upsetplot import plot


## 0. Input arguments
hp_data_file	= 'data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_rt_training_root_loss_normalized_comparison_by_dataset.tsv'	# name of input normalized model performance file  
heat_plot_file	= 'plot/compound_target_probability_tox21_implementation/hyperparameter_comparison/compound_target_probability_tox21_implementation_rt_training_root_loss_normalized_comparison_by_dataset.pdf'	# name of output heatmap plot file 
hp_onto_file	= 'https://raw.githubusercontent.com/yhao-compbio/ontology/master/data/reactome/root_file_map.tsv'	# name of input root pathway setting file  
upset_plot_file	= 'plot/compound_target_probability_tox21_implementation/hyperparameter_comparison/compound_target_probability_tox21_implementation_rt_upset.pdf' # name of output upsetplot file 

## 1. Visualize normalized model performance of Tox21 datasets across root pathway settings by heatmap
# read in normalized model performance of Tox21 datasets across root pathway settings as data frame  
hp_data_df = pd.read_csv(hp_data_file, sep = '\t', header = 0)
## specify figure and font size 
fig = plt.figure(figsize = (12, 6))
plt.rc('font', size = 20)
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)
# make heatmap showing the dataset ~ root pathway setting relationships
cg = sns.heatmap(hp_data_df,
	center = 0,
	cmap = 'PiYG',
	linewidths = 1,
	linecolor = 'black',
	xticklabels = False,
	cbar_kws = {'label': 'Normalized VNN performance'}
)
# save boxplot
fig.tight_layout()
plt.savefig(heat_plot_file)
plt.close()

## 2. Visualize pathway membership of root pathway settings by upsetplot 
# read in root pathway setting info as data frame  
hp_onto_df = pd.read_csv(hp_onto_file, sep = '\t', header = 0)
# obtain the pathway membership of each root pathway setting, convert it to upset format for plotting 
hp_onto_name = []
for i in range(0, hp_onto_df.shape[0]):
	i_onto = hp_onto_df.root_names.values[i].split(',')
	hp_onto_name.append(i_onto)
hp_upset = from_memberships(hp_onto_name)
# specify figure and font size  
fig = plt.figure(figsize = (10, 2))
plt.rc('font', size = 20)
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)
# make upsetplot showing the pathway membership of each root pathway setting  
plot(hp_upset, fig = fig, element_size = 32.5, intersection_plot_elements = 0)
# save boxplot
fig.tight_layout()
plt.savefig(upset_plot_file, bbox_inches = 'tight')
plt.close()

