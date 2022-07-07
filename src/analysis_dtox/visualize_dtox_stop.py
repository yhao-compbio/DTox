# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script uses line charts to visualize the relative efficiency/performance of DTox under alternative settings of early stopping criterion 


## Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## 0. Input arguments 
stop_summary_file	= 'data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_early_stop_summary_by_patience.tsv'
plot_folder		= 'plot/compound_target_probability_tox21_implementation/early_stop/compound_target_fingerprint_maccs_probability_tox21_implementation_early_stop_summary_by_patience'	# name of output line chart folder  

## 1. Visualize relative efficiency/performance of DTox under alternative patience settings 
# read in computed relative efficiency/performance under alternative patience settings as data frame 
stop_summary_df = pd.read_csv(stop_summary_file, sep = '\t', header = 0)
# specify figure and font size 
fig = plt.figure(figsize = (12, 6))
plt.rc('font', size = 22)
plt.rc('axes', titlesize = 22)
plt.rc('axes', labelsize = 22)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)
# plot relative efficiency ~ alternative patience settings, with 95% confidence interval as error bar   
ax1 = fig.add_subplot(121)
run_time_error = stop_summary_df.loc[:, ['run_time_ratio_95ci_lower', 'run_time_ratio_95ci_upper']].T.values - stop_summary_df.run_time_ratio.values
run_time_error = np.absolute(run_time_error) * 100
ax1_x = stop_summary_df.patience.values 
ax1_y = stop_summary_df.run_time_ratio.values * 100
ax1.errorbar(ax1_x, ax1_y, yerr = run_time_error, marker = 'o', capsize = 4, color = 'tab:cyan')
# set axis labels, ticks, ranges 
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.set_xlabel('Patience (default = 20)')
ax1.set_ylabel('Model running time (% of default)')
ax1.set_xticks([1, 5, 10, 15, 20])
ax1.set_ylim(0, 102)
# plot relative performance ~ alternative patience settings, with 95% confidence interval as error bar 
ax2 = fig.add_subplot(122)
perf_error = stop_summary_df.loc[:, ['performance_ratio_95ci_lower', 'performance_ratio_95ci_upper']].T.values - stop_summary_df.performance_ratio.values
perf_error = np.absolute(perf_error) * 100
ax2_x = stop_summary_df.patience.values
ax2_y = stop_summary_df.performance_ratio.values * 100
ax2.errorbar(ax2_x, ax2_y, yerr = perf_error, marker = 'o', capsize = 4, color = 'tab:pink')
# set axis labels, ticks, ranges 
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Patience (default = 20)')
ax2.set_ylabel('Model performance (% of default)')
ax2.set_xticks([1, 5, 10, 15, 20])
ax2.set_ylim(70, 102)
# save barplot      
plt.tight_layout()
plt.savefig(plot_folder + '_line_chart.pdf')
plt.close()

