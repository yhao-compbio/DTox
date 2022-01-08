# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that run fully connected MLP neural network models on Tox21 datasets, which are built with the same number of hidden layer/neuron as matched DTox models  


## functions
source("src/functions.R");


## 0. Input arguments
Args			<- commandArgs(T);
input_optimal_file	<- Args[1];		# name of optimal DTox model performance file  
input_data_folder	<- Args[2];		# folder name of input training-testing/validation data files 
config_folder		<- Args[3];		# folder name of trained DTox model configuration/performance files
output_folder		<- Args[4];		# folder name of output files   
N_cores			<- Args[5];		# number of CPUs
job_name		<- Args[6];		# job name 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome  

## 1. Process input dataset and DTox model configuration files 
# read in optimal DTox model performance file (contains Tox21 dataset name, optimal hyperparameter setting for each dataset, and other info)  
optimal_df <- read.delim(file = input_optimal_file, sep = "\t", header = T);
# obtain disk location of each Tox21 dataset included in the DTox optimal performance data frame 
input_data_files <- sapply(optimal_df$dataset_name, function(oddn) paste(input_data_folder, oddn, "whole_data.tsv", sep = "_"));
whole_data_files <- sapply(input_data_files, function(idf){
	idf_s <- strsplit(idf, "/", fixed = T)[[1]];
	is_len <- length(idf_s);
	return(idf_s[[is_len]]);
});
# obtain disk location of DTox optimal model configuration/performance file for each Tox21 dataset (named after optimal hyperparameter setting in optimal_df) 
config_files <- mapply(function(oddn, wdf, odhs){
	oddn_s <- strsplit(oddn, "tox21-")[[1]][[2]];
	oddn_config <- paste(config_folder, oddn_s, "/", wdf, "_", odhs, "_performance.txt", sep = "");
	return(oddn_config);
}, optimal_df$dataset_name, whole_data_files, optimal_df$hyperparameter_setting);
# generate disk location of output file for each Tox21 dataset  
output_files <- sapply(whole_data_files, function(wdf) paste(output_folder, wdf, "_mlp_fully_connected", sep = ""));

## 2. Generate commands for jobs  
commands <- mapply(function(idf, cf, of){
	# combined training-testing dataset file
	idf_train <- paste(idf, "_train.tsv", sep = "");
	# validation dataset file  
	idf_test <- paste(idf, "_test.tsv", sep = "");
	# put together command that includes files above and other specified files  
	idf_command <- paste("python", "src/mlp/mlp.py", idf_train, idf_test, outcome_col, cf, of, sep = " ");
	return(idf_command);
}, input_data_files, config_files, output_files);
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
