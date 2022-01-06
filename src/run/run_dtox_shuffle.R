# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script generates scripts that run DTox model on Tox21 datasets under shuffeld Reactome pathway hierarchy


## functions
source("src/functions.R");


## 0. Input arguments
input_data_folder	<- "/home/yunhao1/project/TTox/data/compound_target_tox21_data/fingerprint_maccs_probability/compound_target_fingerprint_maccs_probability";	# folder name of input training-testing/validation data files 
hierarchy_folder	<- "/home/yunhao1/project/ontology/data/reactome/hierarchy/";	# folder name of input sorted reactome hierarchy files  
shuffle_folder		<- "/home/yunhao1/project/ontology/data/reactome/hierarchy_shuffle/shuffle/";	# folder name of input sorted shuffled reactome hierarchy files   
perf_file		<- "/home/yunhao1/project/DTox/data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";	# name of optimal model performance file 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome 
output_folder		<- "data/compound_target_probability_tox21_shuffle/compound_target_fingerprint_maccs_probability";	# folder name of output files 
N_cores			<- 15;	# number of CPUs
job_name		<- "dtox_compound_target_probability_tox21_shuffle";	# job name

## 1. Obtain optimal hyperparameter setting of trained model for each Tox21 dataset
# read in optimal performance info of Tox21 datasets as data frame  
perf_df <- read.delim(file = perf_file, header = T, sep = "\t");
# obtain hyperparameter setting of trained optimal model 
perf_hp <- mapply(function(pdhs){
	pdhs_s <- strsplit(pdhs, "_")[[1]];
	pdhs_rpr <- paste(pdhs_s[1:6], collapse = "_");
	pdhs_hp <- c(pdhs_rpr, pdhs_s[[4]], pdhs_s[[8]], pdhs_s[[10]], pdhs_s[[12]]);
	return(pdhs_hp);
}, perf_df$hyperparameter_setting);
perf_hp <- t(perf_hp);

## 2. Generate commands for jobs  
commands <- mapply(function(pddn, ph1, ph2, ph3, ph4, ph5){
	# combined training and testing dataset file 
	pddn_train <- paste(input_data_folder, pddn, "whole_data.tsv_train.tsv", sep = "_");
	# validation dataset file
	pddn_test <- paste(input_data_folder, pddn, "whole_data.tsv_test.tsv", sep = "_");
	# root pathway file 
	ph1_root <- paste(hierarchy_folder, ph1, "_st_0_root.tsv", sep = "");
	# shuffled parent/children node connection file 
	ph1_relation <- paste(shuffle_folder, ph1, "_st_0_knowledge_by_node_shuffle.tsv", sep = "");
	# node gene number file  
	ph1_size <- paste(hierarchy_folder, ph1, "_st_0_node_size.tsv", sep = "");
	# node layer number file  
	ph1_layer <- paste(hierarchy_folder, ph1, "_st_0_layer.tsv", sep = "");
	# output file
	pddn_output <- paste(output_folder, pddn, "whole_data.tsv", ph1, sep = "_");
	# put together command that includes files above and learned optimal hyperparameter setting 
	adfv_command <- paste("python", "src/dtox.py", pddn_train, pddn_test, outcome_col, ph1_root, ph1_relation, ph1_size, ph1_layer, ph2, pddn_output, ph3, ph4, ph5, sep = " ");
	return(adfv_command);
}, perf_df$dataset_name, perf_hp[,1], perf_hp[,2], perf_hp[,3], perf_hp[,4], perf_hp[,5]);
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
