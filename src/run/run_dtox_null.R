# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that run DTox model on outcome-shuffled Tox21 datasets under Reactome pathway hierarchy


## functions
source("src/functions.R");


## 0. Input arguments
train_data_folder	<- "/home/yunhao1/project/TTox/data/compound_target_tox21_data/fingerprint_maccs_probability_shuffle/";	# folder name of input outcome-shuffled training-testing data files 
test_data_folder	<- "/home/yunhao1/project/TTox/data/compound_target_tox21_data/fingerprint_maccs_probability/compound_target_fingerprint_maccs_probability";	# folder name of input validation data files 	
hierarchy_folder	<- "/home/yunhao1/project/ontology/data/reactome/hierarchy/";	# folder name of input sorted reactome hierarchy files 
perf_file		<- "/home/yunhao1/project/DTox/data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";	# name of optimal model performance file 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome
output_folder		<- "data/compound_target_probability_tox21_null/";	# folder name of output files 
N_cores			<- 300;	# number of CPUs
job_name		<- "dtox_compound_target_probability_tox21_null";	# job name

## 1. Process input and output files 
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
# create a folder for each dataset to store returned results of DTox model  
data_output_folder <- sapply(perf_df$dataset_name, function(pddn){
	# name the folder after dataset  
	pddn_s <- strsplit(pddn, "tox21-")[[1]][[2]];
        pddn_folder <- paste(output_folder, pddn_s, "/", sep = "");
        system(paste("mkdir", pddn_folder, sep = " "));
        return(pddn_folder);
});
# list all input outcome-shuffled training-testing data files  
all_train_files <- list.files(train_data_folder);

## 2. Generate commands for jobs 
commands <- mapply(function(pddn, dof, ph1, ph2, ph3, ph4, ph5){
	# training file
	atf_id <- sapply(all_train_files, function(atf) length(strsplit(atf, pddn)[[1]]));
	pddn_train <- all_train_files[atf_id == 2];
	# testing file
	pddn_test <- paste(test_data_folder, pddn, "whole_data.tsv_test.tsv", sep = "_");
	# root pathway file 
	ph1_root <- paste(hierarchy_folder, ph1, "_st_0_root.tsv", sep = "");
	# parent/children node connection file 
	ph1_relation <- paste(hierarchy_folder, ph1, "_st_0_knowledge_by_node.tsv", sep = "");
	# node gene number file  
	ph1_size <- paste(hierarchy_folder, ph1, "_st_0_node_size.tsv", sep = "");
	# node layer number file  
	ph1_layer <- paste(hierarchy_folder, ph1, "_st_0_layer.tsv", sep = "");
	# put together command that includes files above and learned optimal hyperparameter setting 
	pddn_commands <- sapply(pddn_train, function(pt){
		pt_train <- paste(train_data_folder, pt, sep = "");
		pt_output <- paste(dof, pt, sep = "");
		pt_command <- paste("python", "src/dtox.py", pt_train, pddn_test, outcome_col, ph1_root, ph1_relation, ph1_size, ph1_layer, ph2, pt_output, ph3, ph4, ph5, sep = " ");
		return(pt_command);
	});
	return(pddn_commands);
}, perf_df$dataset_name, data_output_folder, perf_hp[,1], perf_hp[,2], perf_hp[,3], perf_hp[,4], perf_hp[,5], SIMPLIFY = F);
commands <- unlist(commands);
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
