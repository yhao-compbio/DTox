# !/usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that runs DTox interpretation procedure on optimal models trained for Tox21 datasets


## functions
source("src/functions.R");


## 0. Input arguments
query_data_folder	<- "/home/yunhao1/project/TTox/data/compound_target_tox21_data/fingerprint_maccs_probability/compound_target_fingerprint_maccs_probability";	# folder name of input whole Tox21 datasets  
hierarchy_folder	<- "/home/yunhao1/project/ontology/data/reactome/hierarchy/";	# folder name of input sorted reactome hierarchy files 
model_folder		<- "data/compound_target_probability_tox21_implementation/"	# folder name of trained DTox models using original Tox21 datasets  
null_folder		<- "data/compound_target_probability_tox21_null/"	# folder name of trained DTox models using outcome-shuffled Tox21 datasets   
perf_file		<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";	# name of optimal model performance file 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome
output_folder		<- "data/compound_target_probability_tox21_interpret/";	# folder name of output files 
N_cores			<- 3;	# number of CPUs
job_name		<- "interpret_dtox_compound_target_probability_tox21_implementation";	# job name

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

## 2. Process input trained DTox model files 
# list all files in folder that contains trained DTox models using outcome-shuffled Tox21 datasets 
all_null_files <- list.files(null_folder, recursive = TRUE);
# select trained DTox model files (ending with 'model.pt') 
model_id <- sapply(all_null_files, function(anf) length(strsplit(anf, "model")[[1]]));
all_null_files <- all_null_files[model_id == 2];
# Match trained DTox model files with Tox21 datasets 
null_name_files <- sapply(perf_df$dataset_name, function(pddn){
	# obtain name of Tox21 dataset  
	pddn_s <- strsplit(pddn, "tox21-")[[1]][[2]];
	# select files that contain DTox model trained for the current dataset    
	anf_id <- sapply(all_null_files, function(anf) length(strsplit(anf, pddn_s)[[1]]));
	pddn_files <- all_null_files[anf_id == 3];
	pddn_files <- sapply(pddn_files, function(pf) paste(null_folder, pf, sep = ""));
	# write names of selected files to output file  
	pddn_out <- paste(null_folder, pddn_s, "_null_files.txt", sep = ""); 
	writeLines(pddn_files, pddn_out);
	return(pddn_out);
});

## 3. Generate parts of commands 
part <- NULL 
# generate first part of command that includes input dataset/hierarchy/model file info 
part[[1]] <- mapply(function(ph1, ph2, ph3, ph4, ph5, pddn, pdhs, nnf){
	# root pathway file  
	ph1_root <- paste(hierarchy_folder, ph1, "_st_0_root.tsv", sep = "");
	# parent/children node connection file
	ph1_relation <- paste(hierarchy_folder, ph1, "_st_0_knowledge_by_node.tsv", sep = "");
	# node gene number file 
	ph1_size <- paste(hierarchy_folder, ph1, "_st_0_node_size.tsv", sep = "");
	# node layer number file  
	ph1_layer <- paste(hierarchy_folder, ph1, "_st_0_layer.tsv", sep = "");
	# node name file 
	ph1_node <- paste(hierarchy_folder, ph1, "_st_0_node.tsv", sep = "");
	# trained DTox model file 
	pddn_s <- strsplit(pddn, "tox21-")[[1]][[2]];
	pddn_model <- paste(model_folder, pddn_s, "/compound_target_fingerprint_maccs_probability_tox21-", pddn_s, "_whole_data.tsv_", pdhs, "_model.pt", sep = "");
	# whole Tox21 dataset file
	pddn_data <- paste(query_data_folder, pddn, "whole_data.tsv", sep = "_"); 
	# put together first part command that includes the files above 
	pddn_command <- paste("python", "src/interpret_dtox.py", ph1_root, ph1_relation, ph1_size, ph1_layer, ph2, ph3, ph1_node, pddn_model, nnf, pddn_data, outcome_col, sep = " "); 
	return(pddn_command);
}, perf_hp[,1], perf_hp[,2], perf_hp[,3], perf_hp[,4], perf_hp[,5], perf_df$dataset_name, perf_df$hyperparameter_setting, null_name_files);
# generate second part of command that includes LRP rule and hyperparameter values     
part[[2]] <- c("alpha-beta 1 0", "alpha-beta 2 1", "alpha-beta 3 2", "gamma-epsilon 0.1 0.1", "gamma-epsilon 0.1 0.01", "gamma-epsilon 0.1 0.001", "gamma-epsilon 0.01 0.1", "gamma-epsilon 0.01 0.01", "gamma-epsilon 0.01 0.001", "gamma-epsilon 0.001 0.1", "gamma-epsilon 0.001 0.01", "gamma-epsilon 0.001 0.001");
# create a folder for each rule to store retured interpretation results  
part2_folders <- sapply(part[[2]], function(p2){
	# name the folder after LRP rule paired with hyperparameters  
	p2_s <- strsplit(p2, " ")[[1]]
	p2_folder <- paste(output_folder, paste(p2_s, collapse = "_"), "/", sep = "");
	system(paste("mkdir", p2_folder, sep = " "));
	return(p2_folder);
});

## 4. Generate commands for jobs 
commands1 <- generate.all.possible.hyperparameter.combinations(part);
# add output file name to each command 
commands <- sapply(commands1, function(com){
	# obtain trained DTox model file name of current command  
	com_s <- strsplit(com, " ")[[1]];
	N_cs <- length(com_s);
	com_model_s <- strsplit(com_s[[10]], "/", fixed = T)[[1]];
	com_model <- com_model_s[[length(com_model_s)]];
	# obtain LRP rule and its hyperparameters of current command
	com_rule <- paste(com_s[(N_cs-2):N_cs], collapse = "_");
	# name output file after the trained DTox model file paired with LRP rule 
	com_op_file <- paste(output_folder, com_rule, "/", com_model, "_", com_rule, sep = "");
	# add outut file name to end of current command 
	com_new <- paste(com, com_op_file, sep = " ");
	return(com_new);
});
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
