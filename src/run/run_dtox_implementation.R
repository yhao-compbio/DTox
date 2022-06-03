# !/usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that run DTox model on Tox21 datasets under Reactome pathway hierarchy


## functions
source("src/functions.R");


## 0. Input arguments
Args			<- commandArgs(T);
input_data_folder	<- Args[1];		# folder name of input training-testing/validation data files  
structure_feature	<- Args[2];		# indicator that represents whether input featues are structural (1: Yes, 0: No) 
hierarchy_folder	<- Args[3];		# folder name of input sorted reactome hierarchy files  
output_folder		<- Args[4];		# folder name of output files 
tune_alpha		<- Args[5];		# whether to tune alpha hyperparameter   
N_cores			<- Args[6];		# number of CPUs 
job_name		<- Args[7];		# job name 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome 

## 1. Process Tox21 dataset files 
# list all files in dataset file folder
all_data_files <- list.files(input_data_folder);
# select files that contains whole datasets 
whole_id <- sapply(all_data_files, function(adf) length(strsplit(adf, "whole_data.tsv", fixed = T)[[1]]));
all_data_files <- all_data_files[whole_id == 1];
# exclude data summary files from file list 
sum_id <- sapply(all_data_files, function(adf) length(strsplit(adf, "summary", fixed = T)[[1]]));
all_data_files <- all_data_files[sum_id == 1];
# obtain name of Tox21 dataset from name of data file 
data_file_info <- sapply(all_data_files, function(adf){
	adf_s1 <- strsplit(adf, "tox21-")[[1]][[2]];
	adf_s2 <- strsplit(adf_s1, "_whole_data.tsv")[[1]][[1]];
	return(adf_s2);
});
# create a folder for each dataset to store returned results of DTox model 
data_output_folder <- sapply(data_file_info, function(dfi){
	# name the folder after dataset  
	dfi_folder <- paste(output_folder, dfi, "/", sep = "");
	system(paste("mkdir", dfi_folder, sep = " "));
	return(dfi_folder);
});

## 2. Process sorted Reactome hierarchy files 
# list all files in hierarchy file folder 
all_h_files <- list.files(hierarchy_folder);
# exclude files that use minimum pathay size of 10 as an input parameter for sorting hierarchy 
ps_id <- sapply(all_h_files, function(ahf) length(strsplit(ahf, "ps_10", fixed = T)[[1]]));
all_h_files <- all_h_files[ps_id == 1];
# exclude files that include chemical reaction terms in sorted hierarchy
re_id <- sapply(all_h_files, function(ahf) length(strsplit(ahf, "re_1", fixed = T)[[1]]));
all_h_files <- all_h_files[re_id == 1];
# obtain unique name prefix for each sorted Reactome hierarchy 
all_h_heads <- sapply(all_h_files, function(ahf) strsplit(ahf, "_st")[[1]][[1]]);
h_heads <- unique(all_h_heads);
# add structure feature indicator to name prefix of each sorted Reactome hierarchy  
h_heads_st <- sapply(h_heads, function(hh) paste(hh, "_st_", structure_feature, sep = ""));
# create all possible combinations of dataset file and name prefix of sorted hierarchies   
all_data_files_vec <- rep(all_data_files, each = length(h_heads));
data_output_folder_vec <- rep(data_output_folder, each = length(h_heads));
h_heads_vec <- rep(h_heads, times = length(all_data_files));
h_heads_st_vec <- rep(h_heads_st, times = length(all_data_files));

## 3. Generate parts of commands  
part <- NULL
# put together first part of DTox command
part[[1]] <- mapply(function(adfv, dofv, hhv, hhsv){
	# combined training and testing dataset file 
	adfv_train <- paste(input_data_folder, adfv, "_train.tsv", sep = "");
	# validation dataset file 
	adfv_test <- paste(input_data_folder, adfv, "_test.tsv", sep = "");
	# root pathway file 
	hhsv_root <- paste(hierarchy_folder, hhsv, "_root.tsv", sep = "");
	# parent/children node connection file 
	hhsv_relation <- paste(hierarchy_folder, hhsv, "_knowledge_by_node.tsv", sep = "");
	# node gene number file  
	hhsv_size <- paste(hierarchy_folder, hhsv, "_node_size.tsv", sep = "");
	# node layer number file  
	hhsv_layer <- paste(hierarchy_folder, hhsv, "_layer.tsv", sep = "");
	# minimal size of pathways 
	hhv_s <- strsplit(hhv, "_")[[1]];
	hs_id <- which(hhv_s %in% "ps") + 1;
	min_path <- hhv_s[[hs_id]];
	# output file
	adfv_output <- paste(dofv, adfv, "_", hhv, sep = "");
	# put together first part that includes the files/parameters above
	adfv_command <- paste("python", "src/dtox.py", adfv_train, adfv_test, outcome_col, hhsv_root, hhsv_relation, hhsv_size, hhsv_layer, min_path, adfv_output, sep = " ");
	return(adfv_command);
}, all_data_files_vec, data_output_folder_vec, h_heads_vec, h_heads_st_vec);
# put together second part of DTox command that includes maximal size of node modules, coefficient for auxiliary loss, coefficient for L2 regularization 
if(as.logical(tune_alpha))	part[[2]] <- c("20 0.1 0.0001", "20 0.5 0.0001", "20 1 0.0001")
if(!as.logical(tune_alpha))	part[[2]] <- c("20 0.5 0.0001")

## 4. Generate commands for jobs 
commands <- generate.all.possible.hyperparameter.combinations(part);
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
