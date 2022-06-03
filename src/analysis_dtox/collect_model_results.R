# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script collects machine learning model basic info and performance metrics from performance files


## Runctions
library(parallel);
source("src/functions.R");


## 0. Input arguments
Args		<- commandArgs(T);
input_folder	<- Args[1];     	# folder name of model result file
name_lower_id	<- as.integer(Args[2]);	# number that indicates lower bound index of dataset name (separated by underscore) in the file name
name_upper_id	<- as.integer(Args[3]);	# number that indicates upper bound index of dataset name (separated by underscore) in the file name
hyper_lower_id	<- as.integer(Args[4]);	# number that indicates lower bound index of hyperparameter (separated by underscore) in the file name 
hyper_upper_id	<- as.integer(Args[5]);	# number that indicates upper bound index of hyperparameter (separated by underscore) in the file name 
collect_mode	<- Args[6];		# character that indicates collection mode ('dtox' for collecting DTox model results, 'simple' for collecting simple machine learning model results, 'mlp' for collecting multi-layer perceptron model results) 
output_folder   <- Args[7];		# folder name of output collection files  

## 1. Obtain names of model performance files  
# list all files in model result file folder
all_files <- list.files(input_folder, recursive = T);
# select model performance files 
perf_id <- sapply(all_files, function(af) length(strsplit(af, "performance.", fixed = T)[[1]]));
all_perf_files <- all_files[perf_id == 2];
all_perf_files <- sapply(all_perf_files, function(apf) paste(input_folder, apf, sep = ""));
# extract dataset name and hyperparameter setting from name of each model performance file 
all_file_info <- mapply(function(apf){	
	apf_s <- strsplit(apf, "_")[[1]];
	if(is.na(name_upper_id))	apf_name <- apf_s[[name_lower_id]]
	if(!is.na(name_upper_id))	apf_name <- paste(apf_s[name_lower_id:name_upper_id], collapse = "_");
	apf_hyper <- paste(apf_s[hyper_lower_id: hyper_upper_id], collapse = "_");
	apf_info <- c(apf_name, apf_hyper);	
	names(apf_info) <- c("dataset_name", "hyperparameter_setting");
	return(apf_info);
}, all_perf_files);
all_file_info <- t(all_file_info);

## 2. Collect info from according to the specified mode  
# extract simple machine learning model basic info and performance metrics from performance files  
if(collect_mode == "simple"){
	file_perf_list <- mclapply(all_perf_files, function(apf) read.simple.performance.files(apf), mc.cores = 10)
}
# extract DTox model basic info and performance metrics from performance files   
if(collect_mode == "dtox"){
	file_perf_list <- mclapply(all_perf_files, function(apf) read.dtox.performance.files(apf), mc.cores = 10)
}
# extract MLP model basic info and performance metrics from performance files   
if(collect_mode == "mlp"){
	file_perf_list <- lapply(all_perf_files, function(apf) read.mlp.performance.files(apf))
}
# aggregate extracted info in data frame form  
file_perf_df <- do.call(rbind, file_perf_list);
combine_result_df <- cbind(all_file_info, file_perf_df);
# write collection results to output file  
combine_file <- paste(output_folder, "performance_summary.tsv", sep = "_"); 
write.table(combine_result_df, file = combine_file, sep = "\t", col.names = T, row.names = F, quote = F);

