# !/usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that run simple machine learning models on Tox21 datasets under different hyperparameter settings 


## functions
source("src/functions.R");


## 0. Input arguments
Args			<- commandArgs(T);
input_data_folder	<- Args[1];		# folder name of input training/validation data files 
output_folder		<- Args[2];		# folder name of output files 	
method			<- Args[3];		# name of classification method to be used: 'RandomForest', 'XGBoost'
N_cores			<- Args[4];		# number of CPUs 
job_name		<- Args[5];		# job name 
outcome_col             <- "assay_outcome";	# name of column that contains Tox21 assay outcome  

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
# create a folder for each dataset to store returned results of learning model 
data_output_folder <- sapply(data_file_info, function(dfi){
	# name the folder after dataset  
	dfi_folder <- paste(output_folder, dfi, "/", sep = "");
	system(paste("mkdir", dfi_folder, sep = " "));
	return(dfi_folder);
});

## 2. Specify hyperparameter values for model training  
hp_list <- NULL 
hp_list[[1]] <- c("n_estimators:50", "n_estimators:100") 
# specify hyperparameter values for random forest model  
if(method == "RandomForest"){
	hp_list[[2]] <- c("criterion:gini", "criterion:entropy")
	hp_list[[3]] <- sapply(1:10*0.1, function(x) paste("max_features", x, sep = ":"))
	hp_list[[4]] <- sapply(1:7*3-1, function(x) paste("min_samples_split", x, sep = ":"))
	hp_list[[5]] <- sapply(1:5*5-4, function(x) paste("min_samples_leaf", x, sep = ":"))
	hp_list[[6]] <- c("bootstrap:True", "bootstrap:False")
}
# specify hyperparameter values for gradient boosting model
if(method == "XGBoost"){
	hp_list[[2]] <- sapply(1:10, function(x) paste("max_depth", x, sep = ":"))
	hp_list[[3]] <- c("learning_rate:1e-3", "learning_rate:1e-2", "learning_rate:1e-1")
	hp_list[[4]] <- sapply(1:10*0.1, function(x) paste("subsample", x, sep = ":"))
	hp_list[[5]] <- sapply(1:5*5-4, function(x) paste("min_child_weight", x, sep = ":"))
}
# generate all possible combinations of different hyperparameter values  
hp_combo <- generate.all.possible.hyperparameter.combinations(hp_list, ",");

## 3. Generate parts of commands  
part <- NULL
# put together first part of simple learning command  
part[[1]] <- mapply(function(adf, dof){
        # training dataset file
        adf_train <- paste(input_data_folder, adf, "_train.tsv", sep = "");
        # validation dataset file 
        adf_test <- paste(input_data_folder, adf, "_test.tsv", sep = "");
        # output file 
        adf_output <- paste(dof, adf, sep = "");
        # put  first part that includes the files above 
        adf_command <- paste("python", "src/simple/simple.py", adf_train, adf_test, outcome_col, adf_output, method, sep = " ");
        return(adf_command);
}, all_data_files, data_output_folder);
# put together second part of simple learning command that includes model hyperparameter setting  
part[[2]] <- hp_combo;

## 4. Generate commands for jobs 
commands <- generate.all.possible.hyperparameter.combinations(part);
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
