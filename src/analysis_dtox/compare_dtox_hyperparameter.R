# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script normalizes DTox model performance across a query hyperparameter for comparison 


## functions
source("src/functions.R");
library(RColorBrewer);
library(gplots);


## 0. Input arguments 
dtox_optimal_file	<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";	# name of optimal DTox model performance (on Tox21 datasets) file  
dtox_result_file	<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_performance_summary.tsv";	# name of DTox model performance (including all hyperparameter settings) file 
tox21_annot_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_assay_info.tsv"	# name of Tox21 assay annotation file 
query_hp		<- "rt";	# code of hyperparameter to be compared ('rt' represents root pathway of DTox model) 
out_file		<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_rt_training_root_loss_normalized_comparison_by_dataset.tsv";	# name of output comparison result file  

## 1. Extract optimal setting values of query hyperparameter as well as other hyperparameters  
# read in optimal DTox model performance as data frame 
dtox_optimal_df <- read.delim(file = dtox_optimal_file, header = T, sep = "\t");
# obtain index of query hyperparameter in hyperparameter setting string (separated by underscore)   
sample_hp <- dtox_optimal_df$hyperparameter_setting[[1]]; 
sh_s <- strsplit(sample_hp, "_")[[1]];
qh_id <- which(sh_s %in% query_hp);
# obtain indexes of all other hyperparameters in hyperparameter setting string (separated by underscore)   
others_id <- setdiff(1:length(sh_s), c(qh_id, qh_id+1));
# use the obtained indexes to extract optimal setting values of query hyperparameter as well as other hyperparameters  
dtox_hps <- mapply(function(dodhs){
	# extract values of query hyperparameter
	dodhs_s <- strsplit(dodhs, "_")[[1]]; 
	dodhs_optimal <- paste(dodhs_s[c(qh_id, qh_id+1)], collapse = "_");
	# extract values of other hyperparameters 
	dodhs_others <- paste(dodhs_s[others_id], collapse = "_");
	return(c(dodhs_optimal, dodhs_others));
}, dtox_optimal_df$hyperparameter_setting);

## 2. Normalize DTox model performance across query hyperparameter setting for each Tox21 dataset
# read in DTox model performance of all hyperparameter settings as data frame 
dtox_result_df <- read.delim(file = dtox_result_file, header = T, sep = "\t");
# Normalize DTox model performance across all query hyperparameter settings, iterate by Tox21 dataset  
dtox_hp_norm_map <- mapply(function(doddn, dh2){
	# obtain DTox model performance of all hyperparameter settings for the current Tox21 dataset  
	doddn_id <- which(dtox_result_df$dataset_name %in% doddn);
	doddn_hps <- dtox_result_df$hyperparameter_setting[doddn_id];
	# use the obtained indexes (in 1) to extract all setting values of query hyperparameter as well as other hyperparameters for the current Tox21 dataset 
	dh_hps <- mapply(function(dh){
		dh_s <- strsplit(dh, "_")[[1]];
		ds_hp <- paste(dh_s[c(qh_id, qh_id+1)], collapse = "_"); 	
		ds_others <- paste(dh_s[others_id], collapse = "_");
		return(c(ds_hp, ds_others));
	}, doddn_hps);
	# select the settings with the same value of all other hyperparameters as optimal setting, except for the query hyperparameter 
	dh2_id1 <- which(dh_hps[2, ] %in% dh2);
	dh2_id <- doddn_id[dh2_id1];
	# obtain the training root loss of selected settings on current Tox21 dataset  
	dh2_scores <- dtox_result_df$training_root_loss[dh2_id];
	# normalize the obtained training root loss by mean and sd (Z transformation)
	dh2_scores_nml <- -(dh2_scores - mean(dh2_scores))/sd(dh2_scores); 
	names(dh2_scores_nml) <- dh_hps[1, dh2_id1];
	return(dh2_scores_nml);
}, dtox_optimal_df$dataset_name, dtox_hps[2, ]);
dtox_hp_norm_map <- t(dtox_hp_norm_map);

## 3. Annotate computed results with Tox21 dataset info, write to output file  
# read in Tox21 assay annotation info as data frame  
tox21_annot_df <- read.delim(file = tox21_annot_file, header = T, sep = "\t");
assay_target <- tox21_annot_df$assay_target;
names(assay_target) <- tox21_annot_df$protocol_name;
# add cell line name after assay name of each Tox21 dataset, make it the row index of normalized performance data frame  
assay_target[["tox21-ar-mda-kb2-luc-antagonist-p2"]] <- "AR-MDA antagonist";
assay_target <- mapply(function(at, tadcl) paste(at, " (", tadcl, ")", sep = ""), assay_target, tox21_annot_df$cell_line);
rownames(dtox_hp_norm_map) <- assay_target[rownames(dtox_hp_norm_map)];
# write normalized performance data frame to output file  
write.table(dtox_hp_norm_map, file = out_file, sep = "\t", col.names = T, row.names = T, quote = F);

