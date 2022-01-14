# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script compares the significant DTox paths detected under different hyperparameter settings of layer-wise relevance propagation rule on Tox21 datasets, and compute Jaccary Index to measure the similarity among distinct hyperparameter settings 

## functions 
source("src/functions.R");


## 0. Input arguments 
interpret_folder	<- "data/compound_target_probability_tox21_interpret/";	# folder name of input DTox model interpretation files  
out_sim_file		<- "data/compound_target_probability_tox21_interpret_analysis/compound_target_fingerprint_maccs_probability_gamma-epsilon_path_similarity.tsv";	# name of output similarity result file  

## 1. Process input DTox model interpretation files  
# list all files in DTox model interpretation result folder  
all_interpret_files <- list.files(interpret_folder, recursive = T);
# select interpretation results from generic rule, which contains 'gamma-epsilon' in file name 
all_rule_id <- sapply(all_interpret_files, function(aip) length(strsplit(aip, "gamma-epsilon")[[1]]));
all_rule_files <- all_interpret_files[all_rule_id > 1];
# select DTox VNN path relevance p-value files from generic rule, which contains 'path_relevance_pv' in file name  
all_path_id <- sapply(all_rule_files, function(aip) length(strsplit(aip, "path_relevance_pv")[[1]]));
all_path_files <- all_rule_files[all_path_id > 1];
# obtain the Tox21 dataset name of each p-value result file, group files by dataset name  
all_path_dataset <- sapply(all_path_files, function(apf) strsplit(apf, "_")[[1]][[8]]);
dataset_path_files <- group.vector.by.categories(all_path_dataset, all_path_files);

## 2. Compute Jaccard Index of significant DTox paths detected under different hyperparameter settings
# obtain significant DTox paths of each compound from path relevance file, iterate by Tox21 dataset 
dataset_path_list <- lapply(dataset_path_files, function(dpf){
	# iterate by generic rule hyperparameter setting (values of gamma and epsilon) of current dataset    
	dpf_list <- lapply(dpf, function(dp){
		# obtain DTox VNN path relevance p-value file of current dataset and setting
		dp_file <- paste(interpret_folder, dp, sep = "");
		# read in path relevance p-value info of compounds as data frame 
		dp_df <- read.delim(file = dp_file, header = T, sep = "\t");
		# group significant DTox paths by compound CID   	
		dp_path_list <- group.vector.by.categories(dp_df$cid, dp_df$path_id);
		return(dp_path_list);
	});
	return(dpf_list);
});
# compute Jaccard Index of significant DTox paths between hyperparameter setting pairs, iterate by Tox21 dataset  
dataset_path_sim_list <- mapply(function(dpf, dpl){
	# iterate by pair of generic rule hyperparameter setting (values of gamma and epsilon), through all possible pair combination  
	dpl_sim_list <- lapply(dpl, function(dp1){
		dp1_sim <- sapply(dpl, function(dp2){
			# find the common compounds between two settings of current pair  
			dp12 <- intersect(names(dp1), names(dp2));
			# iterate by compound, compute the Jaccard Index of significant paths between two settings of current pair  
			dp12_sim <- sapply(dp12, function(dd){
				d12i <- length(intersect(dp1[[dd]], dp2[[dd]]));
				d12u <- length(union(dp1[[dd]], dp2[[dd]]));
				d12_ji <- d12i/d12u;
				return(d12_ji);
			});
			# return median Jaccard Index across all common compounds 
			return(median(dp12_sim));
		});
		return(dp1_sim);
	});	
	# store the computed median Jaccard Index between hyperparameter settings in symmetric matrix form   
	dpl_sim_mat <- matrix(unlist(dpl_sim_list), nrow = length(dpf), byrow = T); 
	# name rows and columns of symmetric matrix after hyperparameter settings
	rownames(dpl_sim_mat) <- colnames(dpl_sim_mat) <- sapply(dpf, function(dp) strsplit(dp, "/", fixed = T)[[1]][[1]]);
	return(dpl_sim_mat);		
}, dataset_path_files, dataset_path_list, SIMPLIFY = F);
# compute average Jaccard Index across all Tox21 datasets  
dataset_path_sim <- Reduce('+', dataset_path_sim_list)/length(dataset_path_sim_list);
# write data frame to output file  
write.table(dataset_path_sim, file = out_sim_file, sep = "\t", row.names = T, col.names = T, quote = F);
