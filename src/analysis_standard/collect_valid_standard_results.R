# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script collects observed outcome and expected probability of compounds from interpretation-validation result files, then compute observed and expected proportion of validated compounds based on collected results


## 0. Input arguments 
valid_folder	<- "data/compound_target_probability_tox21_interpret_standard/validation_result/";	# folder name of DTox interpretation-standard validation files 
output_file	<- "data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_standard_validation"; # folder name of output collection result files
N_sample	<- 1000;	# number of sampling times for computing expected proportion of validated compounds 

## 1. Process DTox interpretation-standard validation files
# list all files in DTox interpretation-standard validation folder 
valid_files <- list.files(valid_folder); 
# extract LRP propagation rule parameters from interpretation result files
valid_rules <- sapply(valid_files, function(vf) strsplit(vf, "_standard")[[1]][[1]]);
valid_files <- sapply(valid_files, function(vf) paste(valid_folder, vf, sep = ""));
 
## 2. Collect observed outcome and expected probability from result files, then analyze collected results
# set seed number for sampling  
set.seed(0);
# iterate by DTox interpretation-standard validation file, collect results and use it to compute observed and expected proportion of validated compounds (identified with at least one standard path among the significant DTox paths) 
valid_collect <- mapply(function(vf, vr){
	# read in standard validation result as data frame, which contains computed expected probablity and observed outcome of compounds 
	vf_df <- read.delim(file = vf, sep = "\t", header = T);
	result_df <- cbind(vr, vf_df);
	colnames(result_df) <- c("rule", colnames(vf_df));
	# iterate by unique Tox21 assay from the data frame, compute observed and expected proportion of validated compounds
	vf_assays <- unique(vf_df$assay);
	vf_sample_compare <- lapply(vf_assays, function(va){
		# select rows from validation result data frame that are relevant to the current Tox21 assay 
		va_id <- which(vf_df$assay %in% va);
		va_df <- vf_df[va_id, ];
		# compute observed proportion of validated compounds 
		N_va <- length(va_id);
		va_observed <- sum(va_df$whether_contain_standard)/N_va;	
		# compute expected proportion using random sampling 
		va_bg_samples <- sapply(1:N_sample, function(ns){
			# iterate by compound, sample with the computed expected probability and obtain sampled outcome for each compound  
			ns_sample <- sapply(va_df$bg_standard_ratio, function(vdbsr){
				indi_vec <- c(1, 0);
				indi_prob <- c(vdbsr, 1-vdbsr);
				indi_s <- sample(indi_vec, 1, prob = indi_prob);
				return(indi_s);
			});
			# compute expected proportion of validated compounds
			ns_ratio <- sum(ns_sample)/N_va;
			return(ns_ratio);
		});
		# sort sampled expected proportion of validated compounds, obtain median and 95% confidence interval 	
		va_bg_samples <- sort(va_bg_samples);
		va_sample_df <- data.frame(va, va_bg_samples);
		colnames(va_sample_df) <- c("assay", "bg_sample");
		va_bg_lower <- va_bg_samples[[as.integer(0.025*N_sample)]];
		va_bg_med <- va_bg_samples[[as.integer(0.5*N_sample)]];
		va_bg_upper <- va_bg_samples[[as.integer(0.975*N_sample)]];
		# store computed observed/expected proportion stats in a data frame
		va_compare_df <- data.frame(va, va_observed, va_bg_med, va_bg_lower, va_bg_upper);
		colnames(va_compare_df) <- c("assay", "standard_observed_ratio", "standard_bg_med", "standard_bg_95ci_lower", "standard_bg_95ci_upper"); 
		return(ls = list(sample = va_sample_df, compare = va_compare_df));
	});	
	# merge sampled expected proportions from multiple assays into one data frame 
	vf_sample_list <- lapply(vf_sample_compare, function(vsc) vsc[[1]]);
	vf_sample_df <- do.call(rbind, vf_sample_list);
	sample_df <- cbind(vr, vf_sample_df);
	colnames(sample_df) <- c("rule", colnames(vf_sample_df));   	
	# merge computed observed/expected proportion stats from multiple assays into one data frame  
	vf_compare_list <- lapply(vf_sample_compare, function(vsc) vsc[[2]]);
	vf_compare_df <- do.call(rbind, vf_compare_list);
	compare_df <- cbind(vr, vf_compare_df);
	colnames(compare_df) <- c("rule", colnames(vf_compare_df));
	return(ls = list(result = result_df, sample = sample_df, compare = compare_df))	
}, valid_files, valid_rules, SIMPLIFY = FALSE);
# merge collected expected probablity and observed outcome of compounds from multiple LRP rules into one data frame, write to output file   
valid_result_list <- lapply(valid_collect, function(vc) vc[[1]]);
valid_result_df <- do.call(rbind, valid_result_list);
write.table(valid_result_df, file = paste(output_file, "_all_result.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# merge sampled expected proportions from multiple LRP rules into one data frame, write to output file   
valid_sample_list <- lapply(valid_collect, function(vc) vc[[2]]);
valid_sample_df <- do.call(rbind, valid_sample_list);
write.table(valid_sample_df, file = paste(output_file, "_bg_sample.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# merge computed observed/expected proportion stats from multiple LRP rules into one data frame, write to output file 
valid_compare_list <- lapply(valid_collect, function(vc) vc[[3]]);
valid_compare_df <- do.call(rbind, valid_compare_list);
write.table(valid_compare_df, file = paste(output_file, "_result_compare_summary.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
