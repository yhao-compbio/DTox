# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2022
## This script analyzes the relative efficiency/performance of DTox under alternative settings of early stopping criterion 


## Functions
# This function generate sorted bootstrap samples of input vector
bootstrap <- function(vec, sampleTime){
	len <- length(vec);
	bootList <- lapply(1:sampleTime,function(i) sample(vec,len,replace=T));
	mean_boot <- sapply(bootList,mean);
	return(sort(mean_boot));
}
# This function generate the 95% confidence interval of input vector by bootstrapping 
CI.boot <- function(vec, sampleTime){
	boot <- bootstrap(vec,sampleTime);
	i1 <- round(0.025*sampleTime);
	i2 <- round(0.975*sampleTime);
	CI <- c(boot[[i1]],boot[[i2]]);
	return(CI);
}


## 0. Input arguments
loss_folder	<- "data/compound_target_probability_tox21_implementation/";
perf_file	<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";   # name of optimal model performance file
tolerance	<- 1:20;
max_epoch	<- 200;
output_folder	<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_early_stop_summary";

## 1. Compute the relative efficiency/performance of DTox under alternative tolerance settings (parameter of early stopping criterion)  
# read in optimal performance info of Tox21 datasets as data frame  
perf_df <- read.delim(file = perf_file, header = T, sep = "\t");
# iterate by Tox21 dataset of interest, compute the relative efficiency/performance of DTox under alternative tolerance settings
dataset_stop_summary <- mapply(function(pddn, pdhs){
	# obtain the loss file name of current Tox21 dataset, read in the training/testing loss tracker as data frame 
	pddn_s <- strsplit(pddn, "tox21-")[[1]][[2]];
	pddn_loss_file <- paste(loss_folder, pddn_s, "/compound_target_fingerprint_maccs_probability_", pddn, "_whole_data.tsv_", pdhs, "_loss.tsv", sep = "");
	pddn_loss_df <- read.delim(file = pddn_loss_file, sep = "\t", header = T);
	# compute the optimal performance point at each epoch  
	pddn_test_loss <- pddn_loss_df$testing_total_loss;
	pddn_optimal <- pddn_test_loss[[1]];
	pddn_optimal_id <- 1;
	pddn_optimal_loc <- rep(1, length(pddn_test_loss));
	for(lptl in 1:length(pddn_test_loss)){
		if(pddn_test_loss[[lptl]] < pddn_optimal){
			pddn_optimal <- pddn_test_loss[[lptl]]
			pddn_optimal_id <- lptl
		}
		pddn_optimal_loc[[lptl]] <- pddn_optimal_id
	}	
	# compute the relative efficiency (running epochs of alternative setting vs default setting) of DTox under alternative tolerance settings 
	pddn_loc_tol <- 1:length(pddn_test_loss) - pddn_optimal_loc;
	tol_run_epochs <- sapply(tolerance, function(tol){
		tol_id <- which(pddn_loc_tol %in% tol)
		if(length(tol_id) == 0)	tol_run <- max_epoch
		else	tol_run <- min(tol_id)
		return(tol_run)
	});
	tol_run_time_ratio <- tol_run_epochs/length(pddn_test_loss);
	# compute the relative performance (testing loss improvement of alternative setting vs default setting) of DTox under alternative tolerance settings
	tol_optimal_epochs <- mapply(function(tre, tol){
		if(tre == max_epoch)	tre_opt <- pddn_optimal_id
		else	tre_opt <- tre - tol 
		return(tre_opt)
	}, tol_run_epochs, tolerance);
	tol_performance_ratio <- sapply(tol_optimal_epochs, function(toe){
		toe_perf <- (pddn_test_loss[[1]] - pddn_test_loss[[toe]])/(pddn_test_loss[[1]] - pddn_test_loss[[pddn_optimal_id]])
	}); 		
	tol_summary_df <- data.frame(pddn, tolerance, tol_optimal_epochs, tol_run_time_ratio, tol_performance_ratio);
	colnames(tol_summary_df) <- c("dataset_name", "patience", "optimal_epoch", "run_time_ratio", "performance_ratio");	
	return(tol_summary_df);
}, perf_df$dataset_name, perf_df$hyperparameter_setting, SIMPLIFY = F);
# aggregate the computed relative efficiency/performance from all datasets, write to output file  
all_dataset_stop_summary <- do.call(rbind, dataset_stop_summary);
all_dataset_stop_summary_file <- paste(output_folder, "_by_dataset.tsv", sep = "");
write.table(all_dataset_stop_summary, file = all_dataset_stop_summary_file, sep = "\t", col.names = T, row.names = F, quote = F);

## 2. Summerize the computed the relative efficiency/performance of DTox over all datasets 
# iterate by alternative tolerance settings, summerize
tol_stop_summary <- mapply(function(tol){
	# compute the average relative efficiency across datasets and its 95% confidence interval  
	tol_id <- which(all_dataset_stop_summary$patience %in% tol);
	tol_run_time_ratio_mean <- mean(all_dataset_stop_summary$run_time_ratio[tol_id]);
	tol_run_time_ratio_ci <- CI.boot(all_dataset_stop_summary$run_time_ratio[tol_id], 1000);
	# compute the average relative performance across datasets and its 95% confidence interval   
	tol_performance_ratio_mean <- mean(all_dataset_stop_summary$performance_ratio[tol_id]);
	tol_performance_ratio_ci <- CI.boot(all_dataset_stop_summary$performance_ratio[tol_id], 1000);
	tol_ratio_vec <- c(tol_run_time_ratio_mean, tol_run_time_ratio_ci, tol_performance_ratio_mean, tol_performance_ratio_ci);
	return(tol_ratio_vec);
}, tolerance);
# aggregate the computed statistics of relative efficiency/performance from all alternative tolerance settings, write to output file  
tol_stop_summary_df <- data.frame(tolerance, t(tol_stop_summary));
colnames(tol_stop_summary_df) <- c("patience", "run_time_ratio", "run_time_ratio_95ci_lower", "run_time_ratio_95ci_upper", "performance_ratio", "performance_ratio_95ci_lower", "performance_ratio_95ci_upper");
tol_stop_summary_file <- paste(output_folder, "_by_patience.tsv", sep = "");
write.table(tol_stop_summary_df, file = tol_stop_summary_file, sep = "\t", col.names = T, row.names = F, quote = F);

