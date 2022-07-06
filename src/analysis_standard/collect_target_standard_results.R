# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2022
## This script collects the validation results by DTox, LIME, and Read-across regarding the interpretation task of connecting active compounds to their respective target receptor in four nuclear receptor assays 


## 0. Input arguments  
dtox_file	<- "data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_standard_validation_all_result.tsv";
lime_folder	<- "data/compound_target_probability_tox21_interpret_standard/lime_interpret_result/compound_target_fingerprint_maccs_probability";
ra_folder	<- "data/compound_target_probability_tox21_interpret_standard/read_across_interpret_result/compound_structure_fingerprint_maccs";
output_file	<- "data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_target_standard_validation_result_compare_summary.tsv";

## 1. Read in validation results by DTox
dtox_df <- read.delim(file = dtox_file, sep = "\t", header = T);
dd_id <- which(dtox_df$rule %in% "gamma-epsilon_0.001_0.1");
dtox_df <- dtox_df[dd_id, ];
all_assays <- unique(dtox_df$assay);
 
## 2. Collect validation results by LIME
# iterate by assay of interest 
dtox_lime_compare <- mapply(function(aa){
	# obtain the name of current assay  
	aa_id <- which(dtox_df$assay %in% aa);
	aa_dtox_df <- dtox_df[aa_id,];
	# obtain the LIME result file for the current assay, read in validation results by LIME 
	aa_lime_file <- paste(lime_folder, aa, "interpret_by_lime_result.tsv", sep = "_");
	aa_lime_df <- read.delim(file = aa_lime_file, sep = "\t", header = T);
	# compute the number of validated compounds by DTox
	aa_dtox_lime_df <- merge(aa_dtox_df, aa_lime_df, by.x = "compound", by.y = "compound_cid");
	aa_dtox_count <- length(which(aa_dtox_lime_df$whether_contain_standard == 1))
	# compute the number of validated compounds by LIME under two feature relevance thresholds (loose/lax: > 0; strict: > average relevance of all features)
	adld_range_id <- grep(">", aa_dtox_lime_df$target_range);
	aa_lime_loose <- length(which(aa_dtox_lime_df$target_relevance[adld_range_id] > 0));
	aa_lime_strict <- length(which(aa_dtox_lime_df$target_relevance[adld_range_id] > aa_dtox_lime_df$average_relevance[adld_range_id]));
	# compute the proportion of validated compounds by DTox and LIME  
	aa_counts <- c(aa_lime_loose, aa_lime_strict, aa_dtox_count);
	aa_ratio <- aa_counts/nrow(aa_dtox_lime_df);
	return(aa_ratio);
}, all_assays);

## 3. Collect validation results by Read-across 
# iterate by assay of interest 
dtox_lime_ra_compare <- lapply(all_assays, function(aa){
	## obtain the Read-across result file for the current assay, read in validation results by Read-across  
	aa_ra_file <- paste(ra_folder, aa, "interpret_by_read_across_result.tsv", sep = "_"); 
	aa_ra_df <- read.delim(file = aa_ra_file, sep = "\t", header = T);
	## combine the computed proportions of validated compounds by DTox, LIME, and Read-across 
	aa_all_df <- cbind(aa_ra_df, t(dtox_lime_compare[, aa]));
	colnames(aa_all_df) <- c(colnames(aa_ra_df), "standard_observed_ratio_lime_loose", "standard_observed_ratio_lime_strict", "standard_observed_ratio_dtox");
	return(aa_all_df);
});

## 4. Output collected validation results  
dtox_lime_ra_compare_df <- do.call(rbind, dtox_lime_ra_compare);
write.table(dtox_lime_ra_compare_df, file = output_file, sep = "\t", row.names = F, col.names = T, quote = F);
