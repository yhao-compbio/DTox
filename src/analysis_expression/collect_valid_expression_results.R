# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script collects computed differential expression proportions of compounds from interpretation-validation result files, and performs t test to compare proportion among significant DTox paths vs among background DTox paths


## 0. Input arguments 
tox21_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_10k_library_info.tsv";	# name of input Tox21 compound library info file 
valid_folder	<- "data/compound_target_probability_tox21_interpret_expression/validation_result/";	# folder name of DTox interpretation-expression validation files   
min_sample	<- 10;	# minimal number of expression-validated compounds in a single assay for the result to be collected   
output_file	<- "data/compound_target_probability_tox21_interpret_expression/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_expression_validation";	# folder name of output collection result files 

## 1. Process Tox21 compound library info data 
# read in Tox21 compound library info as data frame  
tox21_df <- read.delim(file = tox21_file, header = T, sep = "\t");
# add "CID_" before each Pubchem CID of compound, convert compound name characters to lower case  
tox21_df$cid <- sapply(tox21_df$PUBCHEM_CID, function(tdpc) paste("CID", tdpc, sep = "_"));
tox21_df$compound_name <- tolower(tox21_df$SAMPLE_NAME);
# create unique mapping data frame between compound CIDs and names
tox21_df <- tox21_df[, c("cid", "compound_name")];
td_cid_uni <- unique(tox21_df$cid);
tcu_ids <- sapply(td_cid_uni, function(tcu) which(tox21_df$cid %in% tcu)[[1]]);
tox21_df <- tox21_df[tcu_ids, ];

## 2. Process DTox interpretation-expression validation files  
# list all files in DTox interpretation-expression validation folder 
valid_files <- list.files(valid_folder); 
# select files that contain proportion of differential expression comparison results (between significant DTox paths and background DTox paths)
vf_id <- sapply(valid_files, function(vf) length(strsplit(vf, "result")[[1]]));
valid_files <- valid_files[vf_id == 2];
# obtain result file info from file name
valid_file_info	<- mapply(function(vf){
	# LRP propogation rule and hyperparameter values (gamma and epsilon)
	vf_s <- strsplit(vf, "_")[[1]]; 
	rule <- paste(vf_s[1:3], collapse = "_");
	# compound dose of pertubation for generating LINCS gene expression data 
	dose <- vf_s[[7]];
	# time after pertubation for generating LINCS gene expression data 
	time <- vf_s[[8]];
	# FDR threshold for differential expression analysis   
	fdr_cut <- vf_s[[10]];
	return(c(rule, dose, time, fdr_cut));
}, valid_files);
valid_file_info <- t(valid_file_info);
colnames(valid_file_info) <- c("rule", "dose", "time", "fdr_threshold");
valid_files <- sapply(valid_files, function(vf) paste(valid_folder, vf, sep = ""));
 
## 3. Collect computed differential expression proportions of compounds from result files, compare proportion among significant DTox paths vs among background DTox paths by t test 
# iterate by differential expression comparison result file, collect results from each file 
valid_collect <- mapply(function(vf, nvfi){
	# read in differential expression comparison data frame of current result file  
	vf_df <- read.delim(file = vf, sep = "\t", header = T);
	# select compounds with at least one differentially expressed DTox path (background proportion > 0) 
	g0_id <- which(vf_df$background_sig_ratio > 0);
	vf_df <- vf_df[g0_id, ];
	# select Tox21 assay datasets with number of compounds that pass minimal threshold  
	assay_count <- table(vf_df$assay);
	ac_id <- which(assay_count >= min_sample);
	# if no dataset satisfies minimal number requirement, stop 
	if(length(ac_id) == 0)	return(ls = list(result = logical(0), pv = logical(0)))
	else{
		# obtain computed differential expression proportions of all compounds for selected assay datasets   
		vf_assays <- names(assay_count)[ac_id]
		va_id <- which(vf_df$assay %in% vf_assays)
		vf_df <- vf_df[va_id, ]
		result_df <- data.frame(t(valid_file_info[nvfi,]), vf_df)
		# iterate by dataset of Tox21 assay, compare proportion of differential expression between significant DTox paths and background DTox paths 
		compare_pv <- mapply(function(va){
			# obtain proportion among significant DTox paths, and among background DTox paths for current dataset 
			va_id <- which(vf_df$assay %in% va)
			va_ip_ratio <- vf_df$interpret_sig_ratio[va_id]
			va_bg_ratio <- vf_df$background_sig_ratio[va_id]
			# perform one-sided wilcox test to compare proportion among significant DTox paths vs among background DTox paths, obtain P-value of t test    
			va_w_pv <- wilcox.test(va_ip_ratio, va_bg_ratio, paired = T, alternative = "greater")$p.value
			va_w_pv <- round(va_w_pv, 2)
			return(c(length(va_id), va_w_pv))
		}, vf_assays)
		compare_pv <- t(compare_pv)
		pv_df <- data.frame(t(valid_file_info[nvfi,]), vf_assays, compare_pv)
		colnames(pv_df) <- c(colnames(valid_file_info), "assay", "n_compounds", "wilcox_pv")
		return(ls = list(result = result_df, pv = pv_df))	
	}
}, valid_files, 1:nrow(valid_file_info), SIMPLIFY = FALSE);
# extract data frames that contain computed differential expression proportions of compounds, combine data frames from multple assay datasets together  
valid_result_list <- lapply(valid_collect, function(vc) vc[[1]]);
vrl_len <- sapply(valid_result_list, length);
valid_result_list <- valid_result_list[vrl_len > 0];
valid_result_df <- do.call(rbind, valid_result_list);
# merge extracted data frame with Tox21 compound CID-name mapping data frame to show both CIDs and names of compounds  
valid_result_df1 <- merge(valid_result_df, tox21_df, by.x = "compound", by.y = "cid");
valid_result_df1 <- valid_result_df1[, c("rule", "fdr_threshold", "dose", "time", "assay", "compound", "compound_name", "interpret_sig_ratio", "background_sig_ratio")];
# sort differential expression proportion data frame by specified parameters, write sorted data frame to output file 
vrd_od <- order(valid_result_df1$rule, valid_result_df1$fdr_threshold, valid_result_df1$dose, valid_result_df1$time, valid_result_df1$assay);
valid_result_df1 <- valid_result_df1[vrd_od, ];
write.table(valid_result_df1, file = paste(output_file, "_all_result.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# extract data frames that contain computed t test P-value comparing proportion among significant DTox paths vs among background DTox paths, combine data frames from multple assay datasets together   
valid_pv_list <- lapply(valid_collect, function(vc) vc[[2]]);
vpl_len <- sapply(valid_pv_list, length);
valid_pv_list <- valid_pv_list[vpl_len > 0];
valid_pv_df <- do.call(rbind, valid_pv_list);
# sort P-value data frame by specified parameters, write sorted data frame to output file 
vpd_od <- order(valid_pv_df$rule, valid_pv_df$fdr_threshold, valid_pv_df$dose, valid_pv_df$time, valid_pv_df$assay);
valid_pv_df <- valid_pv_df[vpd_od, ];
write.table(valid_pv_df, file = paste(output_file, "_result_compare_summary.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
