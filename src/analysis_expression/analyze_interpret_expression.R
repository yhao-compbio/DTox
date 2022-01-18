# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script analyzes gene expression-validated DTox paths from model interpretation on Tox21 datasets, and identifies recurrent differentially expressed DTox paths among compounds 


## 0. Input arguments 
sig_path_files <- c("data/compound_target_probability_tox21_interpret_expression/validation_result/gamma-epsilon_0.001_0.1_beta_trt_cp_10uM_6h_fdr_0.05_expression_validation_sig_path.tsv", "data/compound_target_probability_tox21_interpret_expression/validation_result/gamma-epsilon_0.001_0.1_beta_trt_cp_1.11uM_24h_fdr_0.05_expression_validation_sig_path.tsv", "data/compound_target_probability_tox21_interpret_expression/validation_result/gamma-epsilon_0.001_0.1_beta_trt_cp_10uM_24h_fdr_0.05_expression_validation_sig_path.tsv");	# vector that contains names of identified differentially expressed DTox path files   
out_file	<- "data/compound_target_probability_tox21_interpret_expression/validation_summary/gamma-epsilon_0.001_0.1_expression_validation_sig_path_recurrent.tsv";	# name of output anallysis result file 
tox21_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_10k_library_info.tsv";	# name of input Tox21 compound library info file 
reactome_file	<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";	# name of Reactome protein-pathway annotation file

## 1. Process Tox21 compound library info data 
# read in Tox21 compound library info as data frame 
tox21_df <- read.delim(file = tox21_file, header = T, sep = "\t");
# add "CID_" before each Pubchem CID of compound, convert compound name characters to lower case 
tox21_df$cid <- sapply(tox21_df$PUBCHEM_CID, function(tdpc) paste("CID", tdpc, sep = "_"));
tox21_df$compound_name <- tolower(tox21_df$SAMPLE_NAME);
# create unique mapping data frame between compound CIDs and names
tox21_df <- unique(tox21_df[, c("cid", "compound_name")]);
tox21_compound_name <- tox21_df$compound_name;
names(tox21_compound_name) <- tox21_df$cid;

## 2. Process Reactome annotation data  
# read in Reactome protein-pathway annotation as data frame 
reactome_df <- read.delim(file = reactome_file, header = F, sep = "\t");
# remove non-human annotation rows in the data frame 
human_id <- which(reactome_df[,6] %in% "Homo sapiens");
reactome_name_df <- unique(reactome_df[human_id, c(2, 4)]);
# create a map between reactome pathway ID and name   
reactome_name <- reactome_name_df[, 2];
names(reactome_name) <- reactome_name_df[, 1];

## 3. Analyze expression-validation results to identify recurrent differentially expressed DTox paths   
# read in identified differentially expressed DTox paths as data frames, merge data frames by row  
sig_path_df_list <- lapply(sig_path_files, function(spf) read.delim(file = spf, header = T, sep = "\t")); 
sig_path_df1 <- do.call(rbind, sig_path_df_list);
# create unique mapping between Tox21 assay, compound, and identified differentially expressed DTox paths 
sig_path_df <- unique(sig_path_df1[, c("assay", "compound", "sig_path")]);
# iterate by Tox21 assay dataset, identify differentially expressed DTox paths that recurrently appear among compounds for each assay 
all_assays <- unique(sig_path_df$assay); 
aa_sig_path <- lapply(all_assays, function(aa){
	# select differentially expressed DTox paths relevant to the current Tox21 assay 
	aa_id <- which(sig_path_df$assay %in% aa);
	aa_df <- sig_path_df[aa_id, ];
	# count the number of occurrence for each unique differentially expressed DTox path, identify recurrent DTox paths (if any)
	aa_path_table <- table(aa_df$sig_path);
	recurrent_id <- which(aa_path_table > 1);
	if(length(recurrent_id) == 0)	return(logical(0))
	else{ 
		aa_recurrent_count <- as.numeric(aa_path_table[recurrent_id])
		aa_recurrent_path <- names(aa_path_table)[recurrent_id]
		# generate strings to represent DTox recurrent paths 
		aa_recurrent_path_name <- sapply(aa_recurrent_path, function(arp){
			# map IDs of pathways along DTox paths to their names, separate names by ';'
			arp_s <- strsplit(arp, "_")[[1]]
			as_name <- reactome_name[arp_s]
			as_p <- paste(as_name, collapse = ";")
			return(as_p)
		})
		# generate strings to represent compounds for which the DTox recurrent paths were identified 	
		aa_recurrent_compounds <- sapply(aa_recurrent_path, function(arp){
			# obtain the compound names for each recurrent path, separate names by ';' 
			arp_compounds <- aa_df$compound[which(aa_df$sig_path %in% arp)]
			ac_names <- tox21_compound_name[arp_compounds]
			ac_p <- paste(ac_names, collapse = ";")
			return(ac_p)
		})
		# store Tox21 assay info, recurrent DTox paths, number of occurrence, compound names in a data frame   
		aa_recurrent_df <- data.frame(rep(aa, length(recurrent_id)), aa_recurrent_path_name, aa_recurrent_count, aa_recurrent_compounds)
		colnames(aa_recurrent_df) <- c("assay", "recurrent_path", "n_compounds", "compound_name")	
		ard_od <- order(aa_recurrent_count, decreasing = T)
		aa_recurrent_df <- aa_recurrent_df[ard_od, ]
		return(aa_recurrent_df)
	}
});
# merge result data frames from multiple Tox21 assays by row, write to output file    
asp_len <- sapply(aa_sig_path, length);
aa_sig_path <- aa_sig_path[asp_len > 0];
sig_path_df <- do.call(rbind, aa_sig_path);
write.table(sig_path_df, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F);
