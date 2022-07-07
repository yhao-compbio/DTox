# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2022
## This script analyzes the DTox-predicted HepG2 cytotoxicity scores by EPA chemical list and DrugBank approval status list 


## 0. Input arguments 
epa_map_file		<- "/home/yunhao1/project/tox_data/data/tox_list/EPA_chemical_lists_dsstox_map.tsv";
drugbank_map_file	<- "/home/yunhao1/project/tox_data/data/tox_list/Drugbank_lists_dsstox_map.tsv";
tox21_sample_file	<- "/home/yunhao1/project/tox_data/data/tox21/tox21_assay_samples.tsv";
query_assay		<- "tox21-rt-viability-hepg2-p1";
tox21_pred_file		<- "data/compound_target_probability_tox21_prediction/rt-viability-hepg2-p1/tox21_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred.tsv";
dsstox_pred_file	<- "data/compound_target_probability_tox21_prediction/rt-viability-hepg2-p1/dsstox_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred.tsv";	# name of DTox HepG2 viability prediction result file  
output_folder		<- "data/compound_target_probability_tox21_prediction/rt-viability-hepg2-p1/dsstox_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred";

## 1. Process input data file to obtain the predicted HepG2 cytotoxicity scores of Tox21 compounds
# read in Tox21 assay screening results as data frame, obtain the active/inactive compounds of query assay  
tox21_sample_df <- read.delim(file = tox21_sample_file, header = T, sep = "\t");
assay_id <- which(tox21_sample_df$assay_name %in% query_assay);
assay_sample_df <- tox21_sample_df[assay_id, ];
# read in DTox HepG2 viability model prediction result on Tox21 compounds as data frame  
tox21_pred_df <- read.delim(file = tox21_pred_file, header = T, sep = "\t");
tox21_pred_df1 <- tox21_pred_df[, c("X", "predicted_outcome")];
colnames(TOX21_pred_df1) <- c("compound_pubchem_cid", "predicted_outcome"); 
tox21_assay_pred_df <- merge(assay_sample_df, tox21_pred_df1, by = "compound_pubchem_cid");
# group the predicted HepG2 cytotoxicity scores by active/inactive compounds
pos_id <- which(tox21_assay_pred_df$assay_outcome == 1);
tox21_pos_pred_df <- tox21_assay_pred_df[pos_id, ];
neg_id <- which(tox21_assay_pred_df$assay_outcome == 0);
tox21_neg_pred_df <- tox21_assay_pred_df[neg_id, ];

## 2. Analyze the predicted HepG2 cytotoxicity scores of compounds by EPA chemical list 
# read in EPA list ~ compound mapping as data frame  
epa_map_df <- read.delim(file = epa_map_file, header = T, sep = "\t");
# read in DTox HepG2 viability model prediction result on DSSTox compounds as data frame  
dsstox_pred_df <- read.delim(file = dsstox_pred_file, header = T, sep = "\t");
dsstox_pred_df1 <- dsstox_pred_df[, c("X", "predicted_outcome")];
colnames(dsstox_pred_df1) <- c("dsstox_substance_id", "predicted_outcome");
# merge prediction data frame with compounds of EPA list 
dsstox_el_pred_df <- merge(epa_map_df, dsstox_pred_df1, by = "dsstox_substance_id");
# obtain the predicted HepG2 cytotoxicity scores of compounds in each EPA list 
epa_list_names <- unique(dsstox_el_pred_df$list_acronym);
epa_list_des <- sapply(epa_list_names, function(eln){
	eln_id <- which(dsstox_el_pred_df$list_acronym %in% eln)[[1]];
	eln_des <- dsstox_el_pred_df$list_description[[eln_id]];
	return(eln_des);
});
# iterate by EPA list of interest, compare the predicted HepG2 cytotoxicity scores of compounds to pos/neg controls  
dsstox_el_pred_compare <- mapply(function(eln){
	# select rows of interest from processed prediction data frame to obtain prediction results of DSSTox compounds associated with current list  
	eln_id <- which(dsstox_el_pred_df$list_acronym %in% eln);
	eln_scores <- dsstox_el_pred_df$predicted_outcome[eln_id];
	es_med <- median(eln_scores);
	es_len <- length(eln_scores);
	# perform one-sided Mann-Whitney U test to compare predicted HepG2 cytotoxicity scores of compounds to negative controls   
	eln_neg_pv <- wilcox.test(eln_scores, tox21_neg_pred_df$predicted_outcome, alternative = "greater")$p.value;
	# perform one-sided Mann-Whitney U test to compare predicted HepG2 cytotoxicity scores of compounds to positive controls 
	eln_pos_pv <- wilcox.test(eln_scores, tox21_pos_pred_df$predicted_outcome, alternative = "less")$p.value;
	return(c(es_len, es_med, eln_neg_pv, eln_pos_pv));
}, epa_list_names);
# store list info, computed wilcoxon test P-value in data frame (only include list with >= 10 compounds), write to output file
epa_list_result_df <- data.frame(epa_list_names, epa_list_des, t(dsstox_el_pred_compare));
colnames(epa_list_result_df) <- c("list_acronym", "list_description", "n_compound", "compound_pred_score_median", "neg_compare_wilcox_pv", "pos_compare_wilcox_pv");
elrd_od <- order(epa_list_result_df$pos_compare_wilcox_pv, epa_list_result_df$neg_compare_wilcox_pv, decreasing = c(T, F));
epa_list_result_df <- epa_list_result_df[elrd_od, ];
elrd_id <- which(epa_list_result_df$n_compound >= 10);
epa_list_result_df <- epa_list_result_df[elrd_id, ];
epa_list_result_file <- paste(output_folder, "_epa_list_compare_result_summary.tsv", sep = "");
write.table(epa_list_result_df, file = epa_list_result_file, sep = "\t", row.names = F, col.names = T, quote = F);
# select lists with no significant difference from pos controls, but significant higher scores than neg controls. output predicted scores of compounds in such lists for plotting    
eop_list_id1 <- which(epa_list_result_df$neg_compare_wilcox_pv < 0.05);
eop_list_id2 <- which(epa_list_result_df$pos_compare_wilcox_pv > 0.05); 
eop_lists <- epa_list_result_df$list_acronym[intersect(eop_list_id1, eop_list_id2)];
eop_lists <- setdiff(eop_lists, c("casmi2017", "tscawp"));
eol_id <- which(dsstox_el_pred_df$list_acronym %in% eop_lists);
eop_list_df <- dsstox_el_pred_df[eol_id, c("list_acronym", "dsstox_substance_id", "predicted_outcome")];
eld_od <- order(eop_list_df$list_acronym);
eop_list_df <- eop_list_df[eld_od, ];

## 3. Analyze the predicted HepG2 cytotoxicity scores of compounds by DrugBank status list  
# read in DrugBank list ~ compound mapping as data frame  
drugbank_map_df <- read.delim(file = drugbank_map_file, header = T, sep = "\t");
# merge prediction data frame with compounds of DrugBank list
dsstox_dl_pred_df <- merge(drugbank_map_df, dsstox_pred_df1, by = "dsstox_substance_id");
# iterate by DrugBank list of interest, compare the predicted HepG2 cytotoxicity scores of compounds to pos/neg controls   
drugbank_list_names <- unique(dsstox_dl_pred_df$list_acronym);
dsstox_dl_pred_compare <- mapply(function(dln){
	# select rows of interest from processed prediction data frame to obtain prediction results of DSSTox compounds associated with current list
	dln_id <- which(dsstox_dl_pred_df$list_acronym %in% dln);
	dln_scores <- dsstox_dl_pred_df$predicted_outcome[dln_id];
	ds_med <- median(dln_scores);
	ds_len <- length(dln_scores);
	# perform one-sided Mann-Whitney U test to compare predicted HepG2 cytotoxicity scores of compounds to negative controls  
	dln_neg_pv <- wilcox.test(dln_scores, tox21_neg_pred_df$predicted_outcome, alternative = "greater")$p.value;
	# perform one-sided Mann-Whitney U test to compare predicted HepG2 cytotoxicity scores of compounds to positive controls 
	dln_pos_pv <- wilcox.test(dln_scores, tox21_pos_pred_df$predicted_outcome, alternative = "less")$p.value;
	return(c(ds_len, ds_med, dln_neg_pv, dln_pos_pv));
}, drugbank_list_names);
# store list info, computed wilcoxon test P-value in data frame (only include list with >= 10 compounds), write to output file
drugbank_list_result_df <- data.frame(drugbank_list_names, t(dsstox_dl_pred_compare));
colnames(drugbank_list_result_df) <- c("list_acronym", "n_compound", "compound_pred_score_median", "neg_compare_wilcox_pv", "pos_compare_wilcox_pv");
dlrd_od <- order(drugbank_list_result_df$pos_compare_wilcox_pv, drugbank_list_result_df$neg_compare_wilcox_pv, decreasing = c(T, F));
drugbank_list_result_df <- drugbank_list_result_df[dlrd_od, ];
dlrd_id <- which(drugbank_list_result_df$n_compound >= 10);
drugbank_list_result_df <- drugbank_list_result_df[dlrd_id, ];
drugbank_list_result_file <- paste(output_folder, "_drugbank_list_compare_result_summary.tsv", sep = "");
write.table(drugbank_list_result_df, file = drugbank_list_result_file, sep = "\t", row.names = F, col.names = T, quote = F);
# output predicted scores of compounds in all DrugBank lists for plotting   
dop_list_df <- dsstox_dl_pred_df[, c("list_acronym", "dsstox_substance_id", "predicted_outcome")];
dld_od <- order(dop_list_df$list_acronym);
dop_list_df <- dop_list_df[dld_od, ];
# include the predicted scores of compounds of pos/neg controls for plotting. write final data frame to output file  
tox21_pos_pred_df$list_name <- "HepG2 +";
tox21_pos_op_df <- tox21_pos_pred_df[, c("list_name", "compound_pubchem_cid", "predicted_outcome")];
tox21_neg_pred_df$list_name <- "HepG2 -";
tox21_neg_op_df <- tox21_neg_pred_df[, c("list_name", "compound_pubchem_cid", "predicted_outcome")];
colnames(tox21_pos_op_df) <- colnames(tox21_neg_op_df) <- colnames(eop_list_df) <- colnames(dop_list_df) <- c("list_name", "compound_id", "predicted_cytotoxicity_score");
score_op_df <- rbind(tox21_pos_op_df, eop_list_df, dop_list_df, tox21_neg_op_df);
score_op_file <- paste(output_folder, "_by_list.tsv", sep = "");
write.table(score_op_df, file = score_op_file, sep = "\t", row.names = F, col.names = T, quote = F);

