# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script analyzes the DTox HepG2 cell viability model prediction results on DSSTox compounds, compares the predicted outcome probability between positive and negative DSSTox compounds associated with drug-induced livery injury (DILI) phenotypes


## 0. Input arguments 
dsstox_dili_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/data/LiverTox/dili_phenotype_dsstox_compounds.tsv";	# name of input DILI DSSTox compounds file
dsstox_pred_file	<- "data/compound_target_probability_tox21_prediction/rt-viability-hepg2-p1/dsstox_compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model_pred.tsv";	# name of DTox HepG2 viability prediction result file  
dili_assoc_file		<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/data/LiverTox/dili_phenotype_tox21-rt-viability-hepg2-p1_association.tsv";	# name of input DILI phenotype-HepG2 viability odds ratio file 
out_file		<- "data/compound_target_probability_tox21_prediction/rt-viability-hepg2-p1/dili_phenotype_tox21-rt-viability-hepg2-p1_dsstox_compound_model_pred";	# name of output result file 

## 1. Process input data file to obtain the predicted outcome probability of DSSTox compounds associated with DILI phenotypes 
# read in DILI phenotype-HepG2 viability association odds ratio as data frame 
dili_assoc_df <- read.delim(file = dili_assoc_file, header = T, sep = "\t");
dili_term_df <- dili_assoc_df[, c("term_name", "term_short")];
# read in DSSTox compounds of DILI phenotypes as data frame 
dsstox_dili_df <- read.delim(file = dsstox_dili_file, header = T, sep = "\t");
ddd_id <- which(dsstox_dili_df$adverse_event_name %in% dili_assoc_df$term_name);
dsstox_dili_df <- dsstox_dili_df[ddd_id, ];
# read in DTox HepG2 viability model prediction result on DSSTox compounds as data frame  
dsstox_pred_df <- read.delim(file = dsstox_pred_file, header = T, sep = "\t");
dsstox_pred_df1 <- dsstox_pred_df[, c("X", "predicted_outcome")];
colnames(dsstox_pred_df1) <- c("dsstox_substance_id", "predicted_outcome");
# merge prediction data frame with DSSTox compounds of DILI phenotypes, to obtain the predicted outcome probability of DSSTox compounds associated with DILI phenotypes   
dsstox_dili_pred_df <- merge(dsstox_dili_df, dsstox_pred_df1, by = "dsstox_substance_id");
dsstox_dili_pred_df1 <- merge(dsstox_dili_pred_df, dili_term_df, by.x = "adverse_event_name", by.y = "term_name");
dsstox_dili_pred_df1 <- dsstox_dili_pred_df1[, c("term_short", "preferred_name", "toxicity_label", "predicted_outcome")];
colnames(dsstox_dili_pred_df1) <- c("adverse_event", "drug", "toxicity_label", "predicted_outcome");
# select DSSTox compounds positively assocaited with DILI phenotypes, label them as "+"
ddpd_pos_id <- which(dsstox_dili_pred_df1$toxicity_label == 1);
dsstox_dili_pred_df1$toxicity_label[ddpd_pos_id] <- "+";
# select DSSTox compounds negatively assocaited with DILI phenotypes, label them as "-" 
ddpd_neg_id <- which(dsstox_dili_pred_df1$toxicity_label == 0);
dsstox_dili_pred_df1$toxicity_label[ddpd_neg_id] <- "-";
# write processed prediction results of DSSTox compounds associated with DILI phenotypes to output file
dsstox_dili_pred_file <- paste(out_file, "_score.tsv", sep = "");
write.table(dsstox_dili_pred_df1, file = dsstox_dili_pred_file, sep = "\t", row.names = F, col.names = T, quote = F);

## 2. Compare the predicted outcome probability of positive vs negative DSSTox compounds of DILI phenotypes
# iterate by DILI phenotype, perform wilcoxon test to compare the predicted outcome probability of positive vs negative DSSTox compounds for each phenotype 
dsstox_dili_pred_compare <- sapply(dili_assoc_df$term_short, function(dadtn){
	# select rows of interest from processed prediction data frame to obtain prediction results of DSSTox compounds associated with current phenotype 
	dadtn_id <- which(dsstox_dili_pred_df1$adverse_event %in% dadtn);
	dadtn_df <- dsstox_dili_pred_df1[dadtn_id, ];
	# obtain predicted outcome probability of DSSTox compounds positively associated with current phenotype  
	dadtn_pos_ids <- which(dadtn_df$toxicity_label %in% "+");
	dadtn_pos_scores <- dadtn_df$predicted_outcome[dadtn_pos_ids];
	# obtain predicted outcome probability of DSSTox compounds negatively associated with current phenotype 
	dadtn_neg_ids <- which(dadtn_df$toxicity_label %in% "-");
	dadtn_neg_scores <- dadtn_df$predicted_outcome[dadtn_neg_ids];
	# perform wilcoxon test to evaluate whether predicted outcome probability of positive DSSTox compounds is significantly higher than that of negative DSSTox compounds, compute P-value 
	dadtn_pv1 <- wilcox.test(dadtn_pos_scores, dadtn_neg_scores, alternative = "greater")$p.value;
	return(dadtn_pv1);
});
# store DILI phenotypes info, HepG2 viability-association odds ratio, and computed wilcoxon test P-value in data frame  
dsstox_dili_pred_compare_df <- data.frame(dili_assoc_df$term_short, dili_assoc_df$viability_odds_ratio, dili_assoc_df$viability_odds_ratio_95ci_lower, dili_assoc_df$viability_odds_ratio_95ci_upper, dsstox_dili_pred_compare);
colnames(dsstox_dili_pred_compare_df) <- c("adverse_event", "viability_odds_ratio", "viability_odds_ratio_95ci_lower", "viability_odds_ratio_95ci_upper", "prediction_compare_greater_pv");
# sort result data frame by Hepg2 viability-association odds ratio, write to output file 
ddpcd_od <- order(dsstox_dili_pred_compare_df$viability_odds_ratio, decreasing = T);
dsstox_dili_pred_compare_df <- dsstox_dili_pred_compare_df[ddpcd_od, ];
dsstox_dili_pred_compare_file <- paste(out_file, "_score_compare.tsv", sep = "");
write.table(dsstox_dili_pred_compare_df, file = dsstox_dili_pred_compare_file, sep = "\t", row.names = F, col.names = T, quote = F);
