# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script analyzes DTox module relevance scores of viability-related pathways in the context of two viability-related assays (CASP3/7 apoptosis and mitochondria toxicity), compares pathway relevance scores between active and inactive compounds, then uses survival plot to visualize the comparison 


## functions
library(survival);


## 0. Input arguments 
query_pathway_file	<- "data/compound_target_probability_tox21_interpret_viability/tox21-rt-viability-hepg2_pathways.tsv";	# name of input viability-related pathway info file  
compound_module_file	<- "data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_module_relevance.tsv";	# name of input compound DTox module relevance file  
tox21_sample_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/data/tox21/tox21_assay_samples_all.tsv";	# name of input Tox21 assay to inactive/active compound file  
out_folder		<- "data/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1";	# folder name of output result files  
plot_folder		<- "plot/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1";	# folder name of output plot files 
casp_assay		<- "tox21-casp3-hepg2-p1";	# name of assay 1 of interest: CASP3/7 apoptosis 
mito_assay		<- "tox21-mitotox-p1";	# name of assay 2 of interest: mitochondria toxicity  

## 1. Process input data files 
# read in viability-related pathways info as data frame, obtain Reactome IDs of viability-related pathways 
query_pathway_df <- read.delim(file = query_pathway_file, header = T, sep = "\t");
query_pathways <- query_pathway_df$pathway_id;
# read in inactive/active compounds of Tox21 assays as data frame   
tox21_sample_df <- read.delim(tox21_sample_file, header = T, sep = "\t");
# read in DTox module relevance of compounds as data frame
compound_module_df <- read.delim(file = compound_module_file, header = T, sep = "\t");
# extract columns of interest from the data frame to obtain module relevance scores of viability-related pathways, convert negative relevance scores to 0 (as in our DTox interpretation method) 
pathway_module <- sapply(query_pathway_df$node_id, function(qpdni) paste("X", qpdni, sep = ""));
compound_pathway_df <- compound_module_df[, pathway_module];
colnames(compound_pathway_df) <- query_pathways;
compound_pathway_df[compound_pathway_df < 0] <- 0;

## 2. Analyze DTox module relevance scores of viability-related pathways in the context of CASP3/7 apoptosis assay 
# select rows of interest from Tox21 assay-compound data frame to obtain inactive/active compounds of CASP3/7 assay 
casp_id <- which(tox21_sample_df$assay_name %in% casp_assay);
casp_sample_df <- tox21_sample_df[casp_id,]; 
# iterate by viability-related pathway, compare module relevance scores of each pathway between active compounds and inactive compounds
casp_module_compare_list <- mapply(function(qp){
	# obtain the module relevance scores of current pathway  
	ccpd_pathway_df <- data.frame(compound_module_df$X, compound_pathway_df[, qp]);
	colnames(ccpd_pathway_df) <- c("compound_pubchem_cid", "pathway_relevance_score");
	# merge pathway relevance score data frame with CASP3/7 assay-compoud data frame, obtain relevance scores of current pathway for active compounds   
	ccpd_relevance_df <- merge(casp_sample_df, ccpd_pathway_df, by = "compound_pubchem_cid"); 
	crd_pos_ids <- which(ccpd_relevance_df$assay_outcome == 1);
	crd_pos_values <- ccpd_relevance_df$pathway_relevance_score[crd_pos_ids];
	# sort obtained relevance scores, compute survival function after sorting  
	crd_pos_value_uni <- sort(unique(crd_pos_values));
	crd_pos_value_uni <- setdiff(crd_pos_value_uni, c(0, max(crd_pos_value_uni)));
	crd_pos_surv <- sapply(crd_pos_value_uni, function(cpvu) length(which(crd_pos_values > cpvu))/length(crd_pos_ids));
	# obtain relevance scores of current pathway for inactive compounds 
	crd_neg_ids <- which(ccpd_relevance_df$assay_outcome == 0);
	crd_neg_values <- ccpd_relevance_df$pathway_relevance_score[crd_neg_ids];
	# sort obtained relevance scores, compute survival function after sorting  
	crd_neg_value_uni <- sort(unique(crd_neg_values));
	crd_neg_value_uni <- setdiff(crd_neg_value_uni, c(0, max(crd_neg_value_uni)));
	crd_neg_surv <- sapply(crd_neg_value_uni, function(cnvu) length(which(crd_neg_values > cnvu))/length(crd_neg_ids));		
	# test survival difference between active compounds and inactive compounds using log-rank test, compute test P-value  
	ccpd_relevance_df$status <- 1;
	crd_surv <- survdiff(Surv(pathway_relevance_score, status) ~ assay_outcome, data = ccpd_relevance_df);	
	crd_lrpv <- 1 - pchisq(crd_surv$chisq, df = 1);		
	return(ls = list(pos_x = crd_pos_value_uni, pos_y = crd_pos_surv, neg_x = crd_neg_value_uni, neg_y = crd_neg_surv, lr_pv = crd_lrpv));
}, query_pathways, SIMPLIFY = F);
# obtain computed log-rank test p-values of all viability-related pathways, perform multiple testing correction by FDR  
casp_module_compare_lrpv <- sapply(casp_module_compare_list, function(cmcl) cmcl[["lr_pv"]]);
casp_module_compare_fdr <- p.adjust(casp_module_compare_lrpv, method = "fdr");
# store pathway P-value and FDR info in a data frame, write to output result file  
casp_compare_result_df <- data.frame(query_pathway_df$pathway_name_short, casp_module_compare_lrpv, casp_module_compare_fdr);
colnames(casp_compare_result_df) <- c("pathway", "log_rank_pvalue", "log_rank_fdr");
casp_compare_file <- paste(out_folder, "_pathway_module_compare_", casp_assay, ".tsv", sep = "");
write.table(casp_compare_result_df, casp_compare_file, sep = "\t", row.names = F, col.names = T, quote = F);

## 3. Visualize DTox module relevance score comparison in the context of CASP3/7 apoptosis assay  
# specify figure and font size of survival plot
casp_plot_file <- paste(plot_folder, "_pathway_module_compare_", casp_assay, ".pdf", sep = "");
pdf(casp_plot_file, width = 9, height = 9, family = "Helvetica"); 
layout(matrix(1:9, 3, 3, byrow = TRUE));
par(cex.lab = 1.8, cex.axis = 1.8, cex.main = 2);
par(mar = c(4, 4.5, 2, 0.5));
# set y-axis range of survival plot 
all_casp_y <- unlist(lapply(casp_module_compare_list, function(cmcl) cmcl[c(2,4)]));
casp_y_min <- floor(min(all_casp_y) * 1000)/1000;
casp_y_max <- ceiling(max(all_casp_y) * 10)/10;
# iterate by viability-related pathway, make survival plot showing the pathway relevance score comparison
for (lcmcl in 1:length(casp_module_compare_list)){
	# set x-axis range 
	cmcl <- casp_module_compare_list[[lcmcl]];
	cmcl_x <- unlist(list(cmcl$pos_x, cmcl$neg_x));
	cmcl_x_min <- 10^(floor(log10(min(cmcl_x))));
	cmcl_x_max <- 10^(ceiling(log10(max(cmcl_x))));
	# make survival plot showing the relevance score comparison of current pathway between active compounds and inactive compounds  
	plot(cmcl$pos_x, cmcl$pos_y, 
		main = query_pathway_df$pathway_name_short[[lcmcl]], 
		log = "xy", 
		xlim = c(cmcl_x_min, cmcl_x_max), 
		ylim = c(casp_y_min, casp_y_max), 
		xlab = "",
		ylab = "",
		axes = F,
		cex = 0.4,
		pch = 16,
		col = "darkorange")
	# set x-axis label and ticks 
	N_xt <- log10(cmcl_x_max/cmcl_x_min)/2;
	xtick <- cmcl_x_min * 10^(0:N_xt*2);
	if(max(xtick) < cmcl_x_max)	xtick <- c(xtick, cmcl_x_max)
	xtick_label <- sapply(xtick, function(x){
		xl <- log10(x);
		if((xl < 0) || (xl > 1))	return(parse(text = paste("10^", xl, sep = "")))
		else	return(x)	
	});
	axis(side = 1, at = xtick, labels = xtick_label);
	mtext("Pathway relevance", side = 1, line = 2.5, cex = 1.35);	
	# set y-axis label and ticks 
	ytick <- c(1e-3, 1e-2, 1e-1, 1);
	ytick_l <- c(expression(10^-3), expression(10^-2), expression(10^-1), 1);
	axis(side = 2, at = ytick, labels = ytick_l);
	points(cmcl$neg_x, cmcl$neg_y, 
		cex = 0.4,
		pch = 16,
		col = "darkgray")
	mtext("Survival function", side = 2, line = 2.7, cex = 1.35);
	# add text specifying the FDR of current pathway, set color to red if significant  
	cmcl_fdr <- casp_module_compare_fdr[[lcmcl]]
	if(cmcl_fdr < 0.01)	fdr_text <- paste("FDR = ", formatC(cmcl_fdr, format = "e", digit = 2), sep = "")
	else	fdr_text <- paste("FDR = ", round(cmcl_fdr, 2), sep = "")
	if(cmcl_fdr < 0.05)	fc <- "red" 
	else	fc <- "black"
	mtext(fdr_text, side = 1, line = -1.5, adj = 0.1, col = fc, cex = 1.2);
	# add figure legend  
	if(lcmcl == 1){
		legend(5e-7, 0.015, c("caspase 3/7+", "caspase 3/7-"), pch = 16, cex = 1.8, col = c("darkorange", "darkgray"), bty = "n")		
	}
}
dev.off();

## 4. Analyze DTox module relevance scores of viability-related pathways in the context of mitochondria toxicity assay
# select rows of interest from Tox21 assay-compound data frame to obtain inactive/active compounds of mitochondria toxicity assay
mito_id <- which(tox21_sample_df$assay_name %in% mito_assay);
mito_sample_df <- tox21_sample_df[mito_id,];
# iterate by viability-related pathway, compare module relevance scores of each pathway between active compounds and inactive compounds
mito_module_compare_list <- mapply(function(qp){
        #  obtain the module relevance scores of current pathway  
	ccpd_pathway_df <- data.frame(compound_module_df$X, compound_pathway_df[, qp]); 
	colnames(ccpd_pathway_df) <- c("compound_pubchem_cid", "pathway_relevance_score");
	# merge pathway relevance score data frame with mitotox assay-compoud data frame, obtain relevance scores of current pathway for active compounds 
	ccpd_relevance_df <- merge(casp_sample_df, ccpd_pathway_df, by = "compound_pubchem_cid"); 
	crd_pos_ids <- which(ccpd_relevance_df$assay_outcome == 1);
	crd_pos_values <- ccpd_relevance_df$pathway_relevance_score[crd_pos_ids];
	# sort obtained relevance scores, compute survival function after sorting 
	crd_pos_value_uni <- sort(unique(crd_pos_values));
	crd_pos_value_uni <- setdiff(crd_pos_value_uni, c(0, max(crd_pos_value_uni)));
	crd_pos_surv <- sapply(crd_pos_value_uni, function(cpvu) length(which(crd_pos_values > cpvu))/length(crd_pos_ids));
	# obtain relevance scores of current pathway for inactive compounds 
	crd_neg_ids <- which(ccpd_relevance_df$assay_outcome == 0);
	crd_neg_values <- ccpd_relevance_df$pathway_relevance_score[crd_neg_ids];
	# sort obtained relevance scores, compute survival function after sorting  
	crd_neg_value_uni <- sort(unique(crd_neg_values));
	crd_neg_value_uni <- setdiff(crd_neg_value_uni, c(0, max(crd_neg_value_uni)));
	crd_neg_surv <- sapply(crd_neg_value_uni, function(cnvu) length(which(crd_neg_values > cnvu))/length(crd_neg_ids));
	# test survival difference between active compounds and inactive compounds using log-rank test, compute test P-value  
	ccpd_relevance_df$status <- 1;
	crd_surv <- survdiff(Surv(pathway_relevance_score, status) ~ assay_outcome, data = ccpd_relevance_df);
	crd_lrpv <- 1 - pchisq(crd_surv$chisq, df = 1);
	return(ls = list(pos_x = crd_pos_value_uni, pos_y = crd_pos_surv, neg_x = crd_neg_value_uni, neg_y = crd_neg_surv, lr_pv = crd_lrpv));
}, query_pathways, SIMPLIFY = F);
# obtain computed log-rank test p-values of all viability-related pathways, perform multiple testing correction by FDR
mito_module_compare_lrpv <- sapply(mito_module_compare_list, function(mmcl) mmcl[["lr_pv"]]);
mito_module_compare_fdr <- p.adjust(mito_module_compare_lrpv, method = "fdr");
# store pathway P-value and FDR info in a data frame, write to output result file
mito_compare_result_df <- data.frame(query_pathway_df$pathway_name_short, mito_module_compare_lrpv, mito_module_compare_fdr);
colnames(mito_compare_result_df) <- c("pathway", "log_rank_pvalue", "log_rank_fdr");
mito_compare_file <- paste(out_folder, "_pathway_module_compare_", mito_assay, ".tsv", sep = "");
write.table(mito_compare_result_df, mito_compare_file, sep = "\t", row.names = F, col.names = T, quote = F); 

## 5. Visualize DTox module relevance score comparison in the context of mitochondria toxicity assay
# specify figure and font size of survival plot
mito_plot_file <- paste(plot_folder, "_pathway_module_compare_", mito_assay, ".pdf", sep = "");
pdf(mito_plot_file, width = 9, height = 9, family = "Helvetica");
layout(matrix(1:9, 3, 3, byrow = TRUE));
par(cex.lab = 1.8, cex.axis = 1.8, cex.main = 2);
par(mar = c(4, 4.5, 2, 0.5));
# set y-axis range of survival plot 
all_mito_y <- unlist(lapply(mito_module_compare_list, function(mmcl) mmcl[c(2,4)]));
mito_y_min <- floor(min(all_mito_y) * 1000)/1000;
mito_y_max <- ceiling(max(all_mito_y) * 10)/10; 
# iterate by viability-related pathway, make survival plot showing the pathway relevance score comparison
for (lmmcl in 1:length(mito_module_compare_list)){
	# set x-axis range 
	mmcl <- mito_module_compare_list[[lmmcl]];
	mmcl_x <- unlist(list(mmcl$pos_x, mmcl$neg_x));
	mmcl_x_min <- 10^(floor(log10(min(mmcl_x))));
	mmcl_x_max <- 10^(ceiling(log10(max(mmcl_x))));
	# make survival plot showing the relevance score comparison of current pathway between active compounds and inactive compounds 
	plot(mmcl$pos_x, mmcl$pos_y,
		main = query_pathway_df$pathway_name_short[[lmmcl]],
		log = "xy",
		xlim = c(mmcl_x_min, mmcl_x_max),
		ylim = c(mito_y_min, mito_y_max),
		xlab = "",
		ylab = "",
		axes = F,
		cex = 0.4,
		pch = 16,
		col = "darkorange")
	# set x-axis label and ticks 
	N_xt <- log10(mmcl_x_max/mmcl_x_min)/2;
	xtick <- mmcl_x_min * 10^(0:N_xt*2);
	if(max(xtick) < mmcl_x_max)	xtick <- c(xtick, mmcl_x_max)
	xtick_label <- sapply(xtick, function(x){
		xl <- log10(x);
		if((xl < 0) || (xl > 1))	return(parse(text = paste("10^", xl, sep = "")))
		else	return(x)
	});
	axis(side = 1, at = xtick, labels = xtick_label);
	mtext("Pathway relevance", side = 1, line = 2.5, cex = 1.35);
	# set y-axis label and ticks 
	ytick <- c(2e-3, 1e-2, 1e-1, 1);
	ytick_l <- c(0.002, expression(10^-2), expression(10^-1), 1);
	axis(side = 2, at = ytick, labels = ytick_l);
	points(mmcl$neg_x, mmcl$neg_y,
		cex = 0.4,
		pch = 16,
		col = "darkgray")
	mtext("Survival function", side = 2, line = 2.7, cex = 1.35);
	# add text specifying the FDR of current pathway, set color to red if significant   
	mmcl_fdr <- mito_module_compare_fdr[[lmmcl]]
	if(mmcl_fdr < 0.01)	fdr_text <- paste("FDR = ", formatC(mmcl_fdr, format = "e", digit = 0), sep = "")
	else	fdr_text <- paste("FDR = ", round(mmcl_fdr, 2), sep = "")
	if(mmcl_fdr < 0.05)	fc <- "red"
	else	fc <- "black"
	mtext(fdr_text, side = 1, line = -1.5, adj = 0.1, col = fc, cex = 1.2);
	# add figure legend
	if(lmmcl == 1){
		legend(5e-7, 0.024, c("mito toxicity+", "mito toxicity-"), pch = 16, cex = 1.8, col = c("darkorange", "darkgray"), bty = "n")
	}
}
dev.off();

