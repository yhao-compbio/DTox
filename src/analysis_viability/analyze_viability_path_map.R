# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script analyzes viability-related DTox paths from model interpretation results in the context of drug-induced liver injury (DILI) adverse events and ATC drug classification, evaluates the enrichment of DILI events/ATC drug classes among compounds identified with viability-related DTox paths, then visualizes the relationships between viability-related DTox paths and DILI events/ATC drug classes by heatmap 


## functions
source("src/functions.R");
library(RColorBrewer);
library(pals);
library(pheatmap);
library(gplots);


## 0. Input arguments 
query_pathway_file	<- "data/compound_target_probability_tox21_interpret_viability/tox21-rt-viability-hepg2_pathways.tsv";	# name of input viability-related pathway info file
compound_path_file	<- "data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_path_relevance_pv.tsv";	# name of input HepG2 cell vibiality DTox model interpretation result file  
compound_dili_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/data/LiverTox/dili_phenotype_tox21_compounds.tsv";	# name of input DILI Tox21 compounds file  
compound_atc_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/data/tox21/tox21_10k_library_atc_code.tsv";	# name of input Tox21 compound ATC code file 
plot_folder		<- "plot/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1";	# folder name of output plot files

## 1. Process input viability-related pathway info and vibiality DTox model interpretation result 
# read in viability-related pathways info as data frame, obtain Reactome IDs of viability-related pathways  
query_pathway_df <- read.delim(file = query_pathway_file, header = T, sep = "\t");
query_pathways <- query_pathway_df$pathway_id;
# read in significant DTox paths of compounds for DTox viability assay model as data frame, group significant DTox path by compounds   
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
compound_path <- group.vector.by.categories(compound_path_df$cid, compound_path_df$path_id);
# iterate by compound, check whether viability-related pathways are along the significant DTox paths of each compound, store result in a 2D matrix map     
compound_path_map <- mapply(function(cp){
	# separate string by underscore to obtain the pathways along each significant DTox paths identified for the current compound
	cp_s <- lapply(cp, function(cpp) strsplit(cpp, "_")[[1]]);
	cp_s <- unique(unlist(cp_s));
	# check overlap with each viability-related pathway 
	qp_inter <- sapply(query_pathways, function(qp) length(intersect(qp, cp_s)));
	names(qp_inter) <- query_pathways;
	return(qp_inter);       
}, compound_path);
rownames(compound_path_map) <- query_pathway_df$pathway_name_short;
path_count <- rowSums(compound_path_map);
N_compound <- ncol(compound_path_map);

## 2. Visualize the relationships between compound and viability-related pathway in DTox model by heatmap 
# specify figure and font size 
heat_file <- paste(plot_folder, "_path_compound_map.pdf", sep = "");
pdf(heat_file, width = 12, height = 2.8, family = "Helvetica");
# make heatmap showing the relationships between compound and viability-related pathway in DTox model 
heat_results <- pheatmap(compound_path_map,
	color = colorRampPalette(brewer.pal(n = 3, name = "Purples"))(2),
	border_color = "black",
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	clustering_distance_cols = "binary",
	clustering_method = "average",
	show_rownames = T,
	show_colnames = F,
	fontsize = 15,
	legend = F
)
dev.off();

## 3. Perform enrichment analysis to study the relationships between DILI events and viability-related pathways 
# read in Tox21 compounds of DILI adverse event as data frame  
compound_dili_df <- read.delim(file = compound_dili_file, header = T, sep = "\t");
# iterate by DILI event, analyze the enrichment of DILI compounds among compounds of viability-related pathways
compound_dili_types <- unique(compound_dili_df$phenotype_name_short);
compound_dili_summary_list <- lapply(compound_dili_types, function(cdt){
	# select rows of interest from DILI compound data frame to obtain Tox21 compounds associated with the current DILI event 
	cdt_ids <- which(compound_dili_df$phenotype_name_short %in% cdt);
	cdt_drugs <- unique(compound_dili_df$cid[cdt_ids]);
	# iterate by viability-related pathway, test whether compounds identified with each pathway are enriched for current DILI event  
	N_cdt <- length(cdt_drugs);
	cdt_summary <- mapply(function(rcpm){
		# build 2*2 contingency table, perform Fisher's exact test to evaluate the enrichment of current DILI event among compounds identified with current pathway   
		N_rcpm_cdt <- sum(compound_path_map[rcpm, cdt_drugs]);
		N_ircpm_cdt <- N_cdt - N_rcpm_cdt;
		N_rcpm_icdt <- path_count[[rcpm]] - N_rcpm_cdt;
		N_ii <- N_compound - N_rcpm_icdt - N_cdt;
		conti_mat <- matrix(c(N_rcpm_cdt, N_ircpm_cdt, N_rcpm_icdt, N_ii), 2, 2, byrow = T);
		conti_test <- fisher.test(conti_mat, alternative = "greater");
		# compute odds ratio and P-value of the Fisher's exact test 
		cdt_rcpm_ratio <- N_rcpm_cdt/N_cdt;
		return(c(cdt_rcpm_ratio, conti_test$estimate, conti_test$p.value));
	}, rownames(compound_path_map));
	# store computed odds ratio and P-value in a data frame, along with DILI events and viability-related pathway info
	cdt_summary_df <- data.frame(cdt, rownames(compound_path_map), t(cdt_summary));
	colnames(cdt_summary_df) <- c("dili_phenotype", "pathway", "compound_path_ratio", "odds_ratio", "fisher_pv");
	return(cdt_summary_df);
});
# aggregate enrichment analysis results from all DILI events into one data frame, perform multiple testing correction by FDR
compound_dili_summary_df <- do.call(rbind, compound_dili_summary_list);
compound_dili_summary_df$fdr <- p.adjust(compound_dili_summary_df$fisher_pv, method = "fdr");
# aggregate computed odds ratio to form relationship map between DILI events and viability-related pathways  
dili_path_map <- mapply(function(cdsl) cdsl$odds_ratio, compound_dili_summary_list);
rownames(dili_path_map) <- rownames(compound_path_map);
colnames(dili_path_map) <- compound_dili_types;
# perform log transformation to the odds ratio map matrix  
dili_path_map1 <- log(dili_path_map);
dpm_max <- max(abs(dili_path_map1[abs(dili_path_map1) != Inf]));
dili_path_map1[dili_path_map1 == Inf] <- dpm_max;
dili_path_map1[dili_path_map1 == -Inf] <- -dpm_max;

## 4. Visualize odds ratio map between DILI events and viability-related pathways by heatmap 
# set x-axis and y-axis range 
x_l <- 0.128
x_u <- 0.792
y_l <- 0.405
y_u <- 0.905
# identify the x and y coordinates for cells representing significant DILI event-pathway pairs (FDR < 0.05) 
dili_sig_loc <- mapply(function(cdsdf, lcdsdf){
	if(cdsdf < 0.05){
		# obtain the column and row indexes of current cell  
		sig_cid <- ceiling(lcdsdf/nrow(dili_path_map1))
		sig_rid <- nrow(dili_path_map1) + 1 - (lcdsdf %% nrow(dili_path_map1))
		if(sig_rid == nrow(dili_path_map1) + 1)	sig_rid <- 1
		# compute x and y coordinates based on indexes 
		sig_x <- x_l + (sig_cid - 0.5) * (x_u - x_l)/ncol(dili_path_map1)	
		sig_y <- y_l + (sig_rid - 0.5) * (y_u - y_l)/nrow(dili_path_map1)	
		return(c(sig_x, sig_y))
	}
	else	return(c(NA, NA))
}, compound_dili_summary_df$fdr, 1:length(compound_dili_summary_df$fdr));
nna_id <- which(is.na(dili_sig_loc[1, ]) == F);
dili_sig_loc <- dili_sig_loc[, nna_id];
# specify figure and font size  
dili_heat_file <- paste(plot_folder, "_dili_path_map.pdf", sep = "");
pdf(dili_heat_file, width = 12, height = 7.5, family = "Helvetica");
# make heatmap showing the relationships between DILI events and viability-related pathways  
heatmap.2(dili_path_map1,
	Rowv = F, 
	Colv = F,
	dendrogram = c("none"),
	density.info = "none",
	trace = "none",
	symm = F,
	symbreaks = T,
	scale = c("none"),
	col = colorRampPalette(c("light blue","white","red"))(n = 200),
	cexRow = 2,
	cexCol = 2,
	adjCol = 1,
	key = T,
	keysize = 0.9,
	key.par = list(cex.lab = 1.5, cex.axis = 1.25),
	key.title = NA,
	symkey = T,
	key.xlab = "Odds Ratio (log)",	
	margins = c(20, 16),
	colsep = 0:ncol(dili_path_map1),
	rowsep = 0:nrow(dili_path_map1),
	sepcolor = "black",
	sepwidth = c(0,0)
)
# add stars to the cells that represent significant DILI event-pathway pairs 
points(dili_sig_loc[1, ], dili_sig_loc[2, ], pch = 8, cex = 1, col = "black", xpd = T);
dev.off();

## 5. Perform enrichment analysis to study the relationships between drug ATC classes and viability-related pathways 
# read in ATC classes of Tox21 compounds as data frame  
compound_atc_df <- read.delim(file = compound_atc_file, header = T, sep = "\t");
# iterate by ATC class, analyze the enrichment of drug classes among compounds of viability-related pathways
compound_atc_types <- unique(compound_atc_df$atc_id1_annotation);
compound_atc_summary_list <- lapply(compound_atc_types, function(catt){
	# select rows of interest from compound ATC class data frame to obtain Tox21 compounds associated with the current drug class 
	catt_ids <- which(compound_atc_df$atc_id1_annotation %in% catt);
	catt_drugs <- unique(compound_atc_df$cid[catt_ids]); 
	# iterate by viability-related pathway, test whether compounds identified with each pathway are enriched for current ATC class
	N_catt <- length(catt_drugs);
	catt_summary <- mapply(function(rcpm){
		# build 2*2 contingency table, perform Fisher's exact test to evaluate the enrichment of current ATC class among compounds identified with current pathway 
		N_rcpm_catt <- sum(compound_path_map[rcpm, catt_drugs]);
		N_ircpm_catt <- N_catt - N_rcpm_catt;
		N_rcpm_icatt <- path_count[[rcpm]] - N_rcpm_catt;
		N_ii <- N_compound - N_rcpm_icatt - N_catt;
		conti_mat <- matrix(c(N_rcpm_catt, N_ircpm_catt, N_rcpm_icatt, N_ii), 2, 2, byrow = T);
		conti_test <- fisher.test(conti_mat, alternative = "greater");
		# compute odds ratio and P-value of the Fisher's exact test 
		catt_rcpm_ratio <- N_rcpm_catt/N_catt;
		return(c(catt_rcpm_ratio, conti_test$estimate, conti_test$p.value));
	}, rownames(compound_path_map));
	# store computed odds ratio and P-value in a data frame, along with ATC classes and viability-related pathway info
	catt_summary_df <- data.frame(catt, rownames(compound_path_map), t(catt_summary));
	colnames(catt_summary_df) <- c("atc_class", "pathway", "compound_path_ratio", "odds_ratio", "fisher_pv");
	return(catt_summary_df);
});
# aggregate enrichment analysis results from all ATC classes into one data frame, perform multiple testing correction by FDR
compound_atc_summary_df <- do.call(rbind, compound_atc_summary_list);
compound_atc_summary_df$fdr <- p.adjust(compound_atc_summary_df$fisher_pv, method = "fdr");
# aggregate computed odds ratio to form relationship map between ATC classes and viability-related pathways
atc_path_map <- mapply(function(casl) casl$odds_ratio, compound_atc_summary_list);
rownames(atc_path_map) <- rownames(compound_path_map);
colnames(atc_path_map) <- compound_atc_types;
# perform log transformation to the odds ratio map matrix  
atc_path_map1 <- log(atc_path_map);
apm_max <- max(abs(atc_path_map1[abs(atc_path_map1) != Inf]));
atc_path_map1[atc_path_map1 == Inf] <- apm_max;
atc_path_map1[atc_path_map1 == -Inf] <- -apm_max;

## 6. Visualize odds ratio map between ATC classes and viability-related pathways by heatmap
# set x-axis and y-axis range 
x_l <- 0.128
x_u <- 0.577
y_l <- 0.405
y_u <- 0.905
# identify the x and y coordinates for cells representing significant ATC class-pathway pairs (FDR < 0.05) 
atc_sig_loc <- mapply(function(casdf, lcasdf){
	if(casdf < 0.05){
		# obtain the column and row indexes of current cell  
		sig_cid <- ceiling(lcasdf/nrow(atc_path_map1))
		sig_rid <- nrow(atc_path_map1) + 1 - (lcasdf %% nrow(atc_path_map1))
		if(sig_rid == nrow(atc_path_map1) + 1)	sig_rid <- 1
		# compute x and y coordinates based on indexes 
		sig_x <- x_l + (sig_cid - 0.5) * (x_u - x_l)/ncol(atc_path_map1)
		sig_y <- y_l + (sig_rid - 0.5) * (y_u - y_l)/nrow(atc_path_map1)
		return(c(sig_x, sig_y))
	}
	else	return(c(NA, NA))
}, compound_atc_summary_df$fdr, 1:length(compound_atc_summary_df$fdr));
nna_id <- which(is.na(atc_sig_loc[1, ]) == F);
atc_sig_loc <- atc_sig_loc[, nna_id];
# specify figure and font size 
atc_heat_file <- paste(plot_folder, "_atc_path_map.pdf", sep = "");
pdf(atc_heat_file, width = 12, height = 7.5, family = "Helvetica");
# make heatmap showing the relationships between ATC classes and viability-related pathways 
heatmap.2(atc_path_map1,
	Rowv = F,
	Colv = F,
	dendrogram = c("none"),
	density.info = "none",
	trace = "none",
	symm = F,
	symbreaks = T,
	scale = c("none"),
	col = colorRampPalette(c("light blue","white","red"))(n = 200),
	cexRow = 2,
	cexCol = 2,
	adjCol = 1,
	key = T,
	keysize = 0.9,
	key.par = list(cex.lab = 1.5, cex.axis = 1.25),
	key.title = NA,
	symkey = T,
	key.xlab = "Odds Ratio (log)",
	margins = c(20, 30),
	colsep = 0:ncol(atc_path_map1),
	rowsep = 0:nrow(atc_path_map1),
	sepcolor = "black",
	sepwidth = c(0,0)
)
# add stars to the cells that represent significant ATC class-pathway pairs
points(atc_sig_loc[1, ], atc_sig_loc[2, ], pch = 8, cex = 1, col = "black", xpd = T);
dev.off();
