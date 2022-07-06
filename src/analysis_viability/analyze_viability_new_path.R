# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2022
## This script identifies the prevalent target proteins and lowest level pathways that are along the paths of cytotoxic compounds not linked to the viability-related pathways by DTox


## functions
source("src/functions.R");


## 0. Input arguments 
query_pathway_file      <- "data/compound_target_probability_tox21_interpret_viability/tox21-rt-viability-hepg2_pathways.tsv";  # name of input viability-related pathway info file
compound_path_file      <- "data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_path_relevance_pv.tsv";   # name of input HepG2 cell vibiality DTox model interpretation result file  
reactome_file		<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";
id_map_file		<- "/home/yunhao1/project/target/downloads/id_map/hgnc_gene_names.tsv";
output_folder		<- "data/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1";

## 1. Process input viability-related pathway info and vibiality DTox model interpretation result 
# read in viability-related pathways info as data frame, obtain Reactome IDs of viability-related pathways  
query_pathway_df <- read.delim(file = query_pathway_file, header = T, sep = "\t");
query_pathways <- query_pathway_df$pathway_id;
# read in significant DTox paths of cytotoxic compounds for HepG2 viability assay model as data frame, group significant DTox path by compounds   
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
compound_path <- group.vector.by.categories(compound_path_df$cid, compound_path_df$path_id);
# iterate by compound, check whether viability-related pathways are along the significant DTox paths of each compound, store result in a 2D matrix map     
compound_path_inter <- sapply(compound_path, function(cp){
	# separate string by underscore to obtain the pathways along each significant DTox paths identified for the current compound
        cp_s <- lapply(cp, function(cpp) strsplit(cpp, "_")[[1]]);
	cp_s <- unique(unlist(cp_s));
	# check overlap with each viability-related pathway 
	qp_inter <- sapply(query_pathways, function(qp) length(intersect(qp, cp_s)));
	qp_inter_count <- sum(qp_inter);
	return(qp_inter_count);
});

## 2. obtain the target proteins and lowest level pathways that are along the paths of cytotoxic compounds not linked to the viability-related pathways
# obtain the cytotoxic compounds that are not linked to the viability-related pathways 
no_inter_compounds <- names(compound_path_inter)[compound_path_inter == 0];
nic_path <- compound_path[no_inter_compounds];
# obtain the target proteins and lowest level pathways that are along the paths of such cytotoxic compounds  
nic_path_lt <- lapply(nic_path, function(np){
	np_lt <- mapply(function(npp){
		npp_s <- strsplit(npp, "_")[[1]];
		ns_len <- length(npp_s);
		ns_lowest <- paste(npp_s[[2]], npp_s[[ns_len-1]], sep = "_");
		ns_target <- paste(npp_s[[2]], npp_s[[ns_len]], sep = "_"); 			
		return(c(ns_lowest, ns_target));
	}, np);
	np_lowest <- unique(np_lt[1, ]);
	np_target <- unique(np_lt[2, ]);
	return(ls = list(np_lowest, np_target));
});
# compute the prevalence of lowest level pathways among such cytotoxic compounds, rank from highest to lowest  
nic_path_lowest <- lapply(nic_path_lt, function(npl) npl[[1]]);
npl_table <- sort(table(unlist(nic_path_lowest)), decreasing = T);
npl_ratio <- npl_table/length(nic_path);
# compute the prevalence of target proteins among such cytotoxic compounds, rank from highest to lowest  
nic_path_target <- lapply(nic_path_lt, function(npl) npl[[2]]);
npt_table <- sort(table(unlist(nic_path_target)), decreasing = T);
npt_ratio <- npt_table/length(nic_path);

## 3. Process the identified lowest level pathways among such cytotoxic compounds
# read in pathway ID ~ name map from Reactome  
reactome_df <- read.delim(file = reactome_file, header = F, sep = "\t");
reactome_df <- unique(reactome_df[ ,c("V2", "V4")]);
reactome_path_name <- reactome_df$V4;
names(reactome_path_name) <- reactome_df$V2;
# map pathway ID of identified lowest level pathways (and its general pathway category) to pathway name 
npl_top_paths <- mapply(function(nt) strsplit(nt, "_")[[1]], names(npl_table));
npl_df <- data.frame(npl_top_paths[2, ], reactome_path_name[npl_top_paths[2, ]], as.integer(npl_table), as.numeric(npl_ratio), npl_top_paths[1, ], reactome_path_name[npl_top_paths[1, ]]);
# output the identified lowest level pathways, along with other info 
colnames(npl_df) <- c("lowest_level_pathway_id", "lowest_level_pathway_name", "path_compound_count", "path_compound_ratio", "general_pathway_id", "general_pathway_name");
npl_file <- paste(output_folder, "_compound_new_pathways_summary.tsv", sep = "");
write.table(npl_df, file = npl_file, sep = "\t", col.names = T, row.names = F, quote = F);

## 4. Process the identified target proteins among such cytotoxic compounds
# read in uniprot-to-symbol mapping data from HGNC 
id_map <- read.delim(file = id_map_file, header = T, sep = "\t");
# remove gene symbols that cannot be mapped to a UniProt ID 
matched_id <- which(id_map$UniProt.ID.supplied.by.UniProt. != "");
id_map <- id_map[matched_id, ];
# obtain mapped symbol IDs of each Uniprot ID 
id_uni_ids <- lapply(id_map$UniProt.ID.supplied.by.UniProt., function(idu) strsplit(idu, ", ")[[1]]);
id_uni_len <- sapply(id_uni_ids, length);
# build full map of Uniprot ID ~ symbol ID
id_uni_sym <- mapply(function(imas, iul){
	rep(imas, iul)
}, id_map$Approved.symbol, id_uni_len, SIMPLIFY = F);
id_full_map_df <- data.frame(unlist(id_uni_sym), unlist(id_uni_ids));
colnames(id_full_map_df) <- c("Approved.symbol", "UniProt.ID.supplied.by.UniProt.");
# map uniprot ID of identified target proteins to gene symbol 
npt_top_targets <- mapply(function(nt) strsplit(nt, "_")[[1]], names(npt_table)); 
ntt_name <- sapply(npt_top_targets[2, ], function(ntt2){
	ntt2_id <- which(id_full_map_df$UniProt.ID.supplied.by.UniProt. %in% ntt2);
	ntt2_name <- id_full_map_df$Approved.symbol[[ntt2_id]];
	return(ntt2_name);
}); 
npt_df <- data.frame(npt_top_targets[2, ], ntt_name, as.integer(npt_table), as.numeric(npt_ratio), npt_top_targets[1, ], reactome_path_name[npt_top_targets[1, ]]);
colnames(npt_df) <- c("target_protein_id", "target_protein_name", "path_compound_count", "path_compound_ratio", "general_pathway_id", "general_pathway_name");
# output the identified target proteins, along with other info 
npt_file <- paste(output_folder, "_compound_new_targets_summary.tsv", sep = "");
write.table(npt_df, file = npt_file, sep = "\t", col.names = T, row.names = F, quote = F);
