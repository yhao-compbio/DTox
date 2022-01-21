# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script uses visNetwork package to visualize the flow of relevance along DTox viability-related paths between query compound, hidden pathway modules, and the HepG2 cell viability outcome   


## functions
library(visNetwork);
library(RColorBrewer);


## 0. Input arguments 
query_pathway_file	<- "data/compound_target_probability_tox21_interpret_viability/tox21-rt-viability-hepg2_pathways.tsv";	# name of input viability-related pathway info file 
all_path_file		<- "data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_all_paths.txt";	# name of input query DTox paths file  
compound_path_file	<- "data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_path_relevance_pv.tsv";	# name of input query compound significant DTox path file 
compound_module_file	<- "data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_module_relevance.tsv";	# name of input query compound DTox module relevance score file 
module_id_file		<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/data/reactome/hierarchy/rt_25_ps_5_re_0_st_0_node.tsv";	# name of input DTox module index file  
network_node_file	<- "data/compound_target_probability_tox21_interpret_viability/tox21-rt-viability-hepg2_network_nodes.tsv";	# name of input DTox module of interest file  
query_compound		<- "CID_55245";	# name of query compound 
plot_file		<- "plot/compound_target_probability_tox21_interpret_viability/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1";	# output network plot file 

## 1. Process input viability-related pathway info and vibiality DTox model interpretation result  
# read in viability-related pathways info as data frame, obtain Reactome IDs of viability-related pathways  
query_pathway_df <- read.delim(file = query_pathway_file, header = T, sep = "\t");
query_pathways <- query_pathway_df$pathway_id;
# read in all DTox paths, iterate by path to check overlap with viability-related pathways and obtain paths of interest for plotting
all_paths <- readLines(all_path_file);
all_path_overlap <- lapply(all_paths, function(ap){
	# separate path name by underscore to obtain the pathways along the current path, check whether viability-related pathways are among them 
	ap_s <- strsplit(ap, "_")[[1]];
	ap_inter <- intersect(ap_s, query_pathways);
	if(length(ap_inter) == 0)	return(logical(0))
	else{
		# if so, build path of interst that contains the overlapped pathway and its ancestor pathways 
		ai_id <- which(ap_s %in% ap_inter)
		ai_path <- paste(ap_s[1:ai_id], collapse = "_")
		return(ai_path)
	}	
});
all_path_overlap <- unique(unlist(all_path_overlap));
# read in significant DTox paths of compounds for DTox viability assay model as data frame, select rows of interest to obtain significant DTox paths identified for the query compound  
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
cpd_id <- which(compound_path_df$cid %in% query_compound);
query_path <- compound_path_df$path_id[cpd_id];
# read in name of DTox module of interest as data frame 
network_node_df <- read.delim(file = network_node_file, header = T, sep = "\t");
# read in index of DTox module as data frame, obtain the DTox model index of viability-related pathways    
module_id_df <- read.delim(file = module_id_file, header = T, sep = "\t");
network_node_module <- sapply(network_node_df$node_id, function(nndni){
	nndni_id <- which(module_id_df$node_name %in% nndni);
	if(length(nndni_id) == 0)	return(NA)
	else	return(module_id_df$node[[nndni_id]])
});
# read in DTox module relevance scores of compounds for HepG2 viability assay model as data frame, select row/columns of interest to obtain module relevance scores of query compound with respect to viability-related pathways 
compound_module_df <- read.delim(file = compound_module_file, header = T, sep = "\t");
qc_id <- which(compound_module_df$X %in% query_compound)[[1]];
network_node_score <- sapply(network_node_module, function(nnm){
	if(is.na(nnm) == T)     nnm_score <- 0
	else{
		nnm_name <- paste("X", nnm, sep = "")
		nnm_score <- compound_module_df[qc_id, nnm_name]
		# convert negative scores to 0 (as in our interpretation method)
		if(nnm_score < 0)	nnm_score <- 0
	}
	return(nnm_score);
});

## 2. Tune network edge parameters for plotting  
# iterate by DTox paths of interest, break each path into network edges for plotting   
all_path_edges <- lapply(all_path_overlap, function(apo){
	# add root and compound nodes to the current path 
	apo_s <- strsplit(apo, "_")[[1]];
	apo_s <- c("root", apo_s, query_compound);
	# break the path into edges that connect adjacent entities 
	as_len <- length(apo_s);
	apo_edges <- mapply(function(al) apo_s[al:(al-1)], as_len:2);
	return(t(apo_edges));
});
# iterate by DTox paths of interest, check whether each path was identified with the query compound     
all_path_query <- sapply(all_path_overlap, function(apo){
	qp_apo <- sapply(query_path, function(qp) length(grep(apo, qp)));
	return(sum(qp_apo));
});
# decide whether to draw each DTox path of interst in dashed or solid lines
all_edge_dash <- mapply(function(apq, ape){
	# if path was not identified with the query compound, draw in dashed lines, otherwise draw in solid lines 
	apq_l <- (apq == 0)
	apeq <- rep(apq_l, nrow(ape));
	return(apeq);
}, all_path_query, all_path_edges);
# aggregate edges from all DTox paths of interest, obtain unique network edge set 
all_path_edges <- do.call(rbind, all_path_edges)
edges1 <- data.frame(from = all_path_edges[,1], to = all_path_edges[,2], dashes = unlist(all_edge_dash));
edges1$from_to <- mapply(function(ef, et) paste(ef, et, sep = ","), edges1$from, edges1$to);
unique_eft <- unique(edges1$from_to);
# specify whether to draw each edge in dashed or solid line 
ue_edges <- mapply(function(ue){
	# if all DTox paths associated with the edge were not identified with the query compound, draw in dashed lines, otherwise draw in solid lines 
	ue_id <- which(edges1$from_to %in% ue);
	ue_dash <- prod(edges1$dashes[ue_id]);
	return(c(ue_id[[1]], ue_dash));
}, unique_eft);
# specify edge color  
edge_col <- sapply(ue_edges[2, ], function(ue2){
	# if all DTox paths associated with the edge were not identified with the query compound, draw in grey, otherwise draw in purple 
	if(ue2 > 0)	return('grey')
	else	return('purple')
})
# store edge connection, shape, and color info in data frame 
edges <- data.frame(from = all_path_edges[ue_edges[1,], 1], to = all_path_edges[ue_edges[1,], 2], dashes = ue_edges[2, ] > 0, color = edge_col);

## 3. Tune network node parameters for plotting 
# adjust node label name to accommodate space between network layers  
network_node_name <- sapply(as.character(network_node_df$node_name), function(nndnn){
	nndnn_s <- strsplit(nndnn, " ")[[1]];
	ns_len <- length(nndnn_s);
	if(ns_len == 1)	return(nndnn)
	else{
		nndnn1 <- paste(nndnn_s[[1]], paste(nndnn_s[2:ns_len], collapse = " "), sep = "\n")
		return(nndnn1)
	}
});
network_node_name[["Cell death (HepG2)"]] <- "Cell death\n(HepG2)"
# define node color based on the scale of pathway relevance scores
nns_min <- min(network_node_score[which(network_node_score > 0)]);
network_node_score[which(network_node_score == 0)] <- nns_min/10;
network_node_score <- log10(network_node_score);
node_pal <- colorRampPalette(c('white', 'purple', 'purple1'));
network_node_col <- node_pal(10)[as.numeric(cut(network_node_score, breaks = 10))];
# define node shape based on node class: triangle for query compound, circle for pathways, square for root outcome  
node_shape <- sapply(network_node_df$node_id, function(nndni){
	if(nndni == query_compound)	return("triangle")
	else if(nndni == "root")	return("box")
	else	return("dot")
});
# define node border color, grey for pathways, black for other nodes  
network_node_border <- rep("grey", length(network_node_df$node_id));
names(network_node_border) <- network_node_df$node_id;
network_node_border[[query_compound]] <- "black";
network_node_border[["root"]] <- "black"
# store node label, color, shape, and border color info in data frame  
nodes <- data.frame(id = network_node_df$node_id, label = network_node_name, color.background = network_node_col, color.border = network_node_border, shape = node_shape);

## 3. Make the network plot based on specified parameters, showing the flow of relevance along nodes/edges in the network 
network_plot_file <- paste(plot_file, query_compound, sep = "_");
network <- visNetwork(nodes, edges,  height = "520px", width = "30%") %>%
	visNodes(size = 25, font = '25px helvetica black') %>%
	visEdges(arrows = "to") %>%
	visHierarchicalLayout(direction = "LR", sortMethod = 'directed', levelSeparation = 180) %>%
	visEvents(stabilizationIterationsDone = "function(){this.setOptions({physics:false});}") %>%
	visExport(type = "pdf", name = network_plot_file)

