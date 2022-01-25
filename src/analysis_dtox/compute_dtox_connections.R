# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script computes total number of parameters in DTox and matched fully connected multi-layer perceptron (MLP) models  


## 0. Input arguments 
dtox_optimal_file	<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv"; # name of optimal DTox model performance (on Tox21 datasets) file 
assay_info_file		<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_assay_info.tsv";	# name of Tox21 assay annotation file
onto_folder		<- "/home/yunhao1/project/ontology/data/reactome/hierarchy/";	# folder name of sorted reactome hierarchy files 
out_file		<- "data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_model_parameter_summary.tsv";	# name of output result file

## 1. Extract the optmial hyperparameter setting and input files under the setting
# read in optimal DTox model performance as data frame  
dtox_optimal_df <- read.delim(file = dtox_optimal_file, header = T, sep = "\t");
# extract the optmial hyperparameter setting and input files under the setting for each Tox21 dataset
dtox_optimal_h_files <- mapply(function(dophs){
	# Tox21 dataset name  
	dophs_s <- strsplit(dophs, "_")[[1]];
	dophs_p <- paste(dophs_s[1:6], collapse = "_");
	# pathway size (node gene number) file 
	dophs_s_file <- paste(onto_folder, dophs_p, "_st_0_node_size.tsv", sep = "");
	# parent/children node connection file 
	dophs_k_file <- paste(onto_folder, dophs_p, "_st_0_knowledge_by_node.tsv", sep = "");
	# root pathway file 
	dophs_r_file <- paste(onto_folder, dophs_p, "_st_0_root.tsv", sep = "");
	# minimal size of pathways  
	dophs_ps <- dophs_s[[4]];
	# maximal size of node modules
	dophs_xs <- dophs_s[[8]];
	return(c(dophs_s_file, dophs_k_file, dophs_r_file, dophs_ps, dophs_xs));	
}, dtox_optimal_df$hyperparameter_setting);

## 2. Compute total number of parameters in trained DTox neural network models
# iterate by Tox21 dataset 
dtox_optimal_n_connect <- mapply(function(dohf1, dohf2, dohf3, dohf4, dohf5){ 
	# read in pathway size info as data frame   
	ds_df <- read.delim(file = dohf1, header = T, sep = "\t");
	# obtain the maxial pathway size 
	dds_max <- max(ds_df$size);
	# compute the size of each node module based on pathway size   
	dohf4 <- as.integer(dohf4);
	dohf5 <- as.integer(dohf5);
	ds_df$module_size <- sapply(ds_df$size, function(dds){
		if(dds < dohf4)	return(1)
		else{
			dds_s <- 1 + (dohf5 - 1) * log(dds/dohf4)/log(dds_max/dohf4)
			return(round(dds_s))
		}
	});
	# read in parent/children node connection info as data frame  
	dk_df <- read.delim(file = dohf2, header = T, sep = "\t");	
	# compute the number of parameters (weights need to be learned) in each node module   
	dk_c_count <- mapply(function(ddn, ddcn){
		# obtain the children nodes of current node  
		ddcn_s <- strsplit(ddcn, ",")[[1]];
		ddcn_s <- as.integer(ddcn_s);
		# obtain the module sizes of children nodes  
		dds_id <- which(ds_df$node %in% ddcn_s);
		dds_ns <- ds_df$module_size[dds_id];
		# obtain the module size of current node  
		ddn_id <- which(ds_df$node %in% ddn);
		ddn_ns <- ds_df$module_size[[ddn_id]];
		# compute number of weight & bias parameters: fully connected DTox linear net module plus linear auxiliary module  
		dd_ns <- (ddn_ns * sum(dds_ns) + ddn_ns) + (ddn_ns + 1);
		return(dd_ns);
	}, dk_df$node, dk_df$children_node);
	# read in root pathway info as data frame  
	dr_df <- read.delim(file = dohf3, header = T, sep = "\t");
	# obtain the moduze sizes of root pathway nodes (root pathway nodes are connected to output by linear auxiliary module, so the size equals the number of weight parameters) 
	drd_id <- which(ds_df$node %in% dr_df$root);
	dr_c_count <- ds_df$module_size[drd_id];
	# sum up the parameters of all modules and bias of root module to get total number of parameters in DTox model 
	c_count <- sum(dk_c_count) + sum(dr_c_count) + 1;
	return(c_count);
}, dtox_optimal_h_files[1, ], dtox_optimal_h_files[2, ], dtox_optimal_h_files[3, ], dtox_optimal_h_files[4, ], dtox_optimal_h_files[5, ]);
# compute the ratio between number of training samples versus number of MLP parameters 
sample_dtox_ratio <- dtox_optimal_df$n_training/dtox_optimal_n_connect;

## 3. Compute total number of parameters in each fully connected MLP model with same number of hidden neurons as matched DTox model
# obtain the total number of hidden pathway modules in each DTox model
dtox_n_module <- sapply(dtox_optimal_df$n_hidden_module, function(dodnhm) sum(as.integer(strsplit(dodnhm, ",")[[1]])));
# iterate by Tox21 dataset  
dtox_mlp_n_connect <- mapply(function(dodnhn, dodnf){
	# obtain the number of hidden neurons in each hidden layer of current DTox model 
	dodnhn_s <- strsplit(dodnhn, ",")[[1]];
	dodnhn_s <- as.integer(dodnhn_s);
	# add number of features in input layer to list
	dodnhn_s <- c(dodnf, dodnhn_s, 1);
	ds_len <- length(dodnhn_s);	
	# compute number of weight and bias parameters between two adjacent layers: fully connected linear net module 
	ds_c_count <- sapply(1:(ds_len-1), function(dl) dodnhn_s[[dl]] * dodnhn_s[[dl+1]] + dodnhn_s[[dl+1]]);
	# sum up the parameters of all layers to get total number of parameters in MLP model 
	c_count <- sum(ds_c_count);
	return(c_count);
}, dtox_optimal_df$n_hidden_neuron, dtox_optimal_df$n_feature);
# compute the ratio between number of training samples versus number of MLP parameters
dtox_mlp_ratio <- dtox_optimal_n_connect/dtox_mlp_n_connect;

## 4. Annotate computed results with Tox21 dataset info, write to output file  
# read in Tox21 assay annotation info as data frame  
assay_info_df <- read.delim(file = assay_info_file, header = T, sep = "\t");
dtox_assay_names <- sapply(dtox_optimal_df$dataset_name, function(doddn){
	# obtain full name of current Tox21 assay dataset 
	doddn_id <- which(assay_info_df$protocol_name %in% doddn);
	doddn_name <- strsplit(assay_info_df$assay_target[[doddn_id]], " (", fixed = T)[[1]][[1]];
	# add cell line name to dataset name if it is a viability assay (to differentiate two viability assays)
	if(doddn_name == 'Cell viability')	doddn_name <- paste(doddn_name, " (", assay_info_df$cell_line[[doddn_id]], ")", sep = "");
	return(doddn_name);	
});
# aggregate all computed results in one data frame 
connect_df <- data.frame(dtox_optimal_df$dataset_name, dtox_assay_names, dtox_optimal_df$n_training, dtox_n_module, dtox_optimal_n_connect, sample_dtox_ratio, dtox_mlp_n_connect, dtox_mlp_ratio);
colnames(connect_df) <- c("dataset_name", "assay_name", "n_training", "n_pathway_module", "n_dtox_parameters", "dtox_training_parameters_ratio", "n_mlp_parameters", "dtox_mlp_parameters_ratio");
# write to output file 
write.table(connect_df, file = out_file, sep = "\t", row.names = F, col.names = T, quote = F);

