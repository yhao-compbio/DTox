# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script contains R functions required for other scripts in the repository


## This function generates all possible combinations for a list of hyperparameters 
generate.all.possible.hyperparameter.combinations <- function(hp_list, hp_sep = " "){
	## 0. Input arguments 
		# hp_list: list of hyperparameter characters  
		# hp_sep: character that seperates different hyperparameter characters

	## 1. Generate all possible hyperparameter combinations
	# iterate by hyperparameter 
	hp_current <- hp_list[[1]];
	for (i in 2:length(hp_list)){
		# identify next hyperparameter
		hp_next <- hp_list[[i]];
		# combine two hyperparmeters  
		hp_current_vec <- rep(hp_current, each = length(hp_next));
		hp_next_vec <- rep(hp_next, times = length(hp_current)); 
		hp_current <- mapply(function(hcv, hnv) paste(hcv, hnv, sep = hp_sep), hp_current_vec, hp_next_vec)
	}
	
	return(hp_current)
}


## This function generates executable shell scripts that will run input commands 
generate.parallel.bash.files <- function(all_commands, N_group, job_name, folder){
	## 0. Input arguments 
		# all_commands: a vector of commands that are to be run
		# N_group: number of groups that commands will be split into  
		# job_name: name of job 
		# folder: name of folder where shell scripts will be written 

	## 1. Assign the indices of commands to each group 
	# obtain the number of commands in each group
	N_group_member <- ceiling(length(all_commands) / N_group);
	# obtain upper bound of index for each group 
	upper_bound <- 1:N_group * N_group_member;
	ub_id <- min(which(upper_bound >= length(all_commands)));
	upper_bound <- upper_bound[1:ub_id];
	upper_bound[[ub_id]] <- length(all_commands);
	# obtain lower bound of index for each group 
	lower_bound <- 0:(ub_id-1) * N_group_member + 1;
	# assign commands to each group (lower bound - upper bound)
	command_list <- mapply(function(lb, ub) all_commands[lb:ub], lower_bound, upper_bound, SIMPLIFY = F);
	
	## 2. write commands into executable shell scripts
	# name executable shell scripts by "job_nameX.sh" (X is the index of script)
	setwd(folder);
	c_file_name <- sapply(1:ub_id, function(gn) paste(job_name, gn, ".sh", sep = ""));
	# iterate by script
	write_sub_files <- mapply(function(cl, cfn){
		# write commands into the script
		writeLines(cl, cfn);
		# make the script executable
		system(paste("chmod", "775", cfn, sep=" "));
		return(1);
	},command_list,c_file_name,SIMPLIFY=F);
	
	## 3. write an executable shell script that runs all the scripts above  
	final_command <- sapply(c_file_name,function(cfn) paste("./", folder, cfn," &", sep = ""));
	# name executable shell scripts by "job_name.sh" 
	final_file <- paste(job_name, ".sh", sep = "");
	writeLines(final_command, final_file);
	# make the script executable
	system(paste("chmod","775", final_file, sep = " "));

	return("DONE");
}


## This function groups a vector by categories of its elements.
group.vector.by.categories <- function(cate, vec){
	# 0. Input arguments 
		# cate: category of vectors  
		# vec: vector
	
	# 1. Sort vector by order of categories
	vec_od <- order(cate);
	cate <- cate[vec_od];
	vec <- vec[vec_od];

	# 2. Group elements of same category together
	# obtain unique categories
	cate_table <- table(cate);
	# obtain lower bound index/upper bound index of each unique category
	lower_ids <- cumsum(cate_table);
	upper_ids <- c(0, lower_ids[-length(lower_ids)]) + 1;
	# return list of vectors 
	vec_list <- mapply(function(li, ui) vec[li:ui], lower_ids, upper_ids, SIMPLIFY=F);
	names(vec_list) <- names(cate_table);

	return(vec_list);
}


## This function processes read-in lines from DTox performance file to extract basic info about DTox model 
process.dtox.config.lines <- function(file_lines){
	## 0. Input arguments
		# file_lines: list of character strings that contain lines read from DTox performance file 
	
	## 1. extract basic info of DTox visual neural network model
	# extract number of hidden layers (line 1), number of features (line 2) 
	vnn_cf1 <- sapply(file_lines[1:2], function(fl) as.integer(strsplit(fl, ":")[[1]][[2]]));
	# extract number of hidden modules in layer (line 3, character format, numbers separated by ','), extract number of hidden modules in layer (line 3, character format, numbers separated by ',')
	vnn_cf2 <- sapply(file_lines[3:4], function(fl) strsplit(fl, ": ")[[1]][[2]]);
	# extract number of training samples (line 5), number of validation samples (lines 7) 
	vnn_cf3 <- sapply(file_lines[c(5,7)], function(fl) as.integer(strsplit(fl, ":")[[1]][[2]]));
	        
	## 2. Store extracted info in data frame form
	config_df <- data.frame(t(vnn_cf1), t(vnn_cf2), t(vnn_cf3));
	colnames(config_df) <- c("n_hidden_layer", "n_feature", "n_hidden_module", "n_hidden_neuron", "n_training", "n_testing");

	return(config_df);
}


## This function processes read-in lines from DTox performance file to extract values of model performance metrics 
process.metric.lines <- function(file_lines, perf_line, sep_char){
	## 0. Input arguments
		# file_lines: list of character strings that contain lines read from DTox performance file 
		# perf_line: line number that contains model performance metrics
		# sep_char: character to be recognized that comes before model performance metrics     
		
	## 1.
	# extract character of performance metrics from the line 
	perf_row <- strsplit(file_lines[[perf_line]], sep_char)[[1]][[2]]
	# separate different performance metrics, extract values of each metric 
	perf_metric_str <- strsplit(perf_row, ",")[[1]]
	perf_metric_mat <- mapply(function(pms) strsplit(pms, ":")[[1]], perf_metric_str) 
	perf_op <- data.frame(t(as.numeric(perf_metric_mat[2,])))
	colnames(perf_op) <- perf_metric_mat[1,]

	return(perf_op);
};


## This function reads in DTox performance file then extracts model basic info and performance metrics  
read.dtox.performance.files <- function(file_name){
	## 0. Input arguments
		# file_name: name of DTox performance file 

	## 1. Read in performance files line by line as character 
	f_lines <- readLines(file_name);

	## 2. Process read-in lines to extract model basic info and performance metrics 
	# extract basic info about DTox model  
	vnn_config <- process.dtox.config.lines(f_lines);
	# extract training performance of DTox model
	training_perf <- process.metric.lines(f_lines, 6, "performance: ");
	colnames(training_perf) <- sapply(colnames(training_perf), function(ctp) paste("training", ctp, sep = "_"));
	# extract validation performance of DTox model
	testing_perf <- process.metric.lines(f_lines, 8, "performance: ");
	colnames(testing_perf) <- sapply(colnames(testing_perf), function(ctp) paste("testing", ctp, sep = "_"));
	# combine extracted info in data frame form	
	f_df <- cbind(vnn_config, training_perf, testing_perf);

	return(f_df); 	
}


## This function reads in simple machine learning performance file then extracts model basic info and performance metrics   
read.simple.performance.files <- function(file_name){
	## 0. Input arguments 
		# file_name: name of simple machine learning performance file
	
	## 1. Read in performance files line by line as character  
	f_lines <- readLines(file_name);

	## 2. Process read-in lines to extract model basic info and performance metrics 
	# extract number of training samples (line 1) and number of validation samples (lines 2) 
	n_samples <- sapply(f_lines[1:2], function(fl) as.integer(strsplit(fl, "instances: ")[[1]][[2]]));
	n_samples_df <- data.frame(t(n_samples));
	colnames(n_samples_df) <- c("n_training", "n_testing");
	# extract training performance of simple machine learning model
	training_perf <- process.metric.lines(f_lines, 4, "performance: ");
	colnames(training_perf) <- sapply(colnames(training_perf), function(ctp) paste("training", ctp, sep = "_"));
	# extract validation performance of simple machine learning model
	testing_perf <- process.metric.lines(f_lines, 5, "performance: ");
	colnames(testing_perf) <- sapply(colnames(testing_perf), function(ctp) paste("testing", ctp, sep = "_"));
	# combine extracted info in data frame form   
	f_df <- cbind(n_samples_df, training_perf, testing_perf);
	
	return(f_df);   
}


## This function reads in multi-layer perceptron (MLP) model performance file then extracts model basic info and performance metrics  
read.mlp.performance.files <- function(file_name){
	## 0. Input arguments
		# file_name: name of MLP model performance file
	
	## 1. Read in performance files line by line as character   
	f_lines <- readLines(file_name);
	
	## 2. Process read-in lines to extract model basic info and performance metrics 
	# extract number of hidden layers (line 1), number of features (line 2)  
	mlp_cf1 <- sapply(f_lines[1:2], function(fl) as.integer(strsplit(fl, ": ")[[1]][[2]]));
	# extract number of hidden modules in layer (line 3, character format, numbers separated by ',')
	mlp_cf2 <- sapply(f_lines[[3]], function(fl) strsplit(fl, ": ")[[1]][[2]]);
	# extract number of training samples (line 4), number of validation samples (lines 6)  
	mlp_cf3 <- sapply(f_lines[c(4,6)], function(fl) as.integer(strsplit(fl, ": ")[[1]][[2]]));
	config_df <- data.frame(t(mlp_cf1), t(mlp_cf2), t(mlp_cf3));
	colnames(config_df) <- c("n_hidden_layer", "n_feature", "n_hidden_neuron", "n_training", "n_testing");
	# extract training performance of MLP model
	training_perf <- process.metric.lines(f_lines, 5, "performance: ");
	colnames(training_perf) <- sapply(colnames(training_perf), function(ctp) paste("training", ctp, sep = "_"));
	# extract validation performance of MLP model
	testing_perf <- process.metric.lines(f_lines, 7, "performance: ");
	colnames(testing_perf) <- sapply(colnames(testing_perf), function(ctp) paste("testing", ctp, sep = "_"));
	# combine extracted info in data frame form 
	f_df <- cbind(config_df, training_perf, testing_perf);

	return(f_df);
}


## This function performs Fisher Exact Test between a set of interesting genes and a list of gene sets.
compute.fisher.p.value <- function(geneSetList, interestingGene, background){
	## 0. Input arguments: 
		# geneSetList: list of gene sets
		# interestingGene: interesting genes 
		# background: background genes

	## 1. Remove genes that are not in the background list from the ones to be tested    
	interestingGene <- intersect(interestingGene,background);
	geneSetList <- lapply(geneSetList,function(x) intersect(x,background));
	
	## 2. Perform Fisher Exact Test for each gene set  
	# iterate by gene set 
	fisher_pval <- sapply(geneSetList,function(gs){
		# count numbers in the 2*2 contigency table  
		GsIn <- length(intersect(gs,interestingGene));
		diffin <- setdiff(background,interestingGene);
		diffgs <- setdiff(background,gs);
		DiffgsIn <- length(intersect(interestingGene,diffgs));
		DiffinGs <- length(intersect(gs,diffin));
		DiffinDiffgs <- length(intersect(diffin,diffgs));
		# compute P value of Fisher Exact Test  
		fisher.test(matrix(c(GsIn,DiffgsIn,DiffinGs,DiffinDiffgs),2,2),alternative="greater")$p.value;
	});
	names(fisher_pval) <- names(geneSetList);
	return(fisher_pval);
}


## This function identifies differentially expressed DTox paths after compound treatment from LINCS data (output 2), then compares proportion of differential expression to background DTox paths (output 1) 
valid.path.by.expression <- function(query_assay, ip_path_df, bg_path, gene_map_df, pathway_gene_list, pathway_fdr_cut, exp_map_df, gene_exp_df){
	## 0. Input argument
		# query_assay: name of query Tox21 assay
		# ip_path_df: data frame that contains identified significant DTox paths of compounds 
		# bg_path: vector that contains all possible paths in DTox neural 
		# gene_map_df: data frame that contains uniprot to entrez ID map of genes measured by LINCS  
		# pathway_gene_list: list that contains annotated genes (entrez ID) of pathways 
		# pathway_fdr_cut: number that indicates FDR threshold for differentially expressed pathway 
		# exp_map_df: data frame that contains processed LINCS instance-profile map  
		# gene_exp_df: data frame that contains processed LINCS gene expression data  

	## 1. Map pathways along DTox pathways to their gene annotations     
	# seperate string by underscore to obtain gene/pathways along each DTox path  
	bg_path_entity <- lapply(bg_path, function(bp) strsplit(bp, "_")[[1]]);
	# check whether pathways along each DTox pathway are included in the pathway annotation list    
	bg_path_entity <- lapply(bg_path_entity, function(bpe){
		# obtain the pathways along current DTox path
		N_bpe <- length(bpe);
		bpe_path_id <- which(names(pathway_gene_list) %in% bpe[1:(N_bpe-1)]);
		# check if the obtained pathways are included in the pathway annotation list 
		if(length(bpe_path_id) == (N_bpe-1)) return(bpe[1:(N_bpe-1)])
		else	return(logical(0))
	});
	# remove DTox paths with pathways that are not included in the pathway annotation list  
	bpe_len <- sapply(bg_path_entity, length);
	bg_path_entity <- bg_path_entity[bpe_len > 0];
	bg_path_entity <- unique(bg_path_entity);
	# name each DTox pathways after pathways along the path (separated by underscore) 
	names(bg_path_entity) <- sapply(bg_path_entity, function(bpe) paste(bpe, collapse = "_"));
	# obtain names of all pathways that appear along DTox paths, as well as their gene annotations  
	all_entity <- unique(unlist(bg_path_entity));
	pathway_gene_list <- pathway_gene_list[all_entity];

	## 2. Identify differentially expressed DTox paths from all compound-path pairs, compare proportion of differential expression to background DTox paths  
	# select rows from processed LINCS instance-profile map that are relevant to the query Tox21 assay
	assay_id <- which(exp_map_df$assay_name %in% query_assay);
	exp_map_df <- exp_map_df[assay_id, ];
	# create a mapping vector bewteen column index of express data and compound CID  
	exp_map <- exp_map_df$col_id;
	names(exp_map) <- exp_map_df$compound_pubchem_cid;
	# check whether identified compounds from interpretation results are in LINCS expression data, only process when such compound exists    
	nem_id <- which(ip_path_df$cid %in% names(exp_map));
	if(length(nem_id) == 0)	return(logical(0))
	else{
		# group identified significant DTox paths by compounds  
		ip_path_list <- group.vector.by.categories(ip_path_df$cid[nem_id], ip_path_df$path_id[nem_id])
		# iterate by compound, obtain the pathways along each identified significant DTox path 
		ip_path_list <- lapply(ip_path_list, function(ipl){
			ipl_sp <- sapply(ipl, function(ip){
				ip_s <- strsplit(ip, "_")[[1]]
				N_is <- length(ip_s)
				ip_p <- paste(ip_s[1:(N_is-1)], collapse = "_")
			})
			ipl_inter <- intersect(unique(ipl_sp), names(bg_path_entity))
			return(ipl_inter)
		});
		ipl_len <- sapply(ip_path_list, length)
		ip_path_list <- ip_path_list[ipl_len > 0]
		# obtain all genes measured by LINCS (used as background list in Fisher's exact test) 
		bg_genes <- rownames(gene_exp_df)
		# filter compound-column ID map vector and gene expression data frame (column) by identified compounds, keep gene expression in data frame form in case of one identified compound-instance 
		exp_map <- exp_map[names(ip_path_list)]
		gene_exp_df <- gene_exp_df[, exp_map]
		if(length(exp_map) == 1){
			gene_exp_df <- data.frame(gene_exp_df)
			rownames(gene_exp_df) <- bg_genes 
		}	
		# iterate by instance (column of expression data frame), identify the differentially expressed genes of each instance (score > 2)  
		sig_genes <- lapply(1:ncol(gene_exp_df), function(nged){
			nged_sig_id <- which(abs(gene_exp_df[, nged]) > 2)
			nged_sig_genes <- rownames(gene_exp_df)[nged_sig_id]
			return(nged_sig_genes)
		})
		# perform pathway enrichment analysis using the identified differentially expressed genes of each instance    
		pathway_enrich <- mapply(function(sg){
			# use Fisher's exact test to compute pathway enrichment P-values for current instance  
			sg_enrich_pv <- compute.fisher.p.value(pathway_gene_list, sg, bg_genes)
			return(sg_enrich_pv)
		}, sig_genes, SIMPLIFY = TRUE)
		# iterate by DTox path, check whether pathways along each DTox path are differentially expressed in each compound-instance, keep result in matrix form in case of one identified compound-instance 
		bg_path_valid <- mapply(function(bpe){
			# obtain computed enrichment P-values of pathways along current DTox path, keep result in matrix form in case of one compound-instance  
			bpe_path_pv <- pathway_enrich[bpe, ]
			if(length(exp_map) == 1){
				bpe_path_pv <- matrix(bpe_path_pv, length(bpe_path_pv))
				rownames(bpe_path_pv) <- bpe
			}
			# perform multile testing adjustment by FDR   
			bpe_path_fdr <- apply(bpe_path_pv, 2, function(bpp) p.adjust(bpp, method = "fdr"))
			# check whether all pathways along current DTox path are differentially expressed in each compound-instance    
			bpe_path_sig <- as.integer(colSums(bpe_path_fdr < pathway_fdr_cut) == length(bpe))
			return(bpe_path_sig)
		}, bg_path_entity, SIMPLIFY = TRUE)
		if(length(exp_map) == 1){
			bg_path_valid <- matrix(bg_path_valid, 1, length(bg_path_valid))
			colnames(bg_path_valid) <- names(bg_path_entity)
		}
		bg_path_valid <- t(bg_path_valid);
		# iterate by identified compound, use mapped instance column index to select the validation result of identified significant DTox paths from all result of DTox paths (from last step)    
		ip_path_valid <- mapply(function(ipl, lipl){
			ipl_valid <- bg_path_valid[ipl, lipl]
			return(ipl_valid)
		}, ip_path_list, 1:length(ip_path_list), SIMPLIFY = FALSE)
	 	# compute the proportion of DTox paths/identified significant DTox paths that are differentially expressed in each compound-instance, store results in data frame form along with compound and assay info 
		bg_path_sig_ratio <- apply(bg_path_valid, 2, mean)
		ip_path_sig_ratio <- sapply(ip_path_valid, mean) 
		valid_df <- data.frame(rep(query_assay, length(exp_map)), names(exp_map), ip_path_sig_ratio, bg_path_sig_ratio)
		colnames(valid_df) <- c("assay", "compound", "interpret_sig_ratio", "background_sig_ratio")
		# obtain the identified significant DTox paths of each compound that are also differentially expressed, store results in data frame form along with compound and assay info 
		ip_sig_path <- lapply(ip_path_valid, function(ipv) names(ipv)[which(ipv > 0)])
		ip_sig_len <- sapply(ip_sig_path, length)
		ip_sig_compounds <- mapply(function(nisp, isl) rep(nisp, isl), names(ip_sig_path), ip_sig_len)
		ip_sig_paths <- unlist(ip_sig_path)
		ip_sig_path_df <- data.frame(rep(query_assay, length(ip_sig_paths)), rep(length(ip_sig_len), length(ip_sig_paths)), unlist(ip_sig_compounds), ip_sig_paths)
		colnames(ip_sig_path_df) <- c("assay", "n_tested_compound", "compound", "sig_path")
		return(ls = list(valid = valid_df, sig_path = ip_sig_path_df))
	}
}


## This function evaluates the expected probablity and observed outcome of finding at least one standard path among the significant DTox paths
valid.path.by.standard <- function(query_assay, ip_path_df, bg_path, query_standard){
	## 0. Input argument
		# query_assay: name of query Tox21 assay 
		# ip_path_df: data frame that contains identified significant DTox paths of compounds 
		# bg_path: vector that contains all possible paths in DTox neural  
		# query_standard: character that contains query standard pathway-receptor pattern

	## 1. Identify DTox paths that match the standard pathway-receptor pattern 
	query_standard_id <- sapply(bg_path, function(bp) length(grep(query_standard, bp)));
	standard_path <- bg_path[query_standard_id == 1]; 
	N_sp <- length(standard_path)

	## 2. Evaluate the expected probablity and observed outcome for each compound
	# group identified significant DTox paths by compounds   
	ip_path_list <- group.vector.by.categories(ip_path_df$cid, ip_path_df$path_id);
	ipl_len <- sapply(ip_path_list, length);
	# iterate by compound, compute the expected probablity of finding at least one standard path among the significant DTox paths  
	N_bp <- length(bg_path);
	ipl_expected <- sapply(ipl_len, function(il){
		il1 <- (N_bp - il):(N_bp - il - N_sp + 1)
		il2 <- N_bp:(N_bp - N_sp + 1)
		il_prob <- 1 - prod(il1)/prod(il2);
		return(il_prob);
	});
	# iterate by compound, obtain the observed outcome of finding at least one standard path among the significant DTox paths  
	ipl_observed <- sapply(ip_path_list, function(ipl){
		ipl_id <- which(ipl %in% standard_path)
		if(length(ipl_id) > 0)	return(1)
		else	return(0)
	});
	# store expected probability and observed outcome in data frame, along with assay and compound info 
	valid_df <- data.frame(rep(query_assay, length(ipl_len)), names(ip_path_list), ipl_observed, ipl_expected);
	colnames(valid_df) <- c("assay", "compound", "whether_contain_standard", "bg_standard_ratio");		
 
	return(valid_df);
}
