# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script uses LINCS pertubation gene expression data to validate whether significant DTox paths (identified from model interpretation) are differentially expressed after compound treatment, and compare the proportion of differential expression to backtround DTox paths 


## functions  
source("src/functions.R");


## 0. Input arguments  
Args			<- commandArgs(T);
path_folder		<- Args[1];		# folder name of input DTox model interpretation result files 
lincs_exp_file		<- Args[2];		# name of processed LINCS gene expression profile file  
lincs_map_file		<- Args[3];		# name of processed LINCS instance-profile mapping file 
fdr_cut			<- as.numeric(Args[4]);	# number that indicates FDR threshold for differentially expressed pathway 
output_file		<- Args[5];		# folder name of output validation result files
uniprot2entrez_file	<- "/home/yunhao1/project/target/downloads/id_map/HUMAN_9606_idmapping_selected.tab";	# name of uniprot to entrez gene ID mapping file   
reactome_annot_file	<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";	# name of Reactome protein-pathway annotation file 
lincs_gene_file		<- "/home/yunhao1/project/tox_data/downloads/LINCS/geneinfo_beta.txt";	# name of LINCS gene ID file    
 
## 1. Process gene ID map data 
# read in uniprot-entrez ID mapping info as data frame  
uniprot2entrez_df <- read.delim(file = uniprot2entrez_file, sep = "\t", header = F)
# build unique mapping data frame between uniprot protein ID and 
valid_id <- which(uniprot2entrez_df$V3 != "");
uniprot2entrez_df <- unique(uniprot2entrez_df[valid_id, c("V1", "V3")]);
# read in info of genes measured by LINCS as data drame  
lincs_gene_df <- read.delim(file = lincs_gene_file, sep = "\t", header = T); 
lincs_genes <- as.character(lincs_gene_df$gene_id);
# create a gene ID map data frame for genes measured by LINCS  
lg_id <- which(uniprot2entrez_df$V3 %in% lincs_genes);
lincs_gene_map <- uniprot2entrez_df[lg_id, ];

## 2. Process Reactome annotation data  
# read in Reactome protein-pathway annotation as data frame 
reactome_df <- read.delim(file = reactome_annot_file, header = F, sep = "\t");
# remove non-human annotation rows in the data frame 
human_id <- which(reactome_df[, 6] %in% "Homo sapiens");
eactome_df <- reactome_df[human_id, c(1, 2)];
reactome_df$V1 <- sapply(reactome_df$V1, function(rdv) strsplit(rdv, "-")[[1]][[1]]);
reactome_df <- unique(reactome_df);
# group protein annotations by pathway   
reactome_uni_list <- group.vector.by.categories(reactome_df[, 2], reactome_df[, 1]);
# iterate by pathway, map unirot protein IDs in each pathway to gene entrez IDs 
reactome_gene_list <- lapply(reactome_uni_list, function(rul){
	# obtain the unique gene entrez IDs of current pathway   
	rul_id <- which(lincs_gene_map$V1 %in% rul);
	rul_gene <- unique(lincs_gene_map$V3[rul_id]);
	return(rul_gene);
});
# remove pathways with no mapped gene entrez ID 
rgl_len <- sapply(reactome_gene_list, length);
reactome_gene_list <- reactome_gene_list[rgl_len > 0];

## 3. Perform validation using LINCS pertubation gene expression data  
# list all DTox model interpretation result files in the folder
path_folder_files <- list.files(path_folder);
# read in processed LINCS gene expression data as data frame  
lincs_exp_df <- read.delim(file = lincs_exp_file, header = F, sep = "\t", row.names = 1);
# read in processed LINCS instance-profile map as data frame, obtain unique names of mapped Tox21 assays   
lincs_map_df <- read.delim(file = lincs_map_file, header = T, sep = "\t");
all_map_assays <- unique(lincs_map_df$assay_name);
# iterate by Tox21 assay  
all_map_assay_valid <- lapply(all_map_assays, function(ama){
	# obtain the DTox model interpretation result files of models built upon data from the current Tox21 assay 
	ama_id <- sapply(path_folder_files, function(pff) length(strsplit(pff, ama)[[1]]));
	ama_files <- path_folder_files[ama_id == 2];
	# stop if no result files can be found 
	if(length(ama_files) == 0)	return(logical(0))
	# 
	else{
		# obtain the name of model interpretation result file that contains identified significant DTox paths  
		interpret_id <- sapply(ama_files, function(af) length(strsplit(af, "path_relevance_pv")[[1]]))
		interpret_file <- ama_files[[which(interpret_id == 2)]]
		interpret_file <- paste(path_folder, interpret_file, sep = "")
		# read in identified significant DTox paths of compounds as data frame 
		interpret_path_df <- read.delim(file = interpret_file, sep = "\t", header = T)
		# obtain the name of model interpretation result file that contains all possible paths in DTox neural network 
		background_id <- sapply(ama_files, function(af) length(strsplit(af, "all_paths")[[1]]))
		background_file <- ama_files[[which(background_id == 2)]]
		background_file <- paste(path_folder, background_file, sep = "")
		# read in all possible paths in DTox neural network line by line 
		background_path <- readLines(background_file)
		# identify differentially expressed DTox paths after compound treatment from LINCS data (validation result 2), compare proportion of differential expression to background DTox paths (validation result 1) 
		ama_valid_df <- valid.path.by.expression(ama, interpret_path_df, background_path, lincs_gene_map, reactome_gene_list, fdr_cut, lincs_map_df, lincs_exp_df)
		return(ama_valid_df)
	}
});
# remove Tox21 assays without any returned alidation results  
amav_len <- sapply(all_map_assay_valid, length);
all_map_assay_valid <- all_map_assay_valid[amav_len > 0];
if(length(all_map_assay_valid) > 0){
	# extract data frame that contains compared proportion of differential expression between significant DTox paths and background DTox paths, write to output file  
	all_map_assay_valid_result <- lapply(all_map_assay_valid, function(amav) amav[[1]])
	all_map_assay_valid_result_df <- do.call(rbind, all_map_assay_valid_result)
	all_map_assay_valid_result_file <- paste(output_file, "_result.tsv", sep = "")
	write.table(all_map_assay_valid_result_df, file = all_map_assay_valid_result_file, sep = "\t", row.names = F, col.names = T, quote = F)
	# extract data frame that contains identified differentially expressed DTox paths after compound treatment from LINCS data, write to output file
	all_map_assay_valid_sig <- lapply(all_map_assay_valid, function(amav) amav[[2]])
	all_map_assay_valid_sig_df <- do.call(rbind, all_map_assay_valid_sig)	
 	all_map_assay_valid_sig_file <- paste(output_file, "_sig_path.tsv", sep = "")
	write.table(all_map_assay_valid_sig_df, file = all_map_assay_valid_sig_file, sep = "\t", row.names = F, col.names = T, quote = F)
}
