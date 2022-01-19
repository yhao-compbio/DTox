# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2021
## This script uses standard Reactome pathway-receptor patterns to validate whether significant DTox paths (identified from model interpretation) contains particular pattern matched with each compound, and compare the observed outcome with expected probability 


## functions 
source("src/functions.R");


## 0. Input arguments 
Args		<- commandArgs(T);
path_folder	<- Args[1];	# folder name of input DTox model interpretation result files 
output_file	<- Args[2]; 	# name of output validation result file 
standard_file	<- "https://raw.githubusercontent.com/yhao-compbio/tox_data/master/downloads/tox21/tox21_assay_reactome_map.tsv";	# name of Tox 21 assay to Reactome pathway-receptor file 

## 1. Obtain mapping between Tox21 assay and standard Reactome pathway-receptor pattern 
# read in Tox21 assay to standard Reactome pathway-receptor pattern map as data frame   
standard_map_df <- read.delim(file = standard_file, header = T, sep = "\t");
# create a mapping vector between Tox21 assay and standard Reactome pathway-receptor pattern
standard_map <- standard_map_df$reactome_id;
names(standard_map) <- standard_map_df$assay;

## 2. Perform validation based on standard Reactome pathway-receptor pattern   
# list all DTox model interpretation result files in the folder  
path_folder_files <- list.files(path_folder);
# iterate by Tox21 assay of interst and its matched standard Reactome pathway-receptor pattern  
all_standard_valid <- lapply(names(standard_map), function(nsm){
	# obtain the DTox model interpretation result files of models built upon data from the current Tox21 assay  
	nsm_id <- sapply(path_folder_files, function(pff) length(strsplit(pff, nsm)[[1]]));
	nsm_files <- path_folder_files[nsm_id == 2];
	# stop if no result files can be found 
	if(length(nsm_files) == 0)      return(logical(0))
	# 
	else{
		# obtain the name of model interpretation result file that contains identified significant DTox paths  
		interpret_id <- sapply(nsm_files, function(nf) length(strsplit(nf, "path_relevance_pv")[[1]]))
		interpret_file <- nsm_files[[which(interpret_id == 2)]]
                interpret_file <- paste(path_folder, interpret_file, sep = "")
		# read in identified significant DTox paths of compounds as data frame 
                interpret_path_df <- read.delim(file = interpret_file, sep = "\t", header = T)
                # obtain the name of model interpretation result file that contains all possible paths in DTox neural network 
                background_id <- sapply(nsm_files, function(nf) length(strsplit(nf, "all_paths")[[1]]))
                background_file <- nsm_files[[which(background_id == 2)]]
                background_file <- paste(path_folder, background_file, sep = "")
		# read in all possible paths in DTox neural network line by line 
                background_path <- readLines(background_file)
                # compute the expected probablity and observed outcome of finding at least one standard path among the significant DTox paths  
                nsm_valid_df <- valid.path.by.standard(nsm, interpret_path_df, background_path, standard_map[[nsm]]);
                return(nsm_valid_df)
        }

});
# remove Tox21 assays without any returned alidation results 
asv_len <- sapply(all_standard_valid, length);
all_standard_valid <- all_standard_valid[asv_len > 0];
if(length(all_standard_valid) > 0){
	# extract data frame that contains the expected probablity and observed outcome of compounds, write to output file  
	all_standard_valid_df <- do.call(rbind, all_standard_valid)
	write.table(all_standard_valid_df, file = output_file, sep = "\t", row.names = F, col.names = T, quote = F)
}
