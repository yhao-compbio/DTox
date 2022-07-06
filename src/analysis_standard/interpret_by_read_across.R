# !usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2022
## This script implements Read-across to connect query compounds with query target based on their chemical similarity to source compounds in Drugbank/ComptoxAI    

## functions 
library(parallel);
source("src/functions.R");


## 0. Input arguments
Args		<- commandArgs(T); 
data_file	<- Args[1];	# name of compound MACCS fingerprint-Tox21 assay outcome dataset file
query_assay	<- Args[2];	# name of query Tox21 assay 
query_target	<- Args[3];	# Uniprot ID of query target protein 
query_gene_id	<- Args[4];	# Entrez (NCBI) ID of query target gene  
output_folder	<- Args[5];	# name of output file 
# 
dtox_file	<- "data/compound_target_probability_tox21_interpret_standard/validation_summary/compound_target_fingerprint_maccs_probability_tox21_interpret_standard_validation_all_result.tsv";
db_target_file	<- "data/compound_target_probability_tox21_interpret_standard/read_across_interpret_result/uniprot\ links.csv";
cai_target_file	<- "/home/yunhao1/project/comptox_ai/data/tox21_connection/tox21_10k_library_comptoxai_gene_connections.tsv"	
compound_file	<- "/home/yunhao1/project/tox_data/downloads/tox21/tox21_10k_library_info.tsv";
sim_cut		<- 16:20/20;

## 1. Process input datasets  
# read in the validation results of DTox, obtain the validated compounds of query assay  
dtox_df <- read.delim(file = dtox_file, sep = "\t", header = T);
dd_id1 <- which(dtox_df$rule %in% "gamma-epsilon_0.001_0.1");
dtox_df <- dtox_df[dd_id1, ];
dd_id2 <- which(dtox_df$assay %in% query_assay);
dtox_df <- dtox_df[dd_id2, ];
# read in compound MACCS fingerprint-Tox21 assay outcome dataset, obtain the MACCS fingerprint of DTox-validated compounds  
data_df <- read.delim(file = data_file, sep = "\t", header = T);
fps <- setdiff(colnames(data_df), c("X", "assay_outcome"));
fp_df <- data_df[, fps];
ddc_id <- sapply(dtox_df$compound, function(ddc) which(data_df$X %in% ddc)[[1]]);
ddc_fp_df <- fp_df[ddc_id, ];
rownames(ddc_fp_df) <- dtox_df$compound;
# read in Tox21 compound info as data frame, convert compound names to lower case, creat map data frame between compound name PubChem ID 
compound_df <- read.delim(file = compound_file, header = T, sep = "\t");
compound_df$cid <- sapply(compound_df$PUBCHEM_CID, function(tdpc) paste("CID", tdpc, sep = "_"));
compound_df$compound_name <- tolower(compound_df$SAMPLE_NAME);
compound_df <- compound_df[, c("cid", "compound_name")];
compound_df <- unique(compound_df);
cd_id <- which(compound_df$cid %in% dtox_df$compound);
compound_df <- compound_df[cd_id, ];

## 2. Perform Read-across interpretation using compound-target relationships from DrugBank 
# read in compound-target relationships from DrugBank, extract connections that involve query target protein
db_target_df <- read.csv(file = db_target_file);
dtd_id <- which(db_target_df$UniProt.ID %in% query_target);
db_target_df <- db_target_df[dtd_id, ];
db_target_df$Name <- tolower(db_target_df$Name);
# obtain the DTox-validated compounds with known query target connection from DrugBank (source compounds)
db_compound_target_df <- merge(compound_df, db_target_df, by.x = "compound_name", by.y = "Name");
db_reference_compound <- unique(db_compound_target_df$cid);
# Implement Read-across, assigning targets of source compounds to all query compounds based on specified similarty threshold  
if(length(db_reference_compound) == 0){
	db_sim_compound_ratio <- rep(0, length(sim_cut))
}
if(length(db_reference_compound) > 0){ 
	# compute pairwise Jaccard similarity between source compounds and all the query DTox-validated compounds
	db_reference_compound_sim <- mapply(function(drc){
		drc_bid <- which(ddc_fp_df[drc, ] == 1) 
		drc_sim <- sapply(rownames(ddc_fp_df), function(rdfd){
			rdfd_bid <- which(ddc_fp_df[rdfd, ] == 1)
			rdfd_inter <- intersect(drc_bid, rdfd_bid)
			rdfd_union <- union(drc_bid, rdfd_bid)
			rdfd_sim <- length(rdfd_inter)/length(rdfd_union)
			return(rdfd_sim)
		});
		return(drc_sim)
	}, db_reference_compound)
	db_reference_compound_sim <- round(db_reference_compound_sim, 3)
	# identify the query DTox-validated compounds with higher Jaccard similarity than the specified threshold 
	db_sim_compound_count <- sapply(sim_cut, function(sc){
		sc_id_list <- lapply(colnames(db_reference_compound_sim), function(cdrcs) which(db_reference_compound_sim[, cdrcs] >= sc))	
		sc_ids <- unique(unlist(sc_id_list))
		return(length(sc_ids))
	})
	# compute the proportion of read-across-validated compounds
	db_sim_compound_ratio <- db_sim_compound_count/nrow(ddc_fp_df)
}

## 3. Perform Read-across interpretation using compound-target relationships from ComptoxAI
# read in compound-target relationships from ComptoxAI, extract connections that involve query target protein
cai_target_df <- read.delim(file = cai_target_file, sep = "\t", header = T);
ctd_id1 <- which(cai_target_df$gene_id %in% as.integer(query_gene_id));
cai_target_df <- cai_target_df[ctd_id1, ];
# obtain the DTox-validated compounds with known query target connection from ComptoxAI (source compounds)
cai_compound_target_df <- merge(compound_df, cai_target_df, by.x = "cid", by.y = "compound_cid");
cai_reference_compound <- unique(cai_compound_target_df$cid);
# Implement Read-across, assigning targets of source compounds to all query compounds based on specified similarty threshold   
if(length(cai_reference_compound) == 0){
	cai_sim_compound_ratio <- rep(0, length(sim_cut))
}
if(length(cai_reference_compound) > 0){
	# compute pairwise Jaccard similarity between source compounds and all the query DTox-validated compounds
	cai_reference_compound_sim <- mapply(function(crc){
		crc_bid <- which(ddc_fp_df[crc, ] == 1)
		crc_sim <- sapply(rownames(ddc_fp_df), function(rdfd){
			rdfd_bid <- which(ddc_fp_df[rdfd, ] == 1)
			rdfd_inter <- intersect(crc_bid, rdfd_bid)
			rdfd_union <- union(crc_bid, rdfd_bid)
			rdfd_sim <- length(rdfd_inter)/length(rdfd_union)
			return(rdfd_sim)
		});
		return(crc_sim)
	}, cai_reference_compound)
	cai_reference_compound_sim <- round(cai_reference_compound_sim, 3)
	# identify the query DTox-validated compounds with higher Jaccard similarity than the specified threshold  
	cai_sim_compound_count <- sapply(sim_cut, function(sc){
		sc_id_list <- lapply(colnames(cai_reference_compound_sim), function(ccrcs) which(cai_reference_compound_sim[, ccrcs] >= sc))
		sc_ids <- unique(unlist(sc_id_list))
		return(length(sc_ids))
	})
	# compute the proportion of read-across-validated compounds
	cai_sim_compound_ratio <- cai_sim_compound_count/nrow(ddc_fp_df)
}

## 4. Output the Read-across validation results using compound-target relationships from DrugBank and ComptoxAI
sim_result_df <- data.frame(query_assay, sim_cut, db_sim_compound_ratio, cai_sim_compound_ratio);
colnames(sim_result_df) <- c("assay", "tc_threshold", "standard_observed_ratio_drugbank", "standard_observed_ratio_comptoxai");
sim_result_file <- paste(output_folder, query_assay, "interpret_by_read_across_result.tsv", sep = "_");
write.table(sim_result_df, file = sim_result_file, sep = "\t", col.names = T, row.names = T, quote = F);
