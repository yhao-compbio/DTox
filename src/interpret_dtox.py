# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains implements layer-wise relevance propagation to evaluate relevance of DTox paths


## Module
import sys
import numpy as np
import pandas as pd
import torch
sys.path.insert(0, 'src/')
import dtox_hierarchy
import dtox_nn
import dtox_data
import dtox_lrp


## Main function
def main(argv):
	## 0. Input arguments 
		# argv 1: input file that contains root pathways of DTox model 
		# argv 2: input file that contains parent/children node connections of DTox model  
		# argv 3: input file that contains number of annotated genes in nodes of DTox model
		# argv 4: input file that contains layer number of nodes in DTox model 
		# argv 5: minimal size of pathways included in DTox model
		# argv 6: maximal size of node modules in DTox neural network model
		# argv 7: input file that contains structure/gene/pathway name of nodes in DTox model
		# argv 8: input file that contains trained DTox model 
		# argv 9: input file that contains names of trained null DTox model file 
		# argv 10: input file that contains combined query feature-response data 
		# argv 11: name of column that contains label/response data
		# argv 12: propagation rule to be implemented ('gamma-epsilon' for generic rule, 'alpha-beta' for αβ rule)  
		# argv 13: first factor in LRP rules ('gamma' in generic rule, 'alpha' in αβ rule)  
		# argv 14: second factor in LRP rules ('epsilon' in generic rule, 'beta' in αβ rule)  
		# argv 15: name of output file 

	## 1. Load trained DTox model 
	# sort the hierarchy of DTox model based on input hierarchy files
	dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_hierarchy_stat, dtox_input_size, dtox_hidden_size = dtox_hierarchy.sort_dtox_hierarchy(argv[1], argv[2], argv[3], argv[4], int(argv[5]), int(argv[6])) 
	# define structure of whole neural network 
	dtox_model = dtox_nn.DToxNet(dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_input_size) 
	# load trained DTox model parameter into defined structure
	model_state = torch.load(argv[8])
	dtox_model.load_state_dict(model_state['model_state_dict'])

	## 2. Format query dataset  
	# read in query dataset, remove samples with negative label (interpretation only focuses on positive outcome) 
	query_data_df = pd.read_csv(argv[10], sep = '\t', header = 0, index_col = 0)
	label_col = argv[11]
	query_data_df = query_data_df[query_data_df[label_col] == 1]
	# format featuer and label data of query dataset 
	X_query, y_query = query_data_df.drop(label_col, axis = 1).values, query_data_df[label_col].values
	query_data = dtox_data.DTox_dataformat(X_query, y_query)
	query_feature_data = query_data[: , :][0]

	## 3. Implement layer-wise relevance propagation (LRP) process to compute the relevance of all neurons/modules in learned DTox model 
	hidden_relevance_df, pathway_relevance_df = dtox_lrp.lrp(dtox_model, query_feature_data, argv[12], float(argv[13]), float(argv[14]), 0, 1)
	# write data frame of neuron relevance scores of compounds to output file  
	hidden_relevance_df.index = query_data_df.index
	hidden_relevance_df.to_csv(argv[15] + '_neuron_relevance.tsv', sep = '\t')
	# write data frame of node module relevance scores of compounds to output file 
	pathway_relevance_df.index = query_data_df.index
	pathway_relevance_df.to_csv(argv[15] + '_module_relevance.tsv', sep = '\t')
	
	## 4. Compute observed relevance scores for all DTox paths connecting input layer to root module  
	path_relevance_df = dtox_lrp.compute_path_relevance_score(dtox_model, pathway_relevance_df)
	# read in data frame that contains node name info 
	node_map_df = pd.read_csv(argv[7], sep = '\t', header = 0)
	# iterate by path, name each path by names of nodes in the path  	
	path_name_dict = {}
	for prdc in path_relevance_df.columns:
		# link node numbers/names of current path by underscores respectively, store in dictionary 
		prdc_node_id = [int(pr) for pr in prdc.split('_')]
		prdc_path_name = '_'.join(node_map_df.node_name[prdc_node_id].values)
		path_name_dict[prdc] = prdc_path_name
	# write all named paths to output txt file  
	path_file = open(argv[15] + '_all_paths.txt', 'w')
	for lpv in list(path_name_dict.values()):
		path_file.write('%s\n' % lpv)
	path_file.close()
	
	## 5. Compute empirical p-values of module relevance scores for all DTox path-compound pairs 
	# read in names of null DTox model files 
	dtox_null_files = np.loadtxt(argv[9], dtype = 'str')
	# iterate by trained null DTox model
	null_relevance_list = []
	for dnf in dtox_null_files:
		# define structure of whole neural network of null model 
		dnf_model = dtox_nn.DToxNet(dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_input_size)
		# load trained DTox model parameter into defined structure 
		dnf_state = torch.load(dnf)
		dnf_model.load_state_dict(dnf_state['model_state_dict']) 
		# implement LRP to compute the relevance of all modules in null DTox model 
		_, dnf_pathway_relevance_df = dtox_lrp.lrp(dnf_model, query_feature_data, argv[12], float(argv[13]), float(argv[14]), 0, 1)
		dnf_pathway_relevance_df.index = query_data_df.index
		# compute null relevance scores for all DTox paths 
		dnf_path_relevance_df = dtox_lrp.compute_path_relevance_score(dnf_model, dnf_pathway_relevance_df)
		null_relevance_list.append(dnf_path_relevance_df)
	# compute empirical p-value of relevance scores for all DTox path-compound pairs
	relevance_fdr_df = dtox_lrp.compute_path_relevance_pvalue(path_relevance_df, null_relevance_list) 
	# iterate by path, substitute node number with node name   
	relevance_path_names = []
	for rpv in relevance_fdr_df.path_id.values:
		relevance_path_names.append(path_name_dict[rpv])
	relevance_fdr_df.path_id = relevance_path_names
	# write data frame that contains p-values of DTox path-compound pairs to output tsv file 
	relevance_fdr_df.to_csv(argv[15] + '_path_relevance_pv.tsv', sep = '\t', float_format = '%.5f', index = False)

	return 1


## Call main function 
main(sys.argv)
