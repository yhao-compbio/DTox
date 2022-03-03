# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions used to process sorted DTox hiearchy files and compute model statistics 


## Module
import numpy as np
import pandas as pd


## This function computes the input and output size of each node module in DTox model, i.e. the number of neurons each module contains 
def compute_module_size(node_size, node_children, min_path, max_module, node_layer, root_node, log_scale = True):
	## 0. Input arguments 
		# node_size: data frame that stores number of annotated genes in nodes of DTox model   
		# node_children: dictionary that stores children nodes of each parent node in DTox model
		# min_path: minimal size of pathways included in DTox model
		# max_module: maximal size of node modules (number of neurons the module contains) in DTox model
		# node_layer: data frame that stores layer number of nodes in DTox model
		# root_node: node number of root pathways in DTox model
		# log_scale: whether to use log scale when adjusting for min/max pathway size
	
	## 1. Compute the output size of each node module
	# obtain maximal size of pathways
	node_size_vec = node_size['size'].values 
	max_node_size = np.max(node_size_vec)
	# iterate by each node module   
	module_size = []
	for nsv in node_size_vec:
		# for gene modules, each module contains 1 neuron  
		if nsv < min_path:
			module_size.append(1)
		# for pathway modules, the number of neurons varies between 1 and 'max_module', depends on the pathway size 
		else:
			# adjust pathway size by min/max in log scale
			if log_scale:
				nsv_module = 1 + (max_module - 1) * np.log(nsv/min_path)/np.log(max_node_size/min_path)
			# adjust pathway size by min/max
			else:
				nsv_module = 1 + (max_module - 1) * (nsv-min_path)/(max_node_size-min_path)	 
			module_size.append(int(round(nsv_module)))
	
	## 2. Compute the input size of each node module using parent/children node relationships, as well as the size of root module 
	# aggregate node information with module size  
	info_df = pd.merge(node_size, node_layer, on = 'node')
	info_df['module_size'] = module_size
	# iterate by each node module  
	in_module = []
	out_module = []
	out_root = 0
	for i in range(0, info_df.shape[0]):
		# ignore nodes in layer 0 (input layer) 
		if info_df.layer_number[i] == 0:
			continue
		# nodes in layer 1+ (hidden layers) 
		else:
			# obtain children nodes of current node
			i_children_info = info_df[info_df.node.isin(node_children[i])]
			# input module size is the sum of output module size of children nodes 
			in_module.append(i_children_info.module_size.sum())
			# output module size 
			out_module.append(info_df.module_size[i])
			# add output module size to root module size if the node contains a root pathway  
			if info_df.node[i] in root_node:
				out_root = out_root + info_df.module_size[i]

	return in_module, out_module, out_root


## This function computes basic statistics of DTox neural network model 
def compute_hierarchy_statistics(node_layer, out_module):
	## 0. Input arguments 
		# node_layer: data frame that contains layer number of nodes in DTox model 
		# out_module: list/array that contains output size of each node in DTox model    

	## 1. Compute basic statistics of DTox model
	# obtain number of hidden layers   
	N_hidden_layer = node_layer.layer_number.max()
	# obtain number of input features 
	N_input_node = node_layer[node_layer.layer_number == 0].shape[0]
	# obtain the total number of gene/pathway node modules in each hidden layer 
	node_layer1 = node_layer[node_layer.layer_number > 0]
	layer_size_df = pd.DataFrame(node_layer1.layer_number.value_counts(sort = False))
	layer_size_df.columns = ['layer_size']
	N_hidden_node = layer_size_df.layer_size.values
	# obtain the total number of gene/pathway node modules in the corresponding hidden layer of each node 
	node_layer_size = pd.merge(node_layer1, layer_size_df, left_on = 'layer_number', right_index = True)	
	N_layer_size = node_layer_size.layer_size.values
	# obtain the total number of neurons in each hidden layer  
	node_layer1.insert(node_layer1.shape[1], 'module_size', out_module)
	N_hidden_neuron = node_layer1.groupby(['layer_number']).module_size.sum().values
	
	## 2. Add computed statistics to output strings  
	h_stat = []
	h_stat.append('Number of hidden layers:' + str(N_hidden_layer))
	h_stat.append('Number of input features: ' + str(N_input_node))
	h_stat.append('Number of hidden modules: ' + ','.join([str(nhn) for nhn in N_hidden_node]))
	h_stat.append('Number of hidden neurons: ' + ','.join([str(nhn) for nhn in N_hidden_neuron]))

	return h_stat, N_input_node, N_layer_size	


## This function reads in sorted hierarchy files, then computes the size of node modules in DTox model  
def sort_dtox_hierarchy(root_file, relation_file, node_size_file, layer_file, min_path_size, max_module_size):
	## 0. Input arguments 
		# root_file: input file that contains root pathways of DTox model 
		# relation_file: input file that contains parent/children node connections of DTox model 
		# node_size_file: input file that contains number of annotated genes in nodes of DTox model 
		# layer_file: input file that contains layer number of nodes in DTox model
		# min_path_size: minimal size of pathways included in DTox model
		# max_module_size: maximal size of node modules in DTox neural network model

	## 1. Obtain node number of root pathways 
	root_df = pd.read_csv(root_file, sep = '\t')
	root = root_df.root.values

	## 2. Obtain children node of each parent node, store information in a dictionary
	relation_df = pd.read_csv(relation_file, sep = '\t')
	# iterate by parent node  
	node_children_dict = {}
	for i in range(0, relation_df.shape[0]):
		# obtain numbers of current parent node and its children node(s) 
		i_node = relation_df.node.values[i] 
		i_children = relation_df.children_node.values[i].split(',')
		node_children_dict[i_node] = [int(ic) for ic in i_children]	

	## 3. Compute the input and output size of each node module in DTox model 
	size_df = pd.read_csv(node_size_file, sep = '\t')
	layer_df = pd.read_csv(layer_file, sep = '\t')
	input_module_size, output_module_size, output_root_size = compute_module_size(size_df, node_children_dict, min_path_size, max_module_size, layer_df, root)	

	# 4. Compute basic statistics of DTox neural network model
	hierarchy_stat, input_layer_size, hidden_layer_size = compute_hierarchy_statistics(layer_df, output_module_size)

	return node_children_dict, input_module_size, output_module_size, root, output_root_size, hierarchy_stat, input_layer_size, hidden_layer_size

