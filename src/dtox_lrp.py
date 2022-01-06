# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions used for implementing LRP to evaluate relevance of DTox paths   


## Module
import numpy as np
import pandas as pd 
import torch
from functools import reduce
from operator import concat
from statsmodels.sandbox.stats.multicomp import multipletests


## This function computes a small increment as a smoothing factor (applied when denominator is 0) for a list of input tensors 
def compute_increment_values(m_list, incre_prop = 0.001):
	## 0. Input argument 
		# m_list: list of tensors 
		# incre_prop: small increment proportion of smallest absolute non-zero element

	## 1. Compute smoothing increment factor
	# iterate by tensor in the list 
	m_abs_min = torch.zeros(len(m_list)).double()
	for ml in range(0, len(m_list)):
		if len(m_list[ml]) > 0:
			# identify non-zero elements of current tensor
			m_list_n0 = m_list[ml][torch.nonzero(m_list[ml], as_tuple=True)]
			# if no non-zero element, move on to next one
			if len(m_list_n0) == 0:
				continue
			# if non-zero elements exist, compute the smallest absolute non-zero element
			else:
				m_abs_min[ml] = m_list_n0.abs().min() 
	# identify the smallest absolute non-zero element among all tensors, then compute increment   
	increment = m_abs_min[torch.nonzero(m_abs_min, as_tuple=True)].abs().min() * incre_prop

	return increment


## This function perform matrix multiplication while using smoothing factors to avoid 0 results. 
def fill_in_zero_matmul(m1, m2, m1_incre, m2_incre, m1_sign, m2_sign):
	## 0. Input arguments:
		# mat1: tensor that represents matrix 1
		# mat2: tensor that represents matrix 2
		# m1_incre: smoothing increment factor of matrix 1 
		# m2_incre: smoothing increment factor of matrix 2
		# m1_sign: sign of adjusting matrix 1 with increment (1: +; -1: -)
		# m2_sign: sign of adjusting matrix 2 with increment (1: +; -1: -)

	## 1. Perform matrix multiplication
	mul = torch.matmul(m1, m2)
	if (mul == 0).sum() > 0:
		# identify row and column id in which the multiplication results are 0 
		mul0_id = torch.where(mul == 0)
		# adjust identified rows of matrix 1 by smoothing increment factor   
		m1_incres = torch.zeros(m1.shape[0], 1).double()
		m1_incres[torch.unique(mul0_id[0]), :] = m1_incre
		m1 = m1 + m1_incres * m1_sign
		# adjust identified columns of matrix 1 by smoothing increment factor  
		m2_incres = torch.zeros(1, m2.shape[1]).double()
		m2_incres[:, torch.unique(mul0_id[1])] = m2_incre
		m2 = m2 + m2_incres * m2_sign
		# perform matrix multiplication again
		mul = torch.matmul(m1, m2)

	return m1, m2, mul


## This function implements generic (combined LRP-γ and LRP-ϵ) rule for propagating relevance from hidden layer k to j
def gamma_epsilon_rule(gamma, epsilon, r_k, a_j, weight_jk, a_incre, w_incre):
	## 0. Input arguments: 
		# gamma: coefficient for positive weights in LRP-γ rule 
		# epsilon: adjusting factor of denominator in LRP-ϵ rule 
		# r_k: tensor that contains relevance of neurons in hidden layer k  
		# a_j: tensor that contains neuron output after activation by hidden layer j
		# weight_jk: tensor that contains weight parameters between hidden layer k and j of learned DTox model 
		# a_incre: smoothing increment factor of a_j
		# w_incre: smoothing increment factor of weight_jk 

	## 1. Clone tensors into double tensors in order to improve the precision of calculation  
	r_k = r_k.clone().detach().double()
	a_j = a_j.clone().detach().double()
	weight_jk = weight_jk.clone().detach().double()

	## 2. Perform step 1 of calculation  
	weight_jk0 = torch.zeros_like(weight_jk)
	pos_weight_jk = torch.where(weight_jk > 0, weight_jk, weight_jk0)
	weight_jk1 = weight_jk + gamma * pos_weight_jk

	## 3. Perform step 2 of calculation
	a_j_new, weight_jk1_new, weight_j_sum = fill_in_zero_matmul(a_j, weight_jk1.T, a_incre, w_incre, 1, 1)
	weight_j_sum1 = weight_j_sum + epsilon * torch.std(weight_j_sum.double())	

	## 4. Perform step 3 of calculation
	weight_j_rk = r_k/weight_j_sum1 

	## 5. Perform step 4 of calculation
	r_j = torch.matmul(weight_j_rk, weight_jk1_new.T) * a_j_new

	return r_j


## This function implements αβ rule for propagating relevance from hidden layer k to j 
def alpha_beta_rule(alpha, beta, r_k, a_j, weight_jk, a_incre, w_incre):
	## 0. Input arguments:
		# alpha: coefficient for positive weights in LRP-αβ rule 
		# beta: coefficient for negative weights in LRP-αβ rule
		# r_k: tensor that contains relevance of neurons in hidden layer k  
		# a_j: tensor that contains neuron output after activation by hidden layer j
		# weight_jk: tensor that contains weight parameters between hidden layer k and j of learned DTox model 
		# a_incre: smoothing increment factor of a_j
		# w_incre: smoothing increment factor of weight_jk 

	## 1. Clone tensors into double tensors in order to improve the precision of calculation 
	r_k = r_k.clone().detach().double()
	a_j = a_j.clone().detach().double()
	weight_jk = weight_jk.clone().detach().double() 

	## 2. Perform step 1 of calculation 
	weight_jk0 = torch.zeros_like(weight_jk)
	pos_weight_jk = torch.where(weight_jk > 0, weight_jk, weight_jk0)
	neg_weight_jk = torch.where(weight_jk <= 0, weight_jk, weight_jk0)

	## 3. Perform step 2 of calculation 
	a_j_pos, pos_weight_jk_new, pos_weight_j_sum = fill_in_zero_matmul(a_j, pos_weight_jk.T, a_incre, w_incre, 1, 1)
	a_j_neg, neg_weight_jk_new, neg_weight_j_sum = fill_in_zero_matmul(a_j, neg_weight_jk.T, a_incre, w_incre, 1, -1)

	## 4. Perform step 3 of calculation 
	pos_weight_j_rk = r_k/pos_weight_j_sum 
	neg_weight_j_rk = r_k/neg_weight_j_sum 

	## 5. Perform step 4 of calculation 
	r_j = alpha * torch.matmul(pos_weight_j_rk, pos_weight_jk_new.T) * a_j_pos - beta * torch.matmul(neg_weight_j_rk, neg_weight_jk_new.T) * a_j_neg

	return r_j


## This function implements special rule for propagating relevance from first hidden layer to input layer 
def input_layer_rule(low, high, r_j, x_i, weight_ij, x_incre, w_incre):
	## 0. Input arguments 
		# low: lower bound of input feature values 
		# high: upper bound of input feature values 
		# r_j: tensor that contains relevance of neurons in first hidden layer (layer j)
		# x_i: tensor that contains feature value of input layer (layer i) 
		# weight_ij: tensor that contains weight parameters between input layer (i) and first hidden layer (j) of learned DTox model 
		# x_incre: smoothing increment factor of x_j
		# w_incre: smoothing increment factor of weight_ij

	## 1. Clone tensors into double tensors in order to improve the precision of calculation 
	r_j = r_j.clone().detach().double()
	x_i = x_i.clone().detach().double()
	weight_ij = weight_ij.clone().detach().double() 

	## 2. Perform step 1 of calculation 
	weight_ij0 = torch.zeros_like(weight_ij)
	pos_weight_ij = torch.where(weight_ij > 0, weight_ij, weight_ij0)
	neg_weight_ij = torch.where(weight_ij <= 0, weight_ij, weight_ij0)	

	## 3. Perform step 2 of calculation 
	x_low = x_i * 0 + low 
	x_high = x_i * 0 + high 
	x_i_new, weight_ij_new, weight_i_sum = fill_in_zero_matmul(x_i, weight_ij.T, x_incre, w_incre, 1, 1)
	weight_i_sum = weight_i_sum - torch.matmul(x_low, pos_weight_ij.T) - torch.matmul(x_high, neg_weight_ij.T)	

	## 4. Perform step 3 of calculation 
	weight_i_rj = r_j/(weight_i_sum + (weight_i_sum == 0).double() * 1e-6)

	## 5. Perform step 4 of calculation 
	r_i = x_i_new * torch.matmul(weight_i_rj, weight_ij_new.T) - x_low * torch.matmul(weight_i_rj, pos_weight_ij) - x_high * torch.matmul(weight_i_rj, neg_weight_ij)

	return r_i


## This function implements layer-wise relevance propagation (LRP) process to compute the relevance of all neurons/modules in learned DTox model 
def lrp(model, x, rule, rule_factor1, rule_factor2, input_low, input_high):
	## 0. Input arguments: 
		# model: DTox model with learned weight parameters 
		# x: tensor that contains input feature values  
		# rule: propagation rule to be implemented ('gamma-epsilon' for generic rule, 'alpha-beta' for αβ rule)
		# rule_factor1: first factor in LRP rules ('gamma' in generic rule, 'alpha' in αβ rule)  
		# rule_factor2: second factor in LRP rules ('epsilon' in generic rule, 'beta' in αβ rule)  
		# input_low: lower bound of input feature values for implementing input layer rule 
		# input_high: lower bound of input feature values for implementing input layer rule 

	## 1. Perform forward propagation based on learned weight parameters of DTox model
	# iterate by order of node module in DTox model, keep track of input/output tensors (value and size) of each module, as well as learned weight tensors 
	layer_children_id = [[] for i in range(model.combine_size+1)]
	layer_children_size = [[] for i in range(model.combine_size+1)]
	layer_input = [[] for i in range(model.combine_size+1)]
	layer_result = [[] for i in range(model.combine_size+1)]
	layer_weight = [[] for i in range(model.combine_size+1)]
	for mcs in range(0, model.combine_size):
		# when node belongs to input layer 
		if mcs < model.input_size:
			layer_result[mcs] = x[:, mcs:(mcs+1)]
		# when node belongs to hidden layer 
		else:
			# obtain children node IDs of current node 
			mcs_children = model.node_children[mcs]
			layer_children_id[mcs] = mcs_children
			# obtain children node sizes of current node  
			mcs_children_output_list = [layer_result[sc] for sc in mcs_children]
			layer_children_size[mcs] = [mcol.shape[1] for mcol in mcs_children_output_list]
			# obtain input tensor of current node  
			layer_input[mcs] = torch.cat(mcs_children_output_list, 1)
			# obtain output tensor of current node  
			mcs_net_id = mcs - model.input_size
			layer_result[mcs] = model.net[mcs_net_id](layer_input[mcs])
			# obtain weight tensor learned for current node 
			layer_weight[mcs] = model.net[mcs_net_id].linear.weight
	# obtain above information for root node  
	root_output_list = [layer_result[mrn] for mrn in model.root_node]
	layer_children_id[mcs+1] = model.root_node
	layer_children_size[mcs+1] = [rol.shape[1] for rol in root_output_list]
	layer_input[mcs+1] = torch.cat(root_output_list, 1)	
	layer_result[mcs+1] = model.output_layer(layer_input[mcs+1])
	layer_weight[mcs+1] = model.output_layer.linear.weight

	## 2. Perform backward propagation based on specified propagation rule
	# compute smoothing increment factor for node value and weight tensors separately    
	a_increment = compute_increment_values(layer_input)
	w_increment = compute_increment_values(layer_weight)
	# iterate by reverse order of node module in DTox model, keep track of computed neuron/module relevance scores 
	N_instance = x.shape[0]
	current_relevance_list = [[] for i in range(model.combine_size+1)]
	current_relevance_list[model.combine_size].append(layer_result[model.combine_size])
	neuron_relevance_list = [[] for i in range(model.combine_size+1)]
	module_relevance_list = [[] for i in range(model.combine_size+1)]
	for mcs in range(0, model.combine_size+1)[::-1]:
		# for nodes that are connected to parent node(s), sum the relevance score propagated from parent nodes of current node as neuron relevance score
		if len(current_relevance_list[mcs]) > 0:
			neuron_relevance_list[mcs] = torch.stack(current_relevance_list[mcs], dim = 0).sum(dim = 0)
		# for nodes in input layer that are connected to no parent node, assign relevance score of 0 
		else:
			neuron_relevance_list[mcs] = torch.zeros(N_instance, 1).double()
		# sum relevance score over all neurons of the current node as module relevance score  
		module_relevance_list[mcs] = neuron_relevance_list[mcs].sum(dim = 1).view(N_instance, 1)
		# for nodes in hidden layers, propagate neuron relevance score of current node backward to its children nodes 
		if mcs >= model.input_size:
			# for nodes in the second hidden layer and beyond, propagate according to the specified rule
			if layer_children_id[mcs][0] >= model.input_size:
				if rule == 'gamma-epsilon':
					mcs_children_relevance = gamma_epsilon_rule(rule_factor1, rule_factor2, neuron_relevance_list[mcs], layer_input[mcs], layer_weight[mcs], a_increment, w_increment)
				if rule == 'alpha-beta':
					mcs_children_relevance = alpha_beta_rule(rule_factor1, rule_factor2, neuron_relevance_list[mcs], layer_input[mcs], layer_weight[mcs], a_increment, w_increment)
			# for nodes in the first hidden layer, propagate according to the special input layer rule
			else:
				mcs_children_relevance = input_layer_rule(input_low, input_high, neuron_relevance_list[mcs], layer_input[mcs], layer_weight[mcs], a_increment, w_increment)
			# obtain the lower and upper column index of each child node in the relevance score matrix 
			mcs_children_len = len(layer_children_size[mcs])
			mcs_children_relevance_id = np.insert(layer_children_size[mcs], 0, 0).cumsum()
			# iterate by children nodes  
			for mcl in range(0, mcs_children_len):
				# assign the propagated neuron relevance score to the current child node
				mcl_id = layer_children_id[mcs][mcl]
				current_relevance_list[mcl_id].append(mcs_children_relevance[: , mcs_children_relevance_id[mcl]:mcs_children_relevance_id[mcl+1]])

	## 3. Aggregate propagated relevance scores of neurons and modules in the DTox model
	# iterate by order of node module in DTox model
	neuron_relevance_name_list = [[] for i in range(model.combine_size)]
	for i in range(0, model.combine_size):
		# name neurons in each module by two numbers separated by underscore: 1st number is the node number, 2nd number is the order of neuron in the node
		neuron_relevance_name_list[i] = [str(i) + '_' + str(j) for j in range(0, neuron_relevance_list[i].shape[1])]
	# aggragate propagated relevance scores of all neurons (including input layer) together, output in data frame form 
	neuron_relevance_df = pd.DataFrame(torch.cat(neuron_relevance_list[0:model.combine_size], dim = 1).detach().numpy())
	neuron_relevance_df.columns = reduce(concat, neuron_relevance_name_list)
	# aggragate propagated relevance scores of all modules (including input layer) together, output in data frame form
	module_relevance_df = pd.DataFrame(torch.cat(module_relevance_list[0:model.combine_size], dim = 1).detach().numpy())

	return neuron_relevance_df, module_relevance_df


## This function computes observed relevance scores for all DTox paths connecting input layer to root module   
def compute_path_relevance_score(model, module_relevance_df):
	## 0. Input arguments 
		# model: DTox model with learned weight parameters
		# module_relevance_df: data frame that contains 

	## 1. Obtain all paths of nodes under the DTox hierarchy  
	# start with root nodes, extend paths by children nodes until all paths reach the input layer 
	path = [[mrn] for mrn in model.root_node] 
	current_pos = 0
	while current_pos < len(path):
		current_path = path[current_pos]
		# if the last node of current path is in input layer, the path is finished, move on to next one 
		if current_path[-1] < model.input_size: 
			current_pos = current_pos + 1
		# if the last node of current path is not in input layer, extend by children of current last node   
		else: 
			# obtain the node numbers of current last node's children 
			current_children = model.node_children[current_path[-1]]
			# iterate by children nodes, update current path 
			for lcc in range(0, len(current_children)):
				# add current children to end of current path  
				new_path = current_path.copy()
				new_path.append(current_children[lcc])
				path.append(new_path)
			# remove the old path from list  
			path.remove(current_path)
	# iterate by path, name each path by nodes in the path  
	path_names = [] 
	for pa in path:
		# use underscore to link all nodes together in a string 
		pa_name = '_'.join([str(p) for p in pa])
		path_names.append(pa_name)

	## 2. Compute observed relevance scores for all DTox paths  
	# convert negative module relevance scores to 0 
	module_relevance_df1 = module_relevance_df.copy()
	module_relevance_df1[module_relevance_df1 < 0] = 0
	# perform log10 transformation  
	log_r_df = np.log10(module_relevance_df1) 
	# iterate by path, compute observed relevance score of each path  
	path_score_list = []
	for pa in path:
		# sum the log10-transformed relevance score of all nodes in the path
		pa_score = log_r_df.iloc[:, pa].sum(axis = 1)
		path_score_list.append(pa_score)
	# aggregate the computed scores of all paths in data frame form
	path_score_df = pd.concat(path_score_list, axis = 1) 
	path_score_df.columns = path_names

	return path_score_df


## This function computes empirical p-values of relevance scores for all DTox path-compound pairs
def compute_path_relevance_pvalue(observe_score_df, null_score_list, sig_only = True, fdr_cut = 0.05):
	## 0. Input arguments 
		# observe_score_df: data frame that contains computed observed scores of compounds (rows) for all DTox paths (columns) 
		# null_score_list: list of data frames that contains computed null scores of compounds (rows) for all DTox paths (columns)  
		# sig_only: whether only keep the significant paths of compounds 
		# fdr_cut: FDR threshold to determine when a path is significant (only applied when sig_only is True)

	## 1. Compute empirical p-values of relevance scores for all DTox path-compound pairs
	# coun the number of null relevance scores that exceed the observed score of each compound-path combination 
	null_compare_list = [(nsl >= observe_score_df) * 1 for nsl in null_score_list] 
	null_compare_df = reduce(lambda a, b: a.add(b, fill_value = 0), null_compare_list)
	pv_df = null_compare_df/len(null_score_list)

	## 2. Convert computed path p-values into data frame sorted by compounds 
	# iterate by path  
	N_osd = observe_score_df.shape[0]
	path_pv_list = []
	for osdc in observe_score_df.columns: 
		# obtain the computed relevance scores and p-values of all compounds for the current path  
		osdc_value = observe_score_df[osdc].values
		osdc_pv = pv_df[osdc].values
		# store compound CID, path name, observed relevance score, and p-value info in data frame form
		osdc_pv_df = pd.DataFrame({'cid': observe_score_df.index, 'path_id': np.repeat(osdc, N_osd), 'observed': osdc_value, 'p_value': osdc_pv})
		path_pv_list.append(osdc_pv_df)	
	# aggregate info data frame of all paths, sort data frame by compound CIDs  
	path_pv_df = pd.concat(path_pv_list, axis = 0)
	path_pv_df = path_pv_df.sort_values(['cid', 'p_value'])
	# correct for multiple testing by FDR, add adjusted p-value to data frame 
	_, path_fdr, _, _ = multipletests(path_pv_df['p_value'].values, method = 'fdr_bh')
	path_pv_df['fdr'] = path_fdr
	# remove non-significant paths of compounds, if specified  
	if sig_only: 
		path_pv_df = path_pv_df[path_pv_df['fdr'] < fdr_cut]

	return path_pv_df
