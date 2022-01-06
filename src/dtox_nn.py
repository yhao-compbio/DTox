# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions used to build basic neural network structure for DTox model


## Module
import torch
from torch import nn

 
## This function defines a root loss-oriented neural network structure for each node module  
class DToxNetModule(nn.Module):
	## 0. Input arguments 
		# M_in: number of input neurons 
		# M_out: number of output neurons 
	
	## 1. Define module parameters  
	def __init__(self, M_in, M_out):
		super(DToxNetModule, self).__init__()
		self.linear = nn.Linear(M_in, M_out)
		self.activate = nn.ReLU()
	
	## 2. Define forward function 	
	def forward(self, x):
		# linear transformation
		m_linear = self.linear(x)
		# nonlinear ReLU activation
		m_active = self.activate(m_linear)
		return m_active


## This function defines an anxiliary loss-oriented neural network structure for each node module
class DToxNetAuxiliary(nn.Module):
	## 0. Input arguments
		# M_in: number of input neurons 

	## 1. Define module parameters  
	def __init__(self, M_in):
		super(DToxNetAuxiliary, self).__init__()
		self.linear = nn.Linear(M_in, 1)
		self.sigmoid = nn.Sigmoid()

	## 2. Define forward function 
	def forward(self, x):
		# linear transformation
		m_linear = self.linear(x)
		# nonlinear sigmoid activation
		m_sigmoid = self.sigmoid(m_linear)
		return m_sigmoid


## This function defines the structure of whole neural network for DTox model
class DToxNet(nn.Module):
	## 0. Input arguments
		# node_children_dict: dictionary that contains children nodes of each parent node in DTox model 
		# input_module_size: list/array that contains the input size of each node module in DTox model 
		# output_module_size: list/array that contains the input size of each node module in DTox model  
		# root: list/array that contains number IDs of root pathways in DTox model
		# output_root_size: size of root module in DTox model 
		# input_feature_size: number of input features 

	## 1. Define neural network parameters  
	def __init__(self, node_children_dict, input_module_size, output_module_size, root, output_root_size, input_feature_size):
		super(DToxNet, self).__init__()
		# iterate by order of node module in DTox model
		self.net = nn.ModuleList()
		self.auxiliary = nn.ModuleList()
		for lims in range(0, len(input_module_size)):
			# define root loss-oriented structure for the current node 
			net_module = DToxNetModule(input_module_size[lims], output_module_size[lims])
			self.net.append(net_module)
			# define auxiliary loss-oriented structure for the current node  
			auxiliary_module = DToxNetAuxiliary(output_module_size[lims])
			self.auxiliary.append(auxiliary_module)
		# define structure for root module
		self.output_layer = DToxNetAuxiliary(output_root_size)
		# define other parameters
		self.node_children = node_children_dict
		self.root_node = root
		self.input_size = input_feature_size
		self.hidden_size = len(input_module_size)
		self.combine_size = input_feature_size + len(input_module_size)
	
	## 2. Define forward function 
	def forward(self, x):
		# iterate by order of input features, assign feature values to result list 
		layer_result = [0] * self.combine_size
		auxiliary_result = [0] * self.hidden_size
		for sis in range(0, self.input_size):
			# assign values of current feature to result list
			layer_result[sis] = x[:, sis:(sis+1)]	
		# iterate by order of hidden node modules  
		for scs in range(self.input_size, self.combine_size):
			# obtain input values for the current node  
			scs_children = self.node_children[scs]
			scs_children_output_list = [layer_result[sc] for sc in scs_children]
			scs_children_output = torch.cat(scs_children_output_list, 1)
			# feed input values into root loss-oriented structure to compute node output 
			scs_net_id = scs - self.input_size
			layer_result[scs] = self.net[scs_net_id](scs_children_output)
			# feed input values into auxiliary loss-oriented structure to compute auxiliary node output
			auxiliary_result[scs_net_id] = self.auxiliary[scs_net_id](layer_result[scs])
		# obtain values for root node  
		root_output_list = [layer_result[srn] for srn in self.root_node]
		root_output = torch.cat(root_output_list, 1)
		# feed values into root structure to compute final output   
		y_pred = self.output_layer(root_output)
		return y_pred, auxiliary_result
