# !/usr/bin/env python
## created by Yun Hao @MooreLab 2019
## This script contains the loss function used in DTox model


## Module
import torch
import torch.nn as nn


## This function defines a hybrid function of DTox model, which combines root and auxiliary loss.
class DToxHybridLoss(nn.Module):
	## 0. Input arguments
		# alpha: coefficient for auxiliary loss in loss function 
		# layer_size: list/array that  contains number of node modules in the corresponding hidden layer of each node 

	## 1. Define loss function parameters 
	def __init__(self, alpha, layer_size):	
		super(DToxHybridLoss, self).__init__()
		self.alpha = alpha
		self.BCELoss = nn.BCELoss()
		self.layer_size = layer_size
	
	## 2. Define loss function
	def forward(self, y, root_pred_y, auxiliary_pred_y):
		# compute the first part of loss function: binary cross entropy loss of final root output 
		root_loss = self.BCELoss(root_pred_y, y)
		# iterate by order of node module to compute the second part of loss function: binary cross entropy loss of auxiliary output
		auxi_loss = 0
		for lapy in range(0, len(auxiliary_pred_y)):
			auxi_loss += 1/float(self.layer_size[lapy]) * self.BCELoss(auxiliary_pred_y[lapy], y) 
		# multiply auxiliary loss by coefficient, add two parts together 
		total_loss = root_loss + self.alpha * auxi_loss
		return total_loss, root_loss, auxi_loss

