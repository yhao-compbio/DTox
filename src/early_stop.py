# !/usr/bin/env python 
## created by Yun Hao @MooreLab 2021
## code inspired by Bjarten's github repo https://github.com/Bjarten/early-stopping-pytorch/blob/master/pytorchtools.py
## This script contains early stop function of DTox model


## Module
import numpy as np
import torch


## This function stops training of neural network if testing loss doesn't improve after a given patience
class stop:
	## 0. Input arguments 
		# patience: (int): how many epochs to wait after last time testing loss improved. Default: 20
		# delta: (float): minimum change in the monitored quantity to qualify as an improvement. Default: 0
		# model_name: (str): name of model to be saved
		# verbose: (bool): If True, prints a message for each testing loss improvement. 

	## 1. Define stop parameters 
	def __init__(self, patience = 20, delta = 0, model_name = 'checkpoint.pt', verbose = False):
		self.patience = patience
		self.delta = delta
		self.model_name = model_name
		self.verbose = verbose
		self.counter = 0
		self.best_score = None
		self.early_stop = False
		self.val_loss_min = np.Inf

	## 2. Check whether early stopping criterion is reached   
	def __call__(self, val_loss, model, optimizer):
		# use negative loss score so that the trend of score to observe is increasing
		score = -val_loss
		# if at first epoch, assign best score, save current model 
		if self.best_score is None:
			self.best_score = score
			self.save_checkpoint(val_loss, model, optimizer)
		# if score does not exceed current best after adjusting by minimum change
		elif score <= self.best_score + self.delta:
			self.counter += 1 
			if self.verbose:
				print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
			# stop training when number of patience is reached   
			if self.counter >= self.patience:
				self.early_stop = True
		# if score exceeds current best, re-assign best score, save current model, reset patience count to 0  
		else:
			self.best_score = score
			self.save_checkpoint(val_loss, model, optimizer)
			self.counter = 0

	## 3. Save current model when testing loss is still improving 
	def save_checkpoint(self, val_loss, model, optimizer):
		if self.verbose:
			print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
		torch.save({'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict()}, self.model_name)
		self.val_loss_min = val_loss
