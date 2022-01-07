# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains deep learning functions used in the DTox model construction


## Module
import sys
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score
sys.path.insert(0, 'src/')
import dtox_hierarchy
import dtox_nn
import dtox_loss
import early_stop


## This function trains a DTox neural network model based on input training data and specified hyperparameters 
def train_dtox_model(dtox_root_file, dtox_relation_file, dtox_node_size_file, dtox_layer_file, dtox_min_path_size, dtox_max_module_size, dtox_train_data_loader, dtox_test_data, dtox_alpha, dtox_lambda, dtox_learning_rate = 0.001, dtox_patience = 20, dtox_max_epoch = 200, dtox_model_name = 'checkpoint.pt', dtox_device = 'cpu'):
	## 0. Input arguments 
		# dtox_root_file: input file that contains root pathways of DTox model 
		# dtox_relation_file: input file that contains parent/children node connections of DTox model  
		# dtox_node_size_file: input file that contains number of annotated genes in nodes of DTox model   
		# dtox_layer_file: input file that contains layer number of nodes in DTox model
		# dtox_min_path_size: minimal size of pathways included in DTox model    
		# dtox_max_module_size: maximal size of node modules in DTox neural network model
		# dtox_train_data_loader: torch data loader that contains formatted training data (see 'dtox_data')
		# dtox_test_data: formatted testing data (see 'dtox_data')
		# dtox_alpha: coefficient for auxiliary loss in loss function (see 'dtox_loss')
		# dtox_lambda: coefficient for L2 regularization 
		# dtox_learning_rate: learning rate for DTox neural network model
		# dtox_patience: maximal number of epochs to run after testing loss stop decreasing, i.e. optimal model is reached 
		# dtox_max_epoch: maximal number of epochs to run before stopping if an optimal testing loss cannot be reached 
		# dtox_model_name: output file that stores DTox neural network model
		# dtox_device: device to train DTox model on, 'cpu' or 'cuda' 
	
	## 1. Construct DTox neural network model using the specified hyperparameters 
	# sort the hierarchy of DTox models, compute module size and basic statistics 
	dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_hierarchy_stat, dtox_input_size, dtox_hidden_size = dtox_hierarchy.sort_dtox_hierarchy(dtox_root_file, dtox_relation_file, dtox_node_size_file, dtox_layer_file, dtox_min_path_size, dtox_max_module_size) 
	# define structure of whole neural network 
	model = dtox_nn.DToxNet(dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_input_size) 
	model.to(dtox_device)
	# define loss function 
	loss_function = dtox_loss.DToxHybridLoss(dtox_alpha, dtox_hidden_size)
	loss_function.to(dtox_device)
	# use Adam optimizer with weight decay (L2 regularization)
	optimizer = torch.optim.Adam(model.parameters(), lr = dtox_learning_rate, weight_decay = dtox_lambda)
	# define early stop function
	stop_function = early_stop.stop(patience = dtox_patience, model_name = dtox_model_name)
	
	## 2. Train DTox model 
	# Perform training by epoch until early stopping criterion is reached  
	epoch = 0
	train_total_loss = []
	train_root_loss = []
	train_auxi_loss = []
	test_total_loss = []
	test_root_loss = []
	test_auxi_loss = []
	while stop_function.early_stop == False:
		epoch += 1
		# set model to training mode
		model.train()
		# iterate by mini-batch, perform forward and backward propogation, keep track of training loss 
		current_train_total_loss = 0
		current_train_root_loss = 0
		current_train_auxi_loss = 0
		for i, batch_data in enumerate(dtox_train_data_loader, 0):
			# get feature and response data of current batch
			batch_feature, batch_label = batch_data[0].to(dtox_device), batch_data[1].to(dtox_device)
			# set the gradients to 0
			optimizer.zero_grad()
			# perform forward propogation to compute root and auxiliary output 
			batch_y_pred, batch_auxiliary_pred = model(batch_feature)
			# compute total, root, and auxiliary loss of current batch
			batch_total_loss, batch_root_loss, batch_auxi_loss = loss_function(batch_label, batch_y_pred, batch_auxiliary_pred)
			current_train_total_loss += float(batch_total_loss.data)
			current_train_root_loss += float(batch_root_loss.data)
			current_train_auxi_loss += float(batch_auxi_loss.data)
			# perform backward propogation
			batch_total_loss.backward()
			# perform optimization
			optimizer.step()
		# average three types of training loss over all mini-batches, store the average  
		current_train_total_loss_mean = current_train_total_loss/(i+1)
		current_train_root_loss_mean = current_train_root_loss/(i+1)
		current_train_auxi_loss_mean = current_train_auxi_loss/(i+1)
		train_total_loss.append(current_train_total_loss_mean)
		train_root_loss.append(current_train_root_loss_mean)
		train_auxi_loss.append(current_train_auxi_loss_mean)
		# set model to evaluation mode
		model.eval()
		# implement current model to testing data, perform forward propogation to compute root and auxiliary output 
		test_feature, test_label = dtox_test_data.features.to(dtox_device), dtox_test_data.labels.to(dtox_device) 
		test_y_pred, test_auxi_pred = model(test_feature)
		# compute total, root, and auxiliary loss of testing output
		current_test_total_loss1, current_test_root_loss1, current_test_auxi_loss1 = loss_function(test_label, test_y_pred, test_auxi_pred)
		current_test_total_loss = float(current_test_total_loss1.data)
		current_test_root_loss = float(current_test_root_loss1.data)
		current_test_auxi_loss = float(current_test_auxi_loss1.data)
		test_total_loss.append(current_test_total_loss)
		test_root_loss.append(current_test_root_loss)
		test_auxi_loss.append(current_test_auxi_loss)
		# check if early stop criterion has been met 
		stop_function(current_test_total_loss, model, optimizer)
		# if so, load the last checkpoint with the best model
		if stop_function.early_stop:
			stop_point_state = torch.load(dtox_model_name)
			model.load_state_dict(stop_point_state['model_state_dict'])
			optimizer.load_state_dict(stop_point_state['optimizer_state_dict'])
			break
		# stop training if the maximum epoch is reached  	
		if epoch == dtox_max_epoch:
			break
	
	## 3. Store training and testing loss of every epoch in data frame form  
	train_epoch = np.arange(1, epoch+1)
	train_summary = pd.DataFrame({'epoch': train_epoch, 'training_total_loss': train_total_loss, 'training_root_loss': train_root_loss, 'training_auxiliary_loss': train_auxi_loss, 'testing_total_loss': test_total_loss, 'testing_root_loss': test_root_loss, 'testing_auxiliary_loss': test_auxi_loss})

	return dtox_hierarchy_stat, model, loss_function, train_summary


## This function computes confidence interval width of metric by bootstrapping
def compute_metric_ci_by_bootsrap(metric_function, label_vec, pred_vec, confidence_interval = 0.95, bootstrap_times = 1000):
	## 0. Input arguments: 
		# metric_function: scikit-learn metric function to evalute classification model
		# label_vec: list/array that conatins true sample labels 
		# pred_vec: list/array that conatins positive label predicted probability of samples 
		# confidence_interval: confidence interval ratio to be computed (number between 0 and 1)
		# bootstrap_times: repeated sampling times for bootstrap

	## 1. Compute confidence interval of mean by bootstrapping
	vec_len = len(pred_vec)
	id_vec = np.arange(0, vec_len)
	# repeat boostrap process
	sample_metrics = []
	np.random.seed(0)
	for sample in range(0, bootstrap_times):
		# sampling with replacement from the input array
		sample_ids = np.random.choice(id_vec, size = vec_len, replace = True)
		sample_ids = np.unique(sample_ids)
		# compute sample metric
		sample_metric = metric_function(label_vec[sample_ids], pred_vec[sample_ids])
		sample_metrics.append(sample_metric)
	# sort means of bootstrap samples 
	sample_metrics = np.sort(sample_metrics)
	# obtain upper and lower index of confidence interval 
	lower_id = int((0.5 - confidence_interval/2) * bootstrap_times) - 1
	upper_id = int((0.5 + confidence_interval/2) * bootstrap_times) - 1
	# compute width of confidence interval
	ci = (sample_metrics[upper_id] - sample_metrics[lower_id])/2
	
	return ci


## This function evaluates the performance of trained DTox model on input validation data  
def evaluate_dtox_model(dtox_model, dtox_loss_function, dtox_eval_data, dtox_device = 'cpu'): 
	## 0. Input arguments
		# dtox_model: trained DTox model
		# dtox_loss_function: defined DTox loss function  
		# dtox_eval_data: formatted validation data (see 'dtox_data')
		# dtox_device: device to train DTox model on, 'cpu' or 'cuda' 

	## 1. Implement trained DTox model on validation data to generate predicted output
	# set model to evaluation mode
	dtox_model.eval()
	# implement trained DTox model to validation data, perform forward propogation to compute root and auxiliary output 
	eval_feature, eval_label = dtox_eval_data.features.to(dtox_device), dtox_eval_data.labels.to(dtox_device)
	eval_y_pred, eval_auxi_pred = dtox_model(eval_feature)
	# convert predicted probability and label of validation data to numpy array (for computing metrics)
	eval_pred = np.array(eval_y_pred.data.cpu()).flatten()
	eval_label1 = np.array(eval_label.cpu()).flatten()
	# convert predicted probability into binary predicted label
	eval_pred_label = []
	for ep in range(0, len(eval_pred)):
		if eval_pred[ep] < 0.5:
			eval_pred_label.append(0)
		else:
			eval_pred_label.append(1)
	eval_pred_label = np.array(eval_pred_label) 
	
	## 2. Compute model evaluation metrics 
	# compute total, root, and auxiliary loss of validation output
	eval_total_loss1, eval_root_loss1, eval_auxi_loss1 = dtox_loss_function(eval_label, eval_y_pred, eval_auxi_pred)
	eval_total_loss = float(eval_total_loss1.data)
	eval_root_loss = float(eval_root_loss1.data)
	eval_auxi_loss = float(eval_auxi_loss1.data)
	# compute AUROC and its 95% confidence interval
	eval_auc = roc_auc_score(eval_label1, eval_pred)
	eval_auc_ci = compute_metric_ci_by_bootsrap(roc_auc_score, eval_label1, eval_pred)
	# compute balanced accuracy and its 95% confidence interval
	eval_bac = balanced_accuracy_score(eval_label1, eval_pred_label)
	eval_bac_ci = compute_metric_ci_by_bootsrap(balanced_accuracy_score, eval_label1, eval_pred_label)
	# compute F1 score and its 95% confidence interval
	eval_f1 = f1_score(eval_label1, eval_pred_label)
	eval_f1_ci = compute_metric_ci_by_bootsrap(f1_score, eval_label1, eval_pred_label)
	# store computed metrics in data frame form
	metric_dict = {'total_loss': eval_total_loss, 'root_loss': eval_root_loss, 'auxiliary_loss': eval_auxi_loss, 'auc': eval_auc, 'auc_ci': eval_auc_ci, 'bac': eval_bac, 'bac_ci': eval_bac_ci, 'f1': eval_f1, 'f1_ci': eval_f1_ci}

	return metric_dict


## This function generates output list of DTox model performance 
def generate_dtox_performance_file(dtox_train_data, dtox_train_perf, dtox_valid_data, dtox_valid_perf, dtox_hierarchy_stat, round_digit = 5):
	## 0. Input arguments 
		# dtox_train_data: formatted training data (see 'dtox_data')
		# dtox_train_perf: dictionary that contains training performance of DTox model
		# dtox_valid_data: formatted validation data (see 'dtox_data')
		# dtox_valid_perf: dictionary that contains validation performance of DTox model
		# dtox_hierarchy_stat: list that contains computed basic staistics of DTox model structure  
		# round_digit: number of decimal places to round to (default: 5) 

	## 1. Generate training performance string  
	# obtain number of training samples  
	N_train = len(dtox_train_data)
	# iterate by metric, convert training performance to string  
	train_perf_str = []
	for k,v in dtox_train_perf.items():
		v = round(v, round_digit)
		train_perf_str.append(k + ':' + str(v))	

	## 2. Generate validation performance string 
	# obtain number of validation samples 
	N_valid = len(dtox_valid_data)
	# iterate by metric, convert validation performance to string  
	valid_perf_str = []
	for k,v in dtox_valid_perf.items():
		v = round(v, round_digit)
		valid_perf_str.append(k + ':' + str(v)) 
		 
	## 3. Generate performance file  
	dtox_perf_stat = dtox_hierarchy_stat.copy()
	dtox_perf_stat.append('Number of training instances: ' + str(N_train))
	dtox_perf_stat.append('Training performance: ' + ','.join(train_perf_str))
	dtox_perf_stat.append('Number of validation instances: ' + str(N_valid))
	dtox_perf_stat.append('Validation performance: ' + ','.join(valid_perf_str))
	
	return dtox_perf_stat
