# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions used in the Multi-Layer Perceptron (MLP) neural network model


## Module
import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data.dataset import Dataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score
sys.path.insert(0, 'src/')
import early_stop


## This function defines formats feature-response array data into tensors used for MLP model training  
class MLP_dataformat(Dataset):
	## 0. Input arguments 
		# feature_data: array that contains input feature data 
		# label_data: array that contains input label/response data

	## 1. Convert feature and label arrays into tensors 
	def __init__(self, feature_data, label_data):
		super(MLP_dataformat, self).__init__()
		self.features = torch.tensor(feature_data, dtype = torch.float)
		labels = torch.tensor(label_data, dtype = torch.float)
		self.labels = labels.view(labels.shape[0], 1)

	## 2. Get feature and label data by index
	def __getitem__(self, index):
		feature = self.features[index]
		label = self.labels[index]
		return feature, label

	## 3. Obtain number of data samples 
	def __len__(self):
		return len(self.labels)


## This function defines the structure of whole neural network for MLP model 
class MLP(nn.Module):
	## 0. Input arguments 
		# input_size: number of input features 
		# hidden_size: list/array that contains specified numbers of hidden neurons (in the order of hidden layer) 
	
	## 1. Define neural network parameters  
	def __init__(self, input_size, hidden_size):
		super(MLP, self).__init__()
		# put size of hidden layer  
		all_size = np.insert(hidden_size, 0, input_size)
		all_layers = []
		# if no hidden layer, apply linear transformation to input feature values (logistic regression) 
		if len(all_size) == 1: 
				all_layers.append(nn.Linear(input_size, 1))
		# iterate by hidden layer
		else:
			for i in range(0, len(all_size)-1):
				# apply linear transformation between current hidden layer and next one 
				all_layers.append(nn.Linear(all_size[i], all_size[i+1]))
				# apply ReLU activation function
				all_layers.append(nn.ReLU())
			# apply linear transformation between last hidden layer and output layer
			all_layers.append(nn.Linear(all_size[-1], 1))
		# apply sigmoid activation function to get final output 
		all_layers.append(nn.Sigmoid())
		# sort all nn functions in the added order  
		self.layers = nn.Sequential(*all_layers)
	
	## 2. Define forward function 	
	def forward(self, x):
		y = self.layers(x)
		return y


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


## This function evaluates the performance of trained MLP model on input validation data 
def evaluate_mlp_model(mlp_model, mlp_loss, mlp_eval_data):
	## 0. Input arguments
		# mlp_model: trained MLP model
		# mlp_loss: defined MLP loss function
		# mlp_eval_data: formatted validation data

	## 1. Implement trained MLP model on validation data to generate predicted output 
	# set model to evaluation mode 
	mlp_model.eval()
	# implement trained MLP model to validation data, perform forward propogation to compute output  
	eval_feature, eval_label = mlp_eval_data.features, mlp_eval_data.labels
	eval_y_pred = mlp_model(eval_feature)
	# convert predicted probability and label of validation data to numpy array (for computing metrics)
	eval_pred = np.array(eval_y_pred.data).flatten()
	eval_label1 = np.array(eval_label).flatten()
	eval_pred_label = []
	# convert predicted probability into binary predicted label
	for ep in range(0, len(eval_pred)):
		if eval_pred[ep] < 0.5:
			eval_pred_label.append(0)
		else:
			eval_pred_label.append(1)
	eval_pred_label = np.array(eval_pred_label)

	## 2. Compute model evaluation metrics 
	# compute loss of validation output   
	eval_loss = mlp_loss(eval_y_pred, eval_label)
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
	metric_dict = {'loss': np.float64(eval_loss.data), 'auc': eval_auc, 'auc_ci': eval_auc_ci, 'bac': eval_bac, 'bac_ci': eval_bac_ci, 'f1': eval_f1, 'f1_ci': eval_f1_ci}
	
	return metric_dict


## This function trains and evaluates an MLP neural network model 
def train_mlp_model(X_combine, X_valid, y_combine, y_valid, N_hidden_layer_nodes, test_prop = 0.125, N_batch = 32, learning_rate = 0.001, patience = 20, max_epoch = 200, model_name = 'checkpoint.pt'):
	## 0. Input argument: 
		# X_combine: array that contains feature values of combined training and testing data 
		# X_valid: array that contains feature values of validation data
		# y_combine: array that contains response/label values of combined training and testing data  
		# y_valid: array that contains response/label values of validation data 
		# N_hidden_layer_nodes: list/array that contains specified numbers of hidden neurons (in the order of hidden layer) 
		# test_prop: proportion of testing samples among combined training and testing data 
		# N_batch: number of mini-batches 
		# learning_rate: learning rate for MLP neural network model 
		# patience: maximal number of epochs to run after testing loss stop decreasing, i.e. optimal model is reached 
		# max_epoch: maximal number of epochs to run before stopping if an optimal testing loss cannot be reached 
		# model_name: output file that stores MLP neural network model

	## 1. Format training, testing, and validation data
	# format feature and label data of combined training and testing data  
	mlp_combine_data = MLP_dataformat(X_combine, y_combine)
	# split combined training and testing data into training and testing according to test_prop
	X_mlp_train, X_mlp_test, y_mlp_train, y_mlp_test = train_test_split(X_combine, y_combine, test_size = test_prop, random_state = 0, stratify = y_combine)
	# format feature and label data of training data, then generate data loader accroding to N_batch 
	mlp_train_data = MLP_dataformat(X_mlp_train, y_mlp_train)
	mlp_train_data_loader = torch.utils.data.DataLoader(mlp_train_data, batch_size = N_batch, shuffle = True)	
	# format feature and label data of testing data  
	mlp_test_data = MLP_dataformat(X_mlp_test, y_mlp_test)
	# format feature and label data of validation data  
	mlp_valid_data = MLP_dataformat(X_valid, y_valid)

	## 2. Construct MLP neural network model using the specified hyperparameters
	# define structure of whole neural network  
	mlp_model = MLP(X_mlp_train.shape[1], N_hidden_layer_nodes)
	# define loss function 
	mlp_loss = nn.BCELoss()
	# use Adam optimizer
	mlp_optim = torch.optim.Adam(mlp_model.parameters(), lr = learning_rate) 
	# define early stop function
	mlp_stop = early_stop.stop(patience = patience, model_name = model_name)	

	## 3. Train MLP model
	# perform training by epoch until early stopping criterion is reached  
	epoch = 0
	train_loss = []
	test_loss = []
	while mlp_stop.early_stop == False:
		epoch += 1
		# set model to training mode 
		mlp_model.train()
		# iterate by mini-batch, perform forward and backward propogation, keep track of training loss 
		current_train_loss = 0 
		for i, batch_data in enumerate(mlp_train_data_loader, 0):
			# get feature and response data of current batch
			batch_feature, batch_label = batch_data
			# set the gradients to 0
			mlp_optim.zero_grad()
			# perform forward propogation to compute predicted output 
			batch_pred = mlp_model(batch_feature)
			# compute loss of current batch
			batch_loss = mlp_loss(batch_pred, batch_label)
			current_train_loss += float(batch_loss.data)	
			# perform backward propogation
			batch_loss.backward()
			# perform optimization
			mlp_optim.step()
		# average computed training loss over all mini-batches, store the average  
		current_train_loss_mean = current_train_loss/(i+1)
		train_loss.append(current_train_loss_mean)
		# implement current model to testing data, perform forward propogation to compute predicted output 
		test_y_pred = mlp_model(mlp_test_data.features)
		# compute loss of testing output
		current_test_loss = mlp_loss(test_y_pred, mlp_test_data.labels)
		test_loss.append(float(current_test_loss.data))
		# check if early stop criterion has been met 
		mlp_stop(float(current_test_loss.data), mlp_model, mlp_optim)
		# if so, load the last checkpoint with the best model
		if mlp_stop.early_stop:	
			stop_point_state = torch.load(model_name)
			mlp_model.load_state_dict(stop_point_state['model_state_dict'])
			mlp_optim.load_state_dict(stop_point_state['optimizer_state_dict'])
			break
		# stop training if the maximum epoch is reached    
		if epoch == max_epoch:
			break
	# store training and testing loss of every epoch in data frame form  	
	train_epoch = np.arange(1, epoch+1)
	train_summary = pd.DataFrame({'epoch': train_epoch, 'training_loss': train_loss,  'testing_loss': test_loss})
	
	## 4. Evaluate trained MLP model on combined trainig-testing and validation data
	combine_perf = evaluate_mlp_model(mlp_model, mlp_loss, mlp_combine_data)  
	valid_perf = evaluate_mlp_model(mlp_model, mlp_loss, mlp_valid_data)
	
	return mlp_model, train_summary, combine_perf, valid_perf


## This function converts query dictionary to a string
def convert_dict_to_string(query_dict, round_digit = 5):
	## 0. Input arguments: 
		# query_dict: dictionary that contains query dictionary
		# round_digit: number of decimal places to round to (default: 5) 

	## 1. Join names and values to build output strings
	# iterate by item in query_dict 
	query_str = []
	for k,v in query_dict.items():
		# round values  
		if type(v) is np.float64:
			v = np.round(v, round_digit)
		# convert values to strings  
		v_str = str(v)
		# join the name 
		query_str.append(k + ':' + v_str)
	# join all item strings together 
	output_str = ','.join(query_str)

	return output_str


## This function generates content for output MLP performance file   
def generate_mlp_performance_file(hidden_size, combine_data, valid_data, combine_metrics, valid_metrics):
	## 0. Input arguments
		# hidden_size: list/array that contains specified numbers of hidden neurons (in the order of hidden layer) 
		# combine_data: array that contains feature values of combined training and testing data 
		# valid_data: array that contains feature values of validation data 
		# combine_metrics: dictionary that contains the model performance on combined training-testing dataset
		# valid_metrics: dictionary that contains the model performance on validation dataset

	## 1. Obtain basic info about MLP model and dataset 
	# obtain number of hidden layers 
	N_hidden_layer = len(hidden_size)
	hidden_char = ','.join([str(hs) for hs in hidden_size])
	# obtain number of features and number of training-testing/validation samples 
	N_feature = combine_data.shape[1]
	N_combine = combine_data.shape[0]
	N_valid = valid_data.shape[0]
	
	## 2. Convert performance metric dictionaries to strings  
	# convert training metric dictionary keys and values to string  
	combine_metric_str = convert_dict_to_string(combine_metrics)		
	# convert validation metric dictionary keys and values to string 
	valid_metric_str = convert_dict_to_string(valid_metrics)

	## 3. Generate list of strings that describe model performance info 
	perf_list = []
	perf_list.append('Number of hidden layers: ' + str(N_hidden_layer))	
	perf_list.append('Number of input features: ' + str(N_feature))
	perf_list.append('Number of hidden neurons: ' + hidden_char)
	perf_list.append('Number of training instances: ' + str(N_combine))
	perf_list.append('Training performance: ' + combine_metric_str)
	perf_list.append('Number of testing instances: ' + str(N_valid))
	perf_list.append('Testing performance: ' + valid_metric_str)

	return perf_list
