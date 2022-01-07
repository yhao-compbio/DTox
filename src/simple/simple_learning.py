# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains functions for building, evaluating, and implementing simple machine learning models


## Modules
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import log_loss
from sklearn.metrics import roc_auc_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score


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


## This function evaluate classification model by computing specified performance metric on input validation dataset 
def evaluate_classifier_by_metrics(classifier_model, X_eval, y_eval):
	## 0. Input arguments 
		# classifier_model: scikit-learn classfication model to be evaluted 
		# X_eval: array that contains feature values of validation dataset 
		# y_eval: array that contains label values of validation dataset
	
	## 1. Implement model to predict binary labels and probabilities of validation samples 
	y_pred_prob = classifier_model.predict_proba(X_eval)[:,1]
	y_pred = classifier_model.predict(X_eval)
	
	## 2. Compute evaluation metrics by comparing predicted label to true label
	# compute binary cross-entroy loss
	y_loss = log_loss(y_eval, y_pred_prob)
	# compute AUROC and its 95% confidence interval
	y_auc = roc_auc_score(y_eval, y_pred_prob)
	y_auc_ci = compute_metric_ci_by_bootsrap(roc_auc_score, y_eval, y_pred_prob)
	# compute balanced accuracy and its 95% confidence interval 
	y_bac = balanced_accuracy_score(y_eval, y_pred)
	y_bac_ci = compute_metric_ci_by_bootsrap(balanced_accuracy_score, y_eval, y_pred)
	# compute F1 score and its 95% confidence interval 
	y_f1 = f1_score(y_eval, y_pred)
	y_f1_ci = compute_metric_ci_by_bootsrap(f1_score, y_eval, y_pred)
	# store all computed performance metrics in a dictionary 
	metric_dict = {'log_loss': y_loss, 'auc': y_auc, 'auc_ci': y_auc_ci, 'bac': y_bac, 'bac_ci': y_bac_ci, 'f1': y_f1, 'f1_ci': y_f1_ci}

	return metric_dict


## This function learns classification model from training data, then evaluates the learned model on validation data 
def build_simple_classifier(X_train, X_test, y_train, y_test, method, hp_settings):
	## 0. Input argument:
		# X_train: array that contains training feature data 
		# X_test: array that contains testing feature data 
		# y_train: array that contains traning response data 
		# y_test: array that contains testing response data 
		# method: classification method to be used: 'RandomForest', 'XGBoost'
		# hp_settings: string that contains hyperparameter names and specified values (seperated by ',') 

	## 1. Obtain specified value of each hyperparameter   
	st_s = hp_settings.split(',')
	# iterate by hyperparameter 
	ss_dict = {}
	for ss in st_s:
		# obtain the specified value of current hyperparameter 
		ss_s = ss.split(':')
		ss_dict[ss_s[0]] = ss_s[1]
		
	## 2. Learn classification model on training data using the specified hyperparameter values  
	# define random forest classification model, specify hyperparameter 
	if method == 'RandomForest':
		simple_classifier = RandomForestClassifier(random_state = 0, n_estimators = np.int(ss_dict['n_estimators']), criterion = ss_dict['criterion'], max_features = np.float(ss_dict['max_features']), min_samples_split = np.int(ss_dict['min_samples_split']), min_samples_leaf = np.int(ss_dict['min_samples_leaf']), bootstrap = np.bool(ss_dict['bootstrap']))
		hp_char = '_ne_' + ss_dict['n_estimators'] + '_ct_' + ss_dict['criterion'] + '_mf_' + ss_dict['max_features'] + '_ms_' + ss_dict['min_samples_split'] + '_ml_' + ss_dict['min_samples_leaf'] + '_bs_' + ss_dict['bootstrap'] 
	# define gradient boosting classification model, specify hyperparameter
	if method == 'XGBoost':
		simple_classifier = xgb.XGBClassifier(random_state = 0, n_estimators = np.int(ss_dict['n_estimators']), max_depth = np.int(ss_dict['max_depth']), learning_rate = np.float(ss_dict['learning_rate']), subsample = np.float(ss_dict['subsample']), min_child_weight = np.int(ss_dict['min_child_weight']))
		hp_char = '_ne_' + ss_dict['n_estimators'] + '_md_' + ss_dict['max_depth'] + '_lr_' + ss_dict['learning_rate'] + '_ss_' + ss_dict['subsample'] + '_mw_' + ss_dict['min_child_weight']
	# learn classification model on training dataset 
	simple_classifier.fit(X_train, y_train)

	## 3. Evaluate model performance  
	# evaluate performance of learned model on training data 
	train_metrics = evaluate_classifier_by_metrics(simple_classifier, X_train, y_train) 
	# evaluate performance of learned model on validation data
	test_metrics = evaluate_classifier_by_metrics(simple_classifier, X_test, y_test)  

	return simple_classifier, hp_char, train_metrics, test_metrics


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


## This function generates content for output performance file   
def generate_simple_performance_file(N_train_instances, N_test_instances, hp_settings, train_metrics, test_metrics):
	## 0. Input arguments 
		# N_train_instanceis: number of samples in training dataset  
		# N_test_instances: number of samples in validation dataset
		# hp_settings: string that contains the hyperparameter setting of interest 
		# train_metrics: dictionary that contains the model performance on training dataset
		# test_metrics: dictionary that contains the model performance on validation dataset

	## 1. Convert performance metric dictionaries to strings  
	# convert training metric dictionary keys and values to string 
	train_metric_str = convert_dict_to_string(train_metrics)
	# convery validation metric dictionary keys and values to string 
	test_metric_str = convert_dict_to_string(test_metrics)		
	
	## 2. Generate list of strings that describes model performance info
	perf_list = []	
	perf_list.append('Number of training instances: ' + str(N_train_instances))
	perf_list.append('Number of testing instances: ' + str(N_test_instances))
	perf_list.append('Hyperparameter setting: ' + hp_settings)
	perf_list.append('Training performance: ' + train_metric_str)
	perf_list.append('Testing performance: ' + test_metric_str)

	return perf_list	
