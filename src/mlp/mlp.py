# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script develops and evaluates a fully connected Multi-Layer Perceptron (MLP) neural network model, otherwise with the same number of hidden layer/neuron as the matched DTox model


## Module
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'src/mlp/')
import mlp_learning


## Main function 
def main(argv):
	## 0. Input arguments 
		# argv 1: input file that contains combined training-testing data 
		# argv 2: input file that contains validation data  
		# argv 3: name of column that contains label/response data 
		# argv 4: input file that contains number of hidden neurons in each layer (model info/performance file of comparable DTox model)   
		# argv 5: name of output file 

	## 1. Read in input training and testing files
	# read in combined training-testing dataset, separate feature and response data  
	outcome_col = argv[3]
	train_data_df = pd.read_csv(argv[1], sep = '\t', header = 0, index_col = 0)
	train_X_data, train_y_data = train_data_df.drop(outcome_col, axis = 1).values, train_data_df[outcome_col].values
	# read in validation dataset, separate feature and response data   
	test_data_df = pd.read_csv(argv[2], sep = '\t', header = 0, index_col = 0)
	test_X_data, test_y_data = test_data_df.drop(outcome_col, axis = 1).values, test_data_df[outcome_col].values
	
	## 2. Develop and evaluate MLP model  
	# read in model info/performance file, go through file line by line until the desired line is reached   
	size_file = open(argv[4], 'r')
	size_lines = size_file.readlines()
	for sl in size_lines:
		sl =  sl.split('\n')[0]
		# check whether the current line stores number of hidden neurons in each layer 
		if sl.startswith('Number of hidden neurons: '):
			# if so, obtain the numbers, then break iteration
			sl_s = sl.split(': ')[1].split(',')
			neuron_size = [int(ss) for ss in sl_s]
			break
	# learn MLP model from combined training-testing data, then evaluate the learned model on training-testing as well as validation data 
	model, loss_df, training_perf, testing_perf = mlp_learning.train_mlp_model(train_X_data, test_X_data, train_y_data, test_y_data, neuron_size, model_name = argv[5] + '_model.pt')
	loss_df.to_csv(argv[5] + '_loss.tsv', sep = '\t', index = False, float_format = '%.5f')
	# generate content for output performance file 
	output_perf_list = mlp_learning.generate_mlp_performance_file(neuron_size, train_X_data, test_X_data, training_perf, testing_perf)
	# write generated model performance info to output file  
	perf_op = open(argv[5] + '_performance.txt', 'w')
	for opl in output_perf_list:
		perf_op.write('%s\n' % opl)
	perf_op.close()

	return 1


## Call main function
main(sys.argv)
