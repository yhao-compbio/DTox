# !/usr/iinfenv python
## created by Yun Hao @MooreLab 2021
## This script develops and evaluates simple machine learning model (random forest or gradient boosting) 


# Module
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'src/simple/')
import simple_learning


## Main function 
def main(argv):
	## 0. Input arguments 
		# argv 1: input file that contains training data 
		# argv 2: input file that contains validation data 
		# argv 3: name of column that contains label/response data 
		# argv 4: name of output file  
		# argv 5: name of classification method to be used: 'RandomForest', 'XGBoost'
		# argv 6: string that contains hyperparameter names and specified values (seperated by ',')  

	## 1. Process training and validation datasets 
	# read in training dataset, separate feature and response data 
	outcome_col = argv[3]
	train_data_df = pd.read_csv(argv[1], sep = '\t', header = 0, index_col = 0)
	train_X_data, train_y_data = train_data_df.drop(outcome_col, axis = 1).values, train_data_df[outcome_col].values
	# read in validation dataset, separate feature and response data  
	test_data_df = pd.read_csv(argv[2], sep = '\t', header = 0, index_col = 0)
	test_X_data, test_y_data = test_data_df.drop(outcome_col, axis = 1).values, test_data_df[outcome_col].values
	
	## 2. Develop and evaluate simple machine learning model 
	# learn classification model from training data, then evaluate the learned model on training as well as validation data 
	simple_learner, hyper_str, train_perf, test_perf = simple_learning.build_simple_classifier(train_X_data, test_X_data, train_y_data, test_y_data, argv[5], argv[6])
	# generate content for output performance file 
	output_perf_list = simple_learning.generate_simple_performance_file(train_X_data.shape[0], test_X_data.shape[0], argv[6], train_perf, test_perf)
	# write generated model performance info to output file 
	output_file = argv[4] + '_md_' + argv[5] + hyper_str + '_performance.txt'
	perf_op = open(output_file, 'w')
	for opl in output_perf_list:
		perf_op.write('%s\n' % opl)
	perf_op.close()

	return 1	


## Call main function
main(sys.argv)
