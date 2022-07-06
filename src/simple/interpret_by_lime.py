# !/usr/iinfenv python
## created by Yun Hao @MooreLab 2022
## This script implements LIME technique to explain sample-level predictions of simple machine learning models 


## Module
import sys
import numpy as np
import pandas as pd
from lime import lime_tabular
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
		# argv 7: name of feature of interest  

	## 1. Process training and validation datasets 
	# read in training dataset, separate feature and response data 
	outcome_col = argv[3]
	train_data_df = pd.read_csv(argv[1], sep = '\t', header = 0, index_col = 0)
	train_X_data, train_y_data = train_data_df.drop(outcome_col, axis = 1).values, train_data_df[outcome_col].values
	train_features = train_data_df.drop(outcome_col, axis = 1).columns
	# read in validation dataset, separate feature and response data  
	test_data_df = pd.read_csv(argv[2], sep = '\t', header = 0, index_col = 0)
	test_X_data, test_y_data = test_data_df.drop(outcome_col, axis = 1).values, test_data_df[outcome_col].values
	# combine training and validation set, extract instances of the positive class 
	combine_data_df = pd.concat([train_data_df, test_data_df], axis = 0)
	pos_data_df = combine_data_df[combine_data_df.assay_outcome == 1]
	pos_feature_data = pos_data_df.drop(outcome_col, axis = 1).values
	
	## 2. Implement LIME technique to explain sample-level predictions made by specified machine learning models   
	# learn the specified classification model from training data 
	simple_learner, hyper_str, train_perf, test_perf = simple_learning.build_simple_classifier(train_X_data, test_X_data, train_y_data, test_y_data, argv[5], argv[6])
	# implement the learned model to predict the probability of positive class for all positive instances   
	pos_data_pred = simple_learner.predict_proba(pos_feature_data)[:,1]
	avg_impt = pos_data_pred/len(train_features)
	# build LIME explainer  
	explainer = lime_tabular.LimeTabularExplainer(train_X_data, feature_names = train_features)
	# iterate by instance
	target_range = []
	target_impt = []
	for i in range(0, len(pos_feature_data)): 
		# implement LIME explainer to compute the relevance score of each feature  
		i_instance = pos_feature_data[i]
		i_exp = explainer.explain_instance(i_instance, simple_learner.predict_proba, num_features = len(train_features), num_samples = 50)
		# extract the relevance score of feature of interest  
		i_exp_list = i_exp.as_list()
		for iel in range(0, len(i_exp_list)):
			if argv[7] in i_exp_list[iel][0]:
				target_range.append(i_exp_list[iel][0])
				target_impt.append(i_exp_list[iel][1])
				break
	# output the computed feature relevance scores of all positive instances, along with other info 
	exp_result_df = pd.DataFrame({'compound_cid':pos_data_df.index, 'compound_predicted_score': pos_data_pred, 'target_range': target_range, 'target_relevance': target_impt, 'average_relevance': avg_impt})
	output_file = argv[4] + '_interpret_by_lime_result.tsv'	
	exp_result_df.to_csv(output_file, sep = '\t', index = False)

	return 1


## Call main function
main(sys.argv)

