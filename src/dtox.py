# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script learns and evaluates DTox model


## Module
import sys
import torch
import numpy as np
import pandas as pd
sys.path.insert(0, 'src/')
import dtox_data
import dtox_learning


# Main function
def main(argv):
	## 0. Input arguments 
		# argv 1: input file that contains combined training and testing data 
		# argv 2: input file that contains validation data 
		# argv 3: name of column that contains label/response data
		# argv 4: input file that contains root pathways of DTox model 
		# argv 5: input file that contains parent/children node connections of DTox model  
		# argv 6: input file that contains number of annotated genes in nodes of DTox model 
		# argv 7: input file that contains layer number of nodes in DTox model 
		# argv 8: minimal size of pathways included in DTox model
		# argv 9: name of output file  
		# argv 10: maximal size of node modules in DTox neural network model
		# argv 11: coefficient for auxiliary loss in loss function (see 'dtox_loss')
		# argv 12: coefficient for L2 regularization 

	## 1. Read in and format input trainig, testing, and validation data 
	combined_data, training_data_loader, testing_data, validation_data = dtox_data.format_dtox_data(argv[1], argv[2], argv[3])

	## 2. Learn DTox model 
	# add specified hyperparmeters into name string of output files 
	output_file = argv[9] + '_xs_' + argv[10] + '_al_' + argv[11] + '_ld_' + argv[12] 
	output_model_name = output_file + '_model.pt'
	# check whether GPU training is avaiable, if not use CPU training
	device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
	# learn DTox model with training and testing data (training data for learning parameters, testing data for implementing early stop), save model to output model file 
	torch.manual_seed(0)
	hierarchy_info, trained_model, loss, training_summary_df = dtox_learning.train_dtox_model(argv[4], argv[5], argv[6], argv[7], int(argv[8]), int(argv[10]), training_data_loader, testing_data, float(argv[11]), float(argv[12]), dtox_model_name = output_model_name, dtox_device = device)
	# write tracked training and testing loss to output loss file
	training_summary_df.to_csv(output_file + '_loss.tsv', sep = '\t', index = False, float_format = '%.5f')
	
	## 3. Evaluate DTox model
	# evaluate learned DTox model on combined training and testing data (training performance) 
	combined_perf = dtox_learning.evaluate_dtox_model(trained_model, loss, combined_data, dtox_device = device)
	# evaluate learned DTox model on validation data (validation performance)
	validation_perf = dtox_learning.evaluate_dtox_model(trained_model, loss, validation_data, dtox_device = device)
	# generate model performance strings, write to output performance file  
	perf_info = dtox_learning.generate_dtox_performance_file(combined_data, combined_perf, validation_data, validation_perf, hierarchy_info)
	perf_op = open(output_file + '_performance.txt', 'w')
	for pi in perf_info:
		perf_op.write('%s\n' % pi)
	perf_op.close()

	return 1

## Call main function
main(sys.argv)

