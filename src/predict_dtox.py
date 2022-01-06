# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script implements trained DTox model to predict outcome probability based on input feature data


## Module
import sys
import numpy as np
import pandas as pd
import torch
sys.path.insert(0, 'src/')
import dtox_hierarchy
import dtox_nn
import dtox_data


## Main function
def main(argv):
	## 0. Input arguments 
		# argv 1: input file that contains root pathways of DTox model 
		# argv 2: input file that contains parent/children node connections of DTox model  
		# argv 3: input file that contains number of annotated genes in nodes of DTox model 
		# argv 4: input file that contains layer number of nodes in DTox model 
		# argv 5: minimal size of pathways included in DTox model
		# argv 6: maximal size of node modules in DTox neural network model
		# argv 7: input file that contains trained DTox model 
		# argv 8: input file that contains feature data for prediction 
		# argv 9: name of output file 

	## 1. Load trained DTox model 
	# sort the hierarchy of DTox model based on input hierarchy files
	dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_hierarchy_stat, dtox_input_size, dtox_hidden_size = dtox_hierarchy.sort_dtox_hierarchy(argv[1], argv[2], argv[3], argv[4], int(argv[5]), int(argv[6])) 
	# define structure of whole neural network  
	dtox_model = dtox_nn.DToxNet(dtox_node_children, dtox_in_module_size, dtox_out_module_size, dtox_root, dtox_out_root_size, dtox_input_size) 
	# load trained DTox model parameter into defined structure  
	model_state = torch.load(argv[7])
	dtox_model.load_state_dict(model_state['model_state_dict'])
	
	## 2. Implement trained  
	# read in feature data for prediction, convert into format for DTox model  
	predict_data_df = pd.read_csv(argv[8], sep = '\t', header = 0, index_col = 0)
	pred_data = dtox_data.DTox_dataformat(predict_data_df.values, np.zeros(predict_data_df.shape[0]))
	# perform forward propagation to predict outcome probability   
	y_pred, _ = dtox_model(pred_data[: , :][0]) 
	# write predicted outcome probability with feature data to output file  
	predict_data_df['predicted_outcome'] = y_pred.data.numpy()	
	predict_data_df.to_csv(argv[9], sep = '\t')

	return 1


## Call main function
main(sys.argv)
