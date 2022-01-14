# This folder contains source code used by the repository.

## R/python scripts 

+ DTox model 
  + [`dtox.py`](dtox.py) learns and evaluates DTox model
    + [`dtox_data.py`](dtox_data.py) contains data-formatting functions used for DTox model training.
    + [`dtox_hierarchy.py`](dtox_hierarchy.py) contains functions used to process sorted DTox hiearchy files and compute model statistics.
    + [`dtox_nn.py`](dtox_nn.py) contains functions used to build basic neural network structure for DTox model.
    + [`dtox_loss.py`](dtox_loss.py) contains the the loss function used in DTox model.
    + [`early_stop.py`](early_stop.py) contains early stop function of DTox model.
    + [`dtox_learning.py`](dtox_learning.py) contains deep learning functions used in the DTox model construction.
    + [`run/run_dtox_implementation.R`](run/run_dtox_implementation.R) generates shell scripts that run DTox model on Tox21 datasets under Reactome pathway hierarchy. [`run/run_dtox_shuffle.R`](run/run_dtox_shuffle.R) generates scripts that run DTox model on Tox21 datasets under shuffeld Reactome pathway hierarchy. [`run/run_dtox_null.R`](run/run_dtox_null.R) generates shell scripts that run DTox model on outcome-shuffled Tox21 datasets under Reactome pathway hierarchy.
  + [`predict_dtox.py`](predict_dtox.py) implements trained DTox model to predict outcome probability based on input feature data.
  + [`interpret_dtox.py`](`interpret_dtox.py`) implements layer-wise relevance propagation to evaluate relevance of DTox paths.
    + [`dtox_lrp.py`](dtox_lrp.py) contains functions used for implementing LRP to evaluate relevance of DTox paths.
    + [`run/run_interpret_dtox.R`](run/run_interpret_dtox.R) generates shell scripts that runs DTox interpretation procedure on optimal models trained for Tox21 datasets.

+ Simple machine learning model
  + [`simple/simple.py`](simple/simple.py) develops and evaluates simple machine learning model (random forest or gradient boosting).
    + [`simple/simple_learning.py`](simple/simple_learning.py) contains functions for building, evaluating, and implementing simple machine learning models.
    + [`run/run_simple.R`](run/run_simple.R) generates shell scripts that run simple machine learning models on Tox21 datasets under different hyperparameter settings.

+ Multi-layer perceptron neural network model
  + [`mlp/mlp.py`](mlp/mlp.py) develops and evaluates a fully connected Multi-Layer Perceptron (MLP) neural network model, otherwise with the same number of hidden layer/neuron as the matched DTox model. 
    + [`mlp/mlp_learning.py`](mlp/mlp_learning.py) contains functions used in the Multi-Layer Perceptron (MLP) neural network model.
    + [`run/run_mlp.R`](run/run_mlp.R) generates shell scripts that run fully connected MLP neural network models on Tox21 datasets, which are built with the same number of hidden layer/neuron as matched DTox models.

+ Model performance analysis, comparison, and visualization 
  + [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) collects machine learning model basic info and performance metrics from performance files.
  + [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) identifies optimal hyperparameter setting of machine learning method implementation, then compares and visualizes model performance across different method implementations.
    + [`analysis_dtox/dtox_analysis.py`](analysis_dtox/dtox_analysis.py) contains functions used in DTox model result anaysis.
    + [`analysis_dtox/dtox_plot.py`](analysis_dtox/dtox_plot.py) contains functions for visualizing DTox model results.
  + [`analysis_dtox/compare_dtox_hyperparameter.R`](analysis_dtox/compare_dtox_hyperparameter.R) normalizes DTox model performance across a query hyperparameter for comparison.
  + [`analysis_dtox/visualize_hyperparameter_comparison.py`](analysis_dtox/visualize_hyperparameter_comparison.py) uses heatmap and upsetplot to visualize normalized model performance of Tox21 datasets across root pathway settings.
  + [`analysis_dtox/compute_dtox_connections.R`](analysis_dtox/compute_dtox_connections.R) computes total number of parameters in DTox and matched fully connected multi-layer perceptron (MLP) models.
  + [`analysis_dtox/visualize_parameter_comparison.py`](analysis_dtox/visualize_parameter_comparison.py) uses barplot to visualize comparison of DTox and MLP model statistics across Tox21 datasets.
  + [`analysis_dtox/visualize_training_loss.py`](analysis_dtox/visualize_training_loss.py) uses line charts to visualize evolution of training/testing loss over epoches during DTox learning process.

+ Model interpretation analysis and comparison
  + [`analysis_interpret/compute_hyperparameter_similarity.R`](analysis_interpret/compute_hyperparameter_similarity.R) compares the significant DTox paths detected under different hyperparameter settings of layer-wise relevance propagation rule on Tox21 datasets, and compute Jaccary Index to measure the similarity among distinct hyperparameter settings.
  + [`analysis_interpret/visualize_hyperparameter_similarity.py`](analysis_interpret/visualize_hyperparameter_similarity.py) uses heatmap to visualize the similarity of significant DTox paths under different hyperparameter settings of layer-wise relevance propagation rule on Tox21 datasets.

## Executable shell scripts

+ DTox model implementation 
  + [`run/run_dtox_implementation_compound_target_tox21.sh`](run/run_dtox_implementation_compound_target_tox21.sh) runs [`run/run_dtox_implementation.R`](run/run_dtox_implementation.R) to generate [`run/dtox_compound_target_probability_tox21_implementation.sh`](run/dtox_compound_target_probability_tox21_implementation.sh). [`run/dtox_compound_target_probability_tox21_implementation.sh`](run/dtox_compound_target_probability_tox21_implementation.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.
  + [`run/dtox_compound_target_probability_tox21_shuffle.sh`](run/dtox_compound_target_probability_tox21_shuffle.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under shuffled Reactome pathway hierarchy.
  + [`run/dtox_compound_target_probability_tox21_null.sh`](run/dtox_compound_target_probability_tox21_null.sh) implements [`dtox.py`](dtox.py) on compound target binding-shuffled Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.

+ DTox model interpretation 
  + [`run/interpret_dtox_compound_target_probability_tox21_implementation.sh`](run/interpret_dtox_compound_target_probability_tox21_implementation.sh) implements [`interpret_dtox.py`](`interpret_dtox.py`) on optimal models trained for compound target binding-Tox21 assay outcome datasets.

+ Simple learning implementation
  + [`run/run_simple_compound_target_tox21.sh`](run/run_simple_compound_target_tox21.sh) runs [`run/run_simple.R`](run/run_simple.R) to generate [`run/simple_compound_target_probability_tox21_randomforest.sh`](run/simple_compound_target_probability_tox21_randomforest.sh) and [`run/simple_compound_target_probability_tox21_xgboost.sh`](run/simple_compound_target_probability_tox21_xgboost.sh). [`run/simple_compound_target_probability_tox21_randomforest.sh`](run/simple_compound_target_probability_tox21_randomforest.sh) implements [`simple/simple.py`](simple/simple.py) to build random forest models on compound target binding-Tox21 assay outcome datasets under different hyperparameter settings. [`run/simple_compound_target_probability_tox21_xgboost.sh`](run/simple_compound_target_probability_tox21_xgboost.sh) implements [`simple/simple.py`](simple/simple.py) to build gradient boosting models on compound target binding-Tox21 assay outcome datasets under different hyperparameter settings.

+ Multi-layer perceptron neural network implementation 
  + [`run/run_mlp.sh`](run/run_mlp.sh) runs [`run/run_mlp.R`](run/run_mlp.R) to generate [`run/mlp_compound_target_probability_tox21_fully_connected.sh`](run/mlp_compound_target_probability_tox21_fully_connected.sh). [`run/mlp_compound_target_probability_tox21_fully_connected.sh`](run/mlp_compound_target_probability_tox21_fully_connected.sh) implements [`mlp/mlp.py`](mlp/mlp.py) on compound target binding-Tox21 assay outcome datasets. 

+ Model performance analysis, comparison, and visualization
  + Result collection 
    + [`run/collect_model_results_compound_target_tox21_implementation.sh`](run/collect_model_results_compound_target_tox21_implementation.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of DTox models built upon compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy. 
    + [`run/collect_model_results_compound_target_tox21_shuffle.sh`](run/collect_model_results_compound_target_tox21_shuffle.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of DTox models built upon compound target binding-Tox21 assay outcome datasets under shuffled Reactome pathway hierarchy.
    + [`run/collect_model_results_compound_target_tox21_simple.sh`](run/collect_model_results_compound_target_tox21_simple.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of simple machine learning models built upon compound target binding-Tox21 assay outcome datasets under different hyperparameter settings.
    + [`run/collect_model_results_compound_target_tox21_mlp.sh`](run/collect_model_results_compound_target_tox21_mlp.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of MLP modles built upon compound target binding-Tox21 assay outcome datasets.
  + Result analysis 
    + [`run/analyze_dtox_results_compound_target_tox21_simple.sh`](run/analyze_dtox_results_compound_target_tox21_simple.sh) implements [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) to identify optimal hyperparameter setting of simple machine learning model implementation on compound target binding-Tox21 assay outcome datasets. 
    + [`run/analyze_dtox_results_compound_target_tox21_dtox.sh`](run/analyze_dtox_results_compound_target_tox21_dtox.sh) implements [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) to identify optimal hyperparameter setting of DTox model implementation on compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy, then compare and visualize model performance across different method implementations.
