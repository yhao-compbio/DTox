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
    + [`simple/interpret_by_lime.py`](simple/interpret_by_lime.py) implements LIME technique to explain sample-level predictions of simple machine learning models. 

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
  + Interpretation result analysis 
     + [`analysis_interpret/compute_hyperparameter_similarity.R`](analysis_interpret/compute_hyperparameter_similarity.R) compares the significant DTox paths detected under different hyperparameter settings of layer-wise relevance propagation rule on Tox21 datasets, and compute Jaccary Index to measure the similarity among distinct hyperparameter settings.
    + [`analysis_interpret/visualize_hyperparameter_similarity.py`](analysis_interpret/visualize_hyperparameter_similarity.py) uses heatmap to visualize the similarity of significant DTox paths under different hyperparameter settings of layer-wise relevance propagation rule on Tox21 datasets.
  + Interpretation validation by gene expression
    + [`analysis_expression/valid_interpret_by_expression.R`](analysis_expression/valid_interpret_by_expression.R) uses LINCS pertubation gene expression data to validate whether significant DTox paths (identified from model interpretation) are differentially expressed after compound treatment, and compare the proportion of differential expression to backtround DTox paths.
    + [`analysis_expression/collect_valid_expression_results.R`](analysis_expression/collect_valid_expression_results.R) collects computed differential expression proportions of compounds from interpretation-validation result files, and performs t test to compare proportion among significant DTox paths vs among background DTox paths.
    + [`analysis_expression/visualize_expression_validation.py`](analysis_expression/visualize_expression_validation.py) uses scatter plots to visualize the gene expression-validation of DTox interpretation results on Tox21 datasets, comparing the compound differential expression proportion compounds among significant DTox paths vs among background DTox paths.
    + [`analysis_expression/analyze_interpret_expression.R`](analysis_expression/analyze_interpret_expression.R) analyzes gene expression-validated DTox paths from model interpretation on Tox21 datasets, and identifies recurrent differentially expressed DTox paths among compounds.
    + [`analysis_expression/visualize_recurrent_path.py`](analysis_expression/visualize_recurrent_path.py) uses barplot to visualize the frequency of recurrent differentially expressed DTox paths from model interpretation results on Tox21 dataset of interest.
  + Interpretation validation by standard pathway-receptor pattern
    + [`analysis_standard/valid_interpret_by_standard.R`](analysis_standard/valid_interpret_by_standard.R) uses standard Reactome pathway-receptor patterns to validate whether significant DTox paths (identified from model interpretation) contains particular pattern matched with each compound, and compare the observed outcome with expected probability. 
    + [`analysis_standard/collect_valid_standard_results.R`](analysis_standard/collect_valid_standard_results.R) collects observed outcome and expected probability of compounds from interpretation-validation result files, then compute observed and expected proportion of validated compounds based on collected results.
    + [`analysis_standard/visualize_standard_validation.py`](analysis_standard/visualize_standard_validation.py)  uses density plot and barplot to visualize the standard pattern-validation of DTox interpretation results on Tox21 datasets, comparing the observed and expected proportion of validated compounds.
    + [`analysis_standard/interpret_by_read_across.R`](analysis_standard/interpret_by_read_across.R) implements Read-across to connect query compounds with query target based on their chemical similarity to source compounds in Drugbank/ComptoxAI.
    + [`analysis_standard/collect_target_standard_results.R`](analysis_standard/collect_target_standard_results.R) collects the validation results by DTox, LIME, and Read-across regarding the interpretation task of connecting active compounds to their respective target receptor in four nuclear receptor assays.
    + [`analysis_standard/visualize_target_standard_validation.py`](analysis_standard/visualize_target_standard_validation.py) uses line charts to visualize the validation performance comparison among DTox, LIME, and Read-across regarding the interpretation task of connecting active compounds to their respective target receptor in four nuclear receptor assays.
  + Interpretation analysis on HepG2 cell viability assay
    + [`analysis_viability/analyze_viability_path_assay.R`](analysis_viability/analyze_viability_path_assay.R) analyzes DTox module relevance scores of viability-related pathways in the context of two viability-related assays (CASP3/7 apoptosis and mitochondria toxicity), compares pathway relevance scores between active and inactive compounds, then uses survival plot to visualize the comparison.
    + [`analysis_viability/analyze_viability_path_map.R`](analysis_viability/analyze_viability_path_map.R) analyzes viability-related DTox paths from model interpretation results in the context of drug-induced liver injury (DILI) adverse events and ATC drug classification, evaluates the enrichment of DILI events/ATC drug classes among compounds identified with viability-related DTox paths, then visualizes the relationships between viability-related DTox paths and DILI events/ATC drug classes by heatmap.
    + [`analysis_viability/analyze_viability_network.R`](analysis_viability/analyze_viability_network.R) uses visNetwork package to visualize the flow of relevance along DTox viability-related paths between query compound, hidden pathway modules, and the HepG2 cell viability outcome.

+ Model prediction 
  + [`analysis_prediction/analyze_prediction_dili.R`](analysis_prediction/analyze_prediction_dili.R) analyzes the DTox HepG2 cell viability model prediction results on DSSTox compounds, compares the predicted outcome probability bewteen positive and negative DSSTox compounds associated with DILI phenotypes.
  + [`analysis_prediction/visualize_prediction_dili.py`](analysis_prediction/visualize_prediction_dili.py) uses boxplot and barplot to visualize the comparison of DTox HepG2 viability prediction results between positive and negative compounds of drug-induced liver injury (DILI) phenotypes.
  + [`analysis_prediction/analyze_prediction_diki.R`](analysis_prediction/analyze_prediction_diki.R) analyzes the DTox HEK293 cell viability model prediction results on DSSTox compounds, compares the predicted outcome probability bewteen positive and negative DSSTox compounds associated with drug-induced kidney injury (DIKI) phenotypes.
  + [`analysis_prediction/visualize_prediction_diki.py`](analysis_prediction/visualize_prediction_diki.py) uses boxplot and barplot to visualize the comparison of DTox HEK293 viability prediction results between positive and negative compounds of drug-induced kidney injury (DIKI) phenotypes. 

+ [`functions.R`](functions.R) contains R functions required for other scripts in the repository.

## Executable shell scripts

+ DTox model implementation 
  + [`run/run_dtox_implementation_compound_target_tox21.sh`](run/run_dtox_implementation_compound_target_tox21.sh) runs [`run/run_dtox_implementation.R`](run/run_dtox_implementation.R) to generate [`run/dtox_compound_target_probability_tox21_implementation.sh`](run/dtox_compound_target_probability_tox21_implementation.sh). [`run/dtox_compound_target_probability_tox21_implementation.sh`](run/dtox_compound_target_probability_tox21_implementation.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.
  + [`run/dtox_compound_target_probability_tox21_shuffle.sh`](run/dtox_compound_target_probability_tox21_shuffle.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under shuffled Reactome pathway hierarchy.
  + [`run/dtox_compound_target_probability_tox21_null.sh`](run/dtox_compound_target_probability_tox21_null.sh) implements [`dtox.py`](dtox.py) on compound target binding-shuffled Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.

+ DTox model interpretation 
  + [`run/interpret_dtox_compound_target_probability_tox21_implementation.sh`](run/interpret_dtox_compound_target_probability_tox21_implementation.sh) implements [`interpret_dtox.py`](`interpret_dtox.py`) on optimal models trained for compound target binding-Tox21 assay outcome datasets.

+ Simple learning implementation
  + [`run/run_simple_compound_target_tox21.sh`](run/run_simple_compound_target_tox21.sh) runs [`run/run_simple.R`](run/run_simple.R) to generate [`run/simple_compound_target_probability_tox21_randomforest.sh`](run/simple_compound_target_probability_tox21_randomforest.sh) and [`run/simple_compound_target_probability_tox21_xgboost.sh`](run/simple_compound_target_probability_tox21_xgboost.sh). [`run/simple_compound_target_probability_tox21_randomforest.sh`](run/simple_compound_target_probability_tox21_randomforest.sh) implements [`simple/simple.py`](simple/simple.py) to build random forest models on compound target binding-Tox21 assay outcome datasets under different hyperparameter settings. [`run/simple_compound_target_probability_tox21_xgboost.sh`](run/simple_compound_target_probability_tox21_xgboost.sh) implements [`simple/simple.py`](simple/simple.py) to build gradient boosting models on compound target binding-Tox21 assay outcome datasets under different hyperparameter settings.
  + [`run/interpret_by_lime.sh`](run/interpret_by_lime.sh) implements [`simple/interpret_by_lime.py`](simple/interpret_by_lime.py) to compute the LIME feature relevance scores of all positive instances in the four nuclear receptor assays

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

+ Model interpretation analysis and comparison
  + Interpretation validation by standard pathway-receptor pattern
    + [`run/interpret_by_read_across.sh`](run/interpret_by_read_across.sh) implements [`analysis_standard/interpret_by_read_across.R`](analysis_standard/interpret_by_read_across.R) to connect query compounds with the four nuclear receptors.  
  + Interpretation validation by gene expression
    + [`run/valid_interpret_by_expression.sh`](run/valid_interpret_by_expression.sh) implements [`analysis_expression/valid_interpret_by_expression.R`](analysis_expression/valid_interpret_by_expression.R) to validate DTox model interpreatation results derived from different LRP rule hyperparameters.
  + Interpretation validation by standard pathway-receptor pattern
    + [`run/valid_interpret_by_standard.sh`](run/valid_interpret_by_standard.sh) implements [`analysis_standard/valid_interpret_by_standard.R`](analysis_standard/valid_interpret_by_standard.R) to validate DTox model interpreatation results derived from different LRP rule hyperparameters.

+ Model prediction
  + [`run/predict_dtox.sh`](run/predict_dtox.sh) implements [`predict_dtox.py`](predict_dtox.py) to predict HepG2 and HEK293 viability outcome using the respective trained optimal DTOx model.
