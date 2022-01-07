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

## Executable shell scripts

+ DTox model implementation 
  + [`run/run_dtox_implementation_compound_target_tox21.sh`](run/run_dtox_implementation_compound_target_tox21.sh) runs [`run/run_dtox_implementation.R`](run/run_dtox_implementation.R) to generate [`run/dtox_compound_target_probability_tox21_implementation.sh`](run/dtox_compound_target_probability_tox21_implementation.sh). [`run/dtox_compound_target_probability_tox21_implementation.sh`](run/dtox_compound_target_probability_tox21_implementation.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.
  + [`run/dtox_compound_target_probability_tox21_shuffle.sh`](run/dtox_compound_target_probability_tox21_shuffle.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under shuffled Reactome pathway hierarchy.
  + [`run/dtox_compound_target_probability_tox21_null.sh`](run/dtox_compound_target_probability_tox21_null.sh) implements [`dtox.py`](dtox.py) on compound target binding-shuffled Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.

+ DTox model interpretation 
  + [`run/interpret_dtox_compound_target_probability_tox21_implementation.sh`](run/interpret_dtox_compound_target_probability_tox21_implementation.sh) implements [`interpret_dtox.py`](`interpret_dtox.py`) on optimal models trained for compound target binding-Tox21 assay outcome datasets.

+ Simple learning implementation
  + [`run/run_simple_compound_target_tox21.sh`](run/run_simple_compound_target_tox21.sh) runs [`run/run_simple.R`](run/run_simple.R) to generate [`run/simple_compound_target_probability_tox21_randomforest.sh`](run/simple_compound_target_probability_tox21_randomforest.sh) and [`run/simple_compound_target_probability_tox21_xgboost.sh`](run/simple_compound_target_probability_tox21_xgboost.sh). [`run/simple_compound_target_probability_tox21_randomforest.sh`](run/simple_compound_target_probability_tox21_randomforest.sh) implements [`simple/simple.py`](simple/simple.py) to build random forest models on compound target binding-Tox21 assay outcome datasets under different hyperparameter settings. [`run/simple_compound_target_probability_tox21_xgboost.sh`](run/simple_compound_target_probability_tox21_xgboost.sh) implements [`simple/simple.py`](simple/simple.py) to build gradient boosting models on compound target binding-Tox21 assay outcome datasets under different hyperparameter settings. 
