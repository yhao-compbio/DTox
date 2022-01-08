#!/bin/bash
#BSUB -q epistasis_normal
#BSUB -J mlp_compound_target_probability_tox21_fully_connected
#BSUB -n 15
#BSUB -o src/run/mlp_compound_target_probability_tox21_fully_connected.%J.out
#BSUB -e src/run/mlp_compound_target_probability_tox21_fully_connected.%J.error
#BSUB -N

module unload python
module load python/3.7

./src/run/mlp_compound_target_probability_tox21_fully_connected1.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected2.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected3.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected4.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected5.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected6.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected7.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected8.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected9.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected10.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected11.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected12.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected13.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected14.sh &
./src/run/mlp_compound_target_probability_tox21_fully_connected15.sh &
