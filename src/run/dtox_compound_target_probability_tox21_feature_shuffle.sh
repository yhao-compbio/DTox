#!/bin/bash
#BSUB -q epistasis_long
#BSUB -J dtox_compound_target_probability_tox21_feature_shuffle
#BSUB -n 3
#BSUB -o src/run/dtox_compound_target_probability_tox21_feature_shuffle.%J.out
#BSUB -e src/run/dtox_compound_target_probability_tox21_feature_shuffle.%J.error
#BSUB -N

module unload python
module load python/3.7

./src/run/dtox_compound_target_probability_tox21_feature_shuffle1.sh &
./src/run/dtox_compound_target_probability_tox21_feature_shuffle2.sh &
./src/run/dtox_compound_target_probability_tox21_feature_shuffle3.sh &
