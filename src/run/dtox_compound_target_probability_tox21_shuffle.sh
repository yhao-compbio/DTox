#!/bin/bash
#BSUB -q epistasis_normal
#BSUB -J dtox_compound_target_probability_tox21_shuffle
#BSUB -n 15
#BSUB -o src/run/dtox_compound_target_probability_tox21_shuffle.%J.out
#BSUB -e src/run/dtox_compound_target_probability_tox21_shuffle.%J.error
#BSUB -N

module unload python
module load python/3.7

./src/run/dtox_compound_target_probability_tox21_shuffle1.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle2.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle3.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle4.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle5.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle6.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle7.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle8.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle9.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle10.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle11.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle12.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle13.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle14.sh &
./src/run/dtox_compound_target_probability_tox21_shuffle15.sh &
