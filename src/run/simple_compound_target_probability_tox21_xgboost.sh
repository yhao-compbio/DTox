#!/bin/bash
#BSUB -q epistasis_long
#BSUB -J simple_compound_target_probability_tox21_xgboost
#BSUB -n 200
#BSUB -o src/run/simple_compound_target_probability_tox21_xgboost.%J.out
#BSUB -e src/run/simple_compound_target_probability_tox21_xgboost.%J.error
#BSUB -N

module unload python
module load python/3.7

./src/run/simple_compound_target_probability_tox21_xgboost1.sh &
./src/run/simple_compound_target_probability_tox21_xgboost2.sh &
./src/run/simple_compound_target_probability_tox21_xgboost3.sh &
./src/run/simple_compound_target_probability_tox21_xgboost4.sh &
./src/run/simple_compound_target_probability_tox21_xgboost5.sh &
./src/run/simple_compound_target_probability_tox21_xgboost6.sh &
./src/run/simple_compound_target_probability_tox21_xgboost7.sh &
./src/run/simple_compound_target_probability_tox21_xgboost8.sh &
./src/run/simple_compound_target_probability_tox21_xgboost9.sh &
./src/run/simple_compound_target_probability_tox21_xgboost10.sh &
./src/run/simple_compound_target_probability_tox21_xgboost11.sh &
./src/run/simple_compound_target_probability_tox21_xgboost12.sh &
./src/run/simple_compound_target_probability_tox21_xgboost13.sh &
./src/run/simple_compound_target_probability_tox21_xgboost14.sh &
./src/run/simple_compound_target_probability_tox21_xgboost15.sh &
./src/run/simple_compound_target_probability_tox21_xgboost16.sh &
./src/run/simple_compound_target_probability_tox21_xgboost17.sh &
./src/run/simple_compound_target_probability_tox21_xgboost18.sh &
./src/run/simple_compound_target_probability_tox21_xgboost19.sh &
./src/run/simple_compound_target_probability_tox21_xgboost20.sh &
./src/run/simple_compound_target_probability_tox21_xgboost21.sh &
./src/run/simple_compound_target_probability_tox21_xgboost22.sh &
./src/run/simple_compound_target_probability_tox21_xgboost23.sh &
./src/run/simple_compound_target_probability_tox21_xgboost24.sh &
./src/run/simple_compound_target_probability_tox21_xgboost25.sh &
./src/run/simple_compound_target_probability_tox21_xgboost26.sh &
./src/run/simple_compound_target_probability_tox21_xgboost27.sh &
./src/run/simple_compound_target_probability_tox21_xgboost28.sh &
./src/run/simple_compound_target_probability_tox21_xgboost29.sh &
./src/run/simple_compound_target_probability_tox21_xgboost30.sh &
./src/run/simple_compound_target_probability_tox21_xgboost31.sh &
./src/run/simple_compound_target_probability_tox21_xgboost32.sh &
./src/run/simple_compound_target_probability_tox21_xgboost33.sh &
./src/run/simple_compound_target_probability_tox21_xgboost34.sh &
./src/run/simple_compound_target_probability_tox21_xgboost35.sh &
./src/run/simple_compound_target_probability_tox21_xgboost36.sh &
./src/run/simple_compound_target_probability_tox21_xgboost37.sh &
./src/run/simple_compound_target_probability_tox21_xgboost38.sh &
./src/run/simple_compound_target_probability_tox21_xgboost39.sh &
./src/run/simple_compound_target_probability_tox21_xgboost40.sh &
./src/run/simple_compound_target_probability_tox21_xgboost41.sh &
./src/run/simple_compound_target_probability_tox21_xgboost42.sh &
./src/run/simple_compound_target_probability_tox21_xgboost43.sh &
./src/run/simple_compound_target_probability_tox21_xgboost44.sh &
./src/run/simple_compound_target_probability_tox21_xgboost45.sh &
./src/run/simple_compound_target_probability_tox21_xgboost46.sh &
./src/run/simple_compound_target_probability_tox21_xgboost47.sh &
./src/run/simple_compound_target_probability_tox21_xgboost48.sh &
./src/run/simple_compound_target_probability_tox21_xgboost49.sh &
./src/run/simple_compound_target_probability_tox21_xgboost50.sh &
./src/run/simple_compound_target_probability_tox21_xgboost51.sh &
./src/run/simple_compound_target_probability_tox21_xgboost52.sh &
./src/run/simple_compound_target_probability_tox21_xgboost53.sh &
./src/run/simple_compound_target_probability_tox21_xgboost54.sh &
./src/run/simple_compound_target_probability_tox21_xgboost55.sh &
./src/run/simple_compound_target_probability_tox21_xgboost56.sh &
./src/run/simple_compound_target_probability_tox21_xgboost57.sh &
./src/run/simple_compound_target_probability_tox21_xgboost58.sh &
./src/run/simple_compound_target_probability_tox21_xgboost59.sh &
./src/run/simple_compound_target_probability_tox21_xgboost60.sh &
./src/run/simple_compound_target_probability_tox21_xgboost61.sh &
./src/run/simple_compound_target_probability_tox21_xgboost62.sh &
./src/run/simple_compound_target_probability_tox21_xgboost63.sh &
./src/run/simple_compound_target_probability_tox21_xgboost64.sh &
./src/run/simple_compound_target_probability_tox21_xgboost65.sh &
./src/run/simple_compound_target_probability_tox21_xgboost66.sh &
./src/run/simple_compound_target_probability_tox21_xgboost67.sh &
./src/run/simple_compound_target_probability_tox21_xgboost68.sh &
./src/run/simple_compound_target_probability_tox21_xgboost69.sh &
./src/run/simple_compound_target_probability_tox21_xgboost70.sh &
./src/run/simple_compound_target_probability_tox21_xgboost71.sh &
./src/run/simple_compound_target_probability_tox21_xgboost72.sh &
./src/run/simple_compound_target_probability_tox21_xgboost73.sh &
./src/run/simple_compound_target_probability_tox21_xgboost74.sh &
./src/run/simple_compound_target_probability_tox21_xgboost75.sh &
./src/run/simple_compound_target_probability_tox21_xgboost76.sh &
./src/run/simple_compound_target_probability_tox21_xgboost77.sh &
./src/run/simple_compound_target_probability_tox21_xgboost78.sh &
./src/run/simple_compound_target_probability_tox21_xgboost79.sh &
./src/run/simple_compound_target_probability_tox21_xgboost80.sh &
./src/run/simple_compound_target_probability_tox21_xgboost81.sh &
./src/run/simple_compound_target_probability_tox21_xgboost82.sh &
./src/run/simple_compound_target_probability_tox21_xgboost83.sh &
./src/run/simple_compound_target_probability_tox21_xgboost84.sh &
./src/run/simple_compound_target_probability_tox21_xgboost85.sh &
./src/run/simple_compound_target_probability_tox21_xgboost86.sh &
./src/run/simple_compound_target_probability_tox21_xgboost87.sh &
./src/run/simple_compound_target_probability_tox21_xgboost88.sh &
./src/run/simple_compound_target_probability_tox21_xgboost89.sh &
./src/run/simple_compound_target_probability_tox21_xgboost90.sh &
./src/run/simple_compound_target_probability_tox21_xgboost91.sh &
./src/run/simple_compound_target_probability_tox21_xgboost92.sh &
./src/run/simple_compound_target_probability_tox21_xgboost93.sh &
./src/run/simple_compound_target_probability_tox21_xgboost94.sh &
./src/run/simple_compound_target_probability_tox21_xgboost95.sh &
./src/run/simple_compound_target_probability_tox21_xgboost96.sh &
./src/run/simple_compound_target_probability_tox21_xgboost97.sh &
./src/run/simple_compound_target_probability_tox21_xgboost98.sh &
./src/run/simple_compound_target_probability_tox21_xgboost99.sh &
./src/run/simple_compound_target_probability_tox21_xgboost100.sh &
./src/run/simple_compound_target_probability_tox21_xgboost101.sh &
./src/run/simple_compound_target_probability_tox21_xgboost102.sh &
./src/run/simple_compound_target_probability_tox21_xgboost103.sh &
./src/run/simple_compound_target_probability_tox21_xgboost104.sh &
./src/run/simple_compound_target_probability_tox21_xgboost105.sh &
./src/run/simple_compound_target_probability_tox21_xgboost106.sh &
./src/run/simple_compound_target_probability_tox21_xgboost107.sh &
./src/run/simple_compound_target_probability_tox21_xgboost108.sh &
./src/run/simple_compound_target_probability_tox21_xgboost109.sh &
./src/run/simple_compound_target_probability_tox21_xgboost110.sh &
./src/run/simple_compound_target_probability_tox21_xgboost111.sh &
./src/run/simple_compound_target_probability_tox21_xgboost112.sh &
./src/run/simple_compound_target_probability_tox21_xgboost113.sh &
./src/run/simple_compound_target_probability_tox21_xgboost114.sh &
./src/run/simple_compound_target_probability_tox21_xgboost115.sh &
./src/run/simple_compound_target_probability_tox21_xgboost116.sh &
./src/run/simple_compound_target_probability_tox21_xgboost117.sh &
./src/run/simple_compound_target_probability_tox21_xgboost118.sh &
./src/run/simple_compound_target_probability_tox21_xgboost119.sh &
./src/run/simple_compound_target_probability_tox21_xgboost120.sh &
./src/run/simple_compound_target_probability_tox21_xgboost121.sh &
./src/run/simple_compound_target_probability_tox21_xgboost122.sh &
./src/run/simple_compound_target_probability_tox21_xgboost123.sh &
./src/run/simple_compound_target_probability_tox21_xgboost124.sh &
./src/run/simple_compound_target_probability_tox21_xgboost125.sh &
./src/run/simple_compound_target_probability_tox21_xgboost126.sh &
./src/run/simple_compound_target_probability_tox21_xgboost127.sh &
./src/run/simple_compound_target_probability_tox21_xgboost128.sh &
./src/run/simple_compound_target_probability_tox21_xgboost129.sh &
./src/run/simple_compound_target_probability_tox21_xgboost130.sh &
./src/run/simple_compound_target_probability_tox21_xgboost131.sh &
./src/run/simple_compound_target_probability_tox21_xgboost132.sh &
./src/run/simple_compound_target_probability_tox21_xgboost133.sh &
./src/run/simple_compound_target_probability_tox21_xgboost134.sh &
./src/run/simple_compound_target_probability_tox21_xgboost135.sh &
./src/run/simple_compound_target_probability_tox21_xgboost136.sh &
./src/run/simple_compound_target_probability_tox21_xgboost137.sh &
./src/run/simple_compound_target_probability_tox21_xgboost138.sh &
./src/run/simple_compound_target_probability_tox21_xgboost139.sh &
./src/run/simple_compound_target_probability_tox21_xgboost140.sh &
./src/run/simple_compound_target_probability_tox21_xgboost141.sh &
./src/run/simple_compound_target_probability_tox21_xgboost142.sh &
./src/run/simple_compound_target_probability_tox21_xgboost143.sh &
./src/run/simple_compound_target_probability_tox21_xgboost144.sh &
./src/run/simple_compound_target_probability_tox21_xgboost145.sh &
./src/run/simple_compound_target_probability_tox21_xgboost146.sh &
./src/run/simple_compound_target_probability_tox21_xgboost147.sh &
./src/run/simple_compound_target_probability_tox21_xgboost148.sh &
./src/run/simple_compound_target_probability_tox21_xgboost149.sh &
./src/run/simple_compound_target_probability_tox21_xgboost150.sh &
./src/run/simple_compound_target_probability_tox21_xgboost151.sh &
./src/run/simple_compound_target_probability_tox21_xgboost152.sh &
./src/run/simple_compound_target_probability_tox21_xgboost153.sh &
./src/run/simple_compound_target_probability_tox21_xgboost154.sh &
./src/run/simple_compound_target_probability_tox21_xgboost155.sh &
./src/run/simple_compound_target_probability_tox21_xgboost156.sh &
./src/run/simple_compound_target_probability_tox21_xgboost157.sh &
./src/run/simple_compound_target_probability_tox21_xgboost158.sh &
./src/run/simple_compound_target_probability_tox21_xgboost159.sh &
./src/run/simple_compound_target_probability_tox21_xgboost160.sh &
./src/run/simple_compound_target_probability_tox21_xgboost161.sh &
./src/run/simple_compound_target_probability_tox21_xgboost162.sh &
./src/run/simple_compound_target_probability_tox21_xgboost163.sh &
./src/run/simple_compound_target_probability_tox21_xgboost164.sh &
./src/run/simple_compound_target_probability_tox21_xgboost165.sh &
./src/run/simple_compound_target_probability_tox21_xgboost166.sh &
./src/run/simple_compound_target_probability_tox21_xgboost167.sh &
./src/run/simple_compound_target_probability_tox21_xgboost168.sh &
./src/run/simple_compound_target_probability_tox21_xgboost169.sh &
./src/run/simple_compound_target_probability_tox21_xgboost170.sh &
./src/run/simple_compound_target_probability_tox21_xgboost171.sh &
./src/run/simple_compound_target_probability_tox21_xgboost172.sh &
./src/run/simple_compound_target_probability_tox21_xgboost173.sh &
./src/run/simple_compound_target_probability_tox21_xgboost174.sh &
./src/run/simple_compound_target_probability_tox21_xgboost175.sh &
./src/run/simple_compound_target_probability_tox21_xgboost176.sh &
./src/run/simple_compound_target_probability_tox21_xgboost177.sh &
./src/run/simple_compound_target_probability_tox21_xgboost178.sh &
./src/run/simple_compound_target_probability_tox21_xgboost179.sh &
./src/run/simple_compound_target_probability_tox21_xgboost180.sh &
./src/run/simple_compound_target_probability_tox21_xgboost181.sh &
./src/run/simple_compound_target_probability_tox21_xgboost182.sh &
./src/run/simple_compound_target_probability_tox21_xgboost183.sh &
./src/run/simple_compound_target_probability_tox21_xgboost184.sh &
./src/run/simple_compound_target_probability_tox21_xgboost185.sh &
./src/run/simple_compound_target_probability_tox21_xgboost186.sh &
./src/run/simple_compound_target_probability_tox21_xgboost187.sh &
./src/run/simple_compound_target_probability_tox21_xgboost188.sh &
./src/run/simple_compound_target_probability_tox21_xgboost189.sh &
./src/run/simple_compound_target_probability_tox21_xgboost190.sh &
./src/run/simple_compound_target_probability_tox21_xgboost191.sh &
./src/run/simple_compound_target_probability_tox21_xgboost192.sh &
./src/run/simple_compound_target_probability_tox21_xgboost193.sh &
./src/run/simple_compound_target_probability_tox21_xgboost194.sh &
./src/run/simple_compound_target_probability_tox21_xgboost195.sh &
./src/run/simple_compound_target_probability_tox21_xgboost196.sh &
./src/run/simple_compound_target_probability_tox21_xgboost197.sh &
./src/run/simple_compound_target_probability_tox21_xgboost198.sh &
./src/run/simple_compound_target_probability_tox21_xgboost199.sh &
./src/run/simple_compound_target_probability_tox21_xgboost200.sh &
