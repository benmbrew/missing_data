#!/bin/bash

home=/home/benbrew/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/missing_data_combat/evaluate_imputation
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt


for i in {1..4}; do # Imputation Method
  for j in {0..49}; do # Seed
     cat "${results}/Clustering/${i}_${j}.txt" >> "${results}/clustering.txt"
     cat "${results}/Imputation/${i}_${j}.txt" >> "${results}/imputation.txt"
	done
done
