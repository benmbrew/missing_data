#!/bin/bash

home=/home/benbrew/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Two_Thousand_Features_Dup_3000/evaluate_original_imputation
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt

for i in {1..5}; do # Data Set
  for j in {1..4}; do # Imputation Method
    cat "${results}/Clustering/${i}_${j}.txt" >> "${results}/clustering.txt"
  done
done
