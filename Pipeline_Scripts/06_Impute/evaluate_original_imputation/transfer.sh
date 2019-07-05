#!/bin/bash

home=/hpf/largeprojects/agoldenb/daniel
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Two_Thousand_Features/evaluate_original_imputation
results=${test}/Results

# Clear output from previous runs
rm ${test}/Results/*.txt

for i in {1..4}; do # Data Set
  for j in {1..4}; do # Imputation Method
    cat "${results}/Clustering/${i}_${j}.txt" >> "${results}/clustering.txt"
  done
done
