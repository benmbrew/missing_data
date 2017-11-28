#!/bin/bash

home=/home/benbrew/Documents
project=${home}/NM_2015_ben
test=${project}/Scripts/06_two_thousand_features/evaluate_original_imputation
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt

for i in {1..4}; do # Data Set
  for j in {1..4}; do # Imputation Method
    cat "${results}/Clustering/${i}_${j}.txt" >> "${results}/clustering.txt"
  done
done
