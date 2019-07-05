#!/bin/bash

home=/hpf/largeprojects/agoldenb/daniel
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Two_Thousand_Features/evaluate_imputation
results=${test}/Results

# Clear output from previous runs
rm ${test}/Results/*.txt

for i in {1..4}; do # Data Set
  for j in {1..4}; do # Imputation Method
    for k in {0..49}; do # Seed
      cat "${results}/Clustering/${i}_${j}_${k}.txt" >> "${results}/clustering.txt"
      cat "${results}/Imputation/${i}_${j}_${k}.txt" >> "${results}/imputation.txt"
    done
  done
done
