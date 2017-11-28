#!/bin/bash

home=/home/benbrew/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Combat_Seed/evaluate_original_imputation
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt

for i in {1..4}; do # Imputation Method
  for j in {1..49}; do # seed 
   cat "${results}/Clustering/${i}_${j}.txt" >> "${results}/clustering.txt"
 done
done

