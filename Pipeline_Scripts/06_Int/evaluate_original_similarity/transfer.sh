#!/bin/bash

home=/home/benbrew/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Int/evaluate_original_similarity
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt

for i in {1..5}; do # Data Set
  cat "${results}/Similarity/${i}.txt" >> "${results}/similarity.txt"
done




