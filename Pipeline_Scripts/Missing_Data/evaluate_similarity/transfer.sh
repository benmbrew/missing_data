#!/bin/bash

home=/home/benbrew/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/Missing_Data/evaluate_similarity
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt

for i in {1..5}; do # Data Set
  for j in {0..49}; do # Seed
    cat "${results}/Similarity/${i}_${j}.txt" >> "${results}/similarity.txt"
  done
done




