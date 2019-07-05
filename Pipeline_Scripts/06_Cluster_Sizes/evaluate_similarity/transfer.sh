#!/bin/bash

home=/hpf/largeprojects/agoldenb/daniel
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Two_Thousand_Features/evaluate_similarity
results=${test}/Results

# Clear output from previous runs
rm ${test}/Results/*.txt

for i in {1..4}; do # Data Set
  for j in {0..49}; do # Seed
    cat "${results}/Similarity/${i}_${j}.txt" >> "${results}/similarity.txt"
  done
done




