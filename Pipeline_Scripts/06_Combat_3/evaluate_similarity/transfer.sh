#!/bin/bash

home=/home/benbrew/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Combat_3/evaluate_similarity
results=${test}/Results

# Clear output from previous runs
#rm ${test}/Results/*.txt

for i in {0..49}; do # Seed
   cat "${results}/Similarity/${i}.txt" >> "${results}/similarity.txt"
done




