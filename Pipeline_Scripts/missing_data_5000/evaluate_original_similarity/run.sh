#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/missing_data/evaluate_original_similarity

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in 3; do # Data Set
  echo "${test}/job.R $i" | qsub -N "${i}" -l gres=localhd:1,vmem=8G,mem=8G,walltime=04:00:00 -o ${test}/Output -e ${test}/Error
  sleep 0.1
done
