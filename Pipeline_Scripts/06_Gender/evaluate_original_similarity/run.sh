#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Gender/evaluate_original_similarity

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in {1..2}; do # Data Set
  echo "Rscript ${test}/job.R $i" | qsub -N "${i}" -l gres=localhd:1,vmem=50G,mem=50G,walltime=03:00:00:00 -o ${test}/Output -e ${test}/Error
  sleep 0.1
done
