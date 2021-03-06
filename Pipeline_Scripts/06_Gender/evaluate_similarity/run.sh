#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Gender/evaluate_similarity

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in {1..2}; do # Data Set
  for j in {0..49}; do # Seed
    echo "Rscript ${test}/job.R $i $j" | qsub -N "${i}_${j}" -l gres=localhd:1,vmem=50G,mem=50G,walltime=03:00:00 -o ${test}/Output -e ${test}/Error
    sleep 0.1
  done
  sleep 0.5
done
