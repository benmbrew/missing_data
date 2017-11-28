#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Percent/evaluate_similarity

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in 1 2; do # Data Set
 for j in {1..3}; do #member
  for k in {1..10}; do #percent
   for m in {0..3}; do # Seed
    echo "${test}/job.R $i $j $k $m" | qsub -N "${i}_${j}_${k}_${m}" -l gres=localhd:1,vmem=8G,mem=8G,walltime=04:00:00 -o ${test}/Output -e ${test}/Error
    sleep 0.1
  done
  sleep 0.5
done
done 
done
