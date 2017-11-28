#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Two_Thousand_Features_Dup/evaluate_imputation

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in 4; do # Data Set
  for j in {1..4}; do # Imputation Method
    for k in {0..49}; do # Seed
      echo "${test}/job.R $i $j $k" | qsub -N "${i}_${j}_${k}" -l nodes=1:ppn=12,gres=localhd:1,vmem=90G,mem=90G,walltime=7:00:00:00 -o ${test}/Output -e ${test}/Error
      sleep 0.1
    done
    sleep 0.5
  done
done
