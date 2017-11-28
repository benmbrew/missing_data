#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Impute/cluster_complete_data

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in {1..4}; do # Data Set
  echo "${test}/job.R $i" | qsub -N "${i}" -l nodes=1:ppn=35,gres=localhd:1,vmem=8G,mem=8G,walltime=2:00:00 -o ${test}/Output -e ${test}/Error
  sleep 0.1
done
