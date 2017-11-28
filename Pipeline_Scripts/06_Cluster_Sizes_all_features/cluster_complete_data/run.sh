#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Cluster_Sizes/cluster_complete_data

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in {1..5}; do # Data Set
 for j in {1..4}; do # Index for cluster sizes
  echo "${test}/job.R $i $j" | qsub -N "${i}_${j}" -l nodes=1:ppn=35,gres=localhd:1,vmem=70G,mem=70G,walltime=4:00:00:00 -o ${test}/Output -e ${test}/Error
  sleep 0.1
done
done