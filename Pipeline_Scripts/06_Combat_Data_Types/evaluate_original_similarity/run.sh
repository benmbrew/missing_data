#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Combat_Data_Types/evaluate_original_similarity

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
echo "Rscript ${test}/job.R" | qsub -N "com_orig_sim" -l gres=localhd:1,vmem=8G,mem=8G,walltime=04:00:00 -o ${test}/Output -e ${test}/Error
