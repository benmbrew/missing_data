#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Combat_Data_Types/cluster_complete_data

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs

echo "Rscript ${test}/job.R" | qsub -N "job.R" -l nodes=1:ppn=35,gres=localhd:1,vmem=70G,mem=70G,walltime=4:00:00:00 -o ${test}/Output -e ${test}/Error
