#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/missing_data_combat/cluster_complete_data

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in {1..4}; do #clustersize 
 echo "Rscript ${test}/job.R $i" | qsub -N "${i}" -l nodes=1:ppn=35,gres=localhd:1,vmem=70G,mem=70G,walltime=1:00:00:00 -o ${test}/Output -e ${test}/Error
done
