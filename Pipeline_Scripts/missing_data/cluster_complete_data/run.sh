#!/bin/bash

home= <ENTER PATH TO cluster_complete_data HERE>

# Run the jobs
for i in 1,2; do # Data Set
  echo "${home}/job.R $i" | qsub -N "${i}" -l nodes=1:ppn=35,gres=localhd:1,
  vmem=70G,mem=70G,walltime=4:00:00:00 -o ${home}/Output -e ${home}/Error
  sleep 0.1
done
