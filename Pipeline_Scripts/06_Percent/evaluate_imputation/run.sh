#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben
project=${home}/Projects/SNF/NM_2015
test=${project}/Scripts/06_Percent/evaluate_imputation

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in 1; do # Data Set
 for j in {1..4}; do # Imputation Method
  for k in {1..3}; do # Member
    for m in {0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5}; do # percent
      for n in {1..25}; do #seed
echo "Rscript ${test}/job.R $i $j $k $m $n" | qsub -N "${i}_${j}_${k}_${m}_${n}" -l nodes=1:ppn=12,gres=localhd:1,vmem=40G,mem=40G,walltime=4:00:00:00 -o ${test}/Output -e ${test}/Error
      sleep 0.1
    done
    sleep 0.5
  done
done
done
done
