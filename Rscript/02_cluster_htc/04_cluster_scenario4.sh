#!/bin/bash
#SBATCH --job-name="04_cluster_scenario4"
#SBATCH -n 60
#SBATCH -t 4-00:00:00 
#SBATCH -o 04_cluster_scenario4.out
cd /public/home/zoujian/2501_PRECISION.ML

module load apps/R/4.3.1

Rscript --vanilla 04_cluster_scenario4.R