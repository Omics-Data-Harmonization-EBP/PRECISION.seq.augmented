#!/bin/bash
#SBATCH --job-name="05_cluster_scenario5"
#SBATCH -n 60
#SBATCH -t 4-00:00:00 
#SBATCH -o 05_cluster_scenario5.out
cd /public/home/zoujian/2501_PRECISION.ML

module load apps/R/4.3.1

Rscript --vanilla 05_cluster_scenario5.R