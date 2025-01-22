#!/bin/bash
#SBATCH --job-name="01_classification_scenario1"
#SBATCH -n 60
#SBATCH -t 4-00:00:00 
#SBATCH -o 01_classification_scenario1.out
cd /public/home/zoujian/2501_PRECISION.ML

module load apps/R/4.3.1

Rscript --vanilla 01_classification_scenario1.R