#!/bin/bash
#
#SBATCH --array=0-503
#SBATCH --job-name=multi_zoo_model
#SBATCH --output=slurm_%a.out

#SBATCH --mem-per-cpu=40G

#SBATCH --time=0-10:00

module load R
R CMD BATCH slurm_run.R
