#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=p800_lowrank_MegaLMM

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_p800_lowrank_MegaLMM_%A.%a.txt
#SBATCH --error=output/errors_p800_lowrank_MegaLMM_%A.%a.txt

#-----------------------------Other information------------------------
#SBATCH --comment='p800_lowrank_MegaLMM'

#-----------------------------Required resources-----------------------
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-5

#-----------------------------Environment, Operations and Job steps----
# Modules
module reset
module load 2024
module load R/4.4.2
module load R-bundle-CRAN/2024
module load CMake
module load GSL

# Converting to unix and running R script
dos2unix ./MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part${SLURM_ARRAY_TASK_ID}.R
chmod +x ./MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part${SLURM_ARRAY_TASK_ID}.R
./MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part${SLURM_ARRAY_TASK_ID}.R


