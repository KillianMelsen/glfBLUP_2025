#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=glfBLUP_p800_V2_MegaLMM

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_p800_V2_MegaLMM_%A.%a.txt
#SBATCH --error=output/errors_p800_V2_MegaLMM_%A.%a.txt

#-----------------------------Other information------------------------
#SBATCH --comment='p800_V2 MegaLMM scripts'

#-----------------------------Required resources-----------------------
#SBATCH --time=4-00:00:00
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
dos2unix ./p800_V2/analyses/MegaLMM_part${SLURM_ARRAY_TASK_ID}.R
chmod +x ./p800_V2/analyses/MegaLMM_part${SLURM_ARRAY_TASK_ID}.R
./p800_V2/analyses/MegaLMM_part${SLURM_ARRAY_TASK_ID}.R


