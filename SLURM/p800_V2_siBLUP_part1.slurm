#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=glfBLUP_p800_V2_siBLUP_part1

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_p800_V2_siBLUP_part1_%j.txt
#SBATCH --error=output/errors_p800_V2_siBLUP_part1_%j.txt

#-----------------------------Other information------------------------
#SBATCH --comment='p800_V2 siBLUP part 1 script'

#-----------------------------Required resources-----------------------
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=80G

#-----------------------------Environment, Operations and Job steps----
# Modules
module reset
module load 2024
module load R/4.4.2
module load R-bundle-CRAN/2024
module load CMake
module load GSL

# Converting to unix and running R script
dos2unix ./p800_V2/analyses/siBLUP_part1.R
chmod +x ./p800_V2/analyses/siBLUP_part1.R
./p800_V2/analyses/siBLUP_part1.R


