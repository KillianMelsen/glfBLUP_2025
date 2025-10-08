# This is one of the 7 masterscripts in the repository.
#
# This particular masterscript sources all R scripts to reproduce the intermediate
# results for the simulated data with random residual structure (p800),
# as well as figure 2 and 5 of the main text.
#
# The scripts sourced below should ideally be run in that order if you are interested
# in reproducing all intermediate results. For reproducibility spot checks, make sure
# that the repository either includes all of the datafiles in the folder
# `p800/datasets/`, or run the following scripts below:
# - Data generation, analyses, and plotting, Step I: data generation
# - Data generation, analyses, and plotting, Step II: train/test division
#
# For more information, please refer to the readme, the short comments before each
# step in this file, or the larger comments at the top of each of the scripts sourced
# below.

# Preliminaries ================================================================
# Getting the project working dir:
wd <- getwd()

# Function to reset the environment:
reset <- function() {
  rm(list = setdiff(ls(".GlobalEnv"), c("wd", "reset")), pos = ".GlobalEnv")
  tictoc::tic.clearlog()
  setwd(wd)
}

# Checking if the old package with some helper functions is installed:
if (!("gfBLUPold" %in% installed.packages()[,1])) {
  install.packages("gfBLUPold_1.3.1.tar.gz", type = "source", repos = NULL)
}



# Data generation, analyses, and plotting ======================================
## Step I: data generation =====================================================
# This script simulates all datasets using a random residual structure.
source("p800/data_generation/generate_sim_p800_datasets.R"); reset()

## Step II: train/test division ================================================
# This script samples the training and test sets for all replications.
source("p800/data_generation/traintest_sim_p800_datasets.R"); reset()

## Step III: univariate analysis ===============================================
# This script runs all replications for the univariate model on all the simulated
# p800 data.
source("p800/analyses/univariate.R"); reset()

## Step IV: benchmark analysis =================================================
# This script runs all replications for the benchmark model on all the simulated
# p800 data.
source("p800/analyses/benchmark.R"); reset()

## Step V: glfBLUP analysis ====================================================
# This script runs all replications for the glfBLUP model on all the simulated
# p800 data.
source("p800/analyses/glfBLUP.R"); reset()

## Step VI: MegaLMM analysis ===================================================
# The 5 scripts below run all replications for the MegaLMM model on all the simulated
# p800 data.
source("p800/analyses/MegaLMM_part1.R"); reset()
source("p800/analyses/MegaLMM_part2.R"); reset()
source("p800/analyses/MegaLMM_part3.R"); reset()
source("p800/analyses/MegaLMM_part4.R"); reset()
source("p800/analyses/MegaLMM_part5.R"); reset()
# This script simply merges the intermediate results produced by the 5 scripts above.
source("p800/misc/MegaLMM_merge.R"); reset()

## Step VII: lsBLUP analysis ===================================================
# This script runs all replications for the lsBLUP model on all the simulated
# p800 data.
source("p800/analyses/lsBLUP.R"); reset()

## Step VIII: siBLUP analysis ==================================================
# The 2 scripts below run all replications for the siBLUP model on all the simulated
# p800 data.
source("p800/analyses/siBLUP_part1.R"); reset()
source("p800/analyses/siBLUP_part2.R"); reset()
# This script simply merges the intermediate results produced by the 2 scripts above.
source("p800/misc/siBLUP_merge.R"); reset()

## Step IX: MultiMLP analysis ==================================================
# The 2 scripts below run all CV1 replications for the multiMLP model on all the
# simulated p800 data.
source("p800/analyses/CV1_multiMLP_part1.R"); reset()
source("p800/analyses/CV1_multiMLP_part2.R"); reset()
# The 2 scripts below run all CV2 replications for the multiMLP model on all the
# simulated p800 data.
source("p800/analyses/CV2_multiMLP_part1.R"); reset()
source("p800/analyses/CV2_multiMLP_part2.R"); reset()
# The 2 scripts below simply merge the intermediate results produced by the scripts
# above.
source("p800/misc/CV1_multiMLP_merge.R"); reset()
source("p800/misc/CV2_multiMLP_merge.R"); reset()

## Step X: Result Plotting =====================================================
# This script uses all p800 intermediate results to produce figure 2 and 5 of
# the main text.
source("p800/plot_p800_results.R"); reset()


