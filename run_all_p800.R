# This is an overarching script that sources all .R files required to reproduce
# the p800 results and plot. If running in WSL2 using the oneMKL BLAS and LAPACK
# libraries the safest option is to set MKL_DYNAMIC to FALSE and MKL_NUM_THREADS
# to 1 to avoid issues with parallelization.

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

# Step I: data generation ======================================================
source("p800/data_generation/generate_sim_p800_datasets.R"); reset()

# Step II: train/test division =================================================
source("p800/data_generation/traintest_sim_p800_datasets.R"); reset()

# Step III: univariate analysis ================================================
source("p800/analyses/univariate.R"); reset()

# Step IV: benchmark analysis ==================================================
source("p800/analyses/benchmark.R"); reset()

# Step V: gfBLUP analysis ======================================================
source("p800/analyses/gfBLUP_part1.R"); reset()
source("p800/analyses/gfBLUP_part2.R"); reset()
source("p800/misc/gfBLUP_merge.R"); reset()

# Step VI: MegaLMM analysis ====================================================
source("p800/analyses/MegaLMM_part1.R"); reset()
source("p800/analyses/MegaLMM_part2.R"); reset()
source("p800/analyses/MegaLMM_part3.R"); reset()
source("p800/analyses/MegaLMM_part4.R"); reset()
source("p800/misc/MegaLMM_merge.R"); reset()

# Step VII: lsBLUP analysis ====================================================
source("p800/analyses/lsBLUP.R"); reset()

# Step VIII: siBLUP analysis ===================================================
source("p800/analyses/siBLUP_part1.R"); reset()
source("p800/analyses/siBLUP_part2.R"); reset()
source("p800/misc/siBLUP_merge.R"); reset()

# Step IX: MultiMLP analysis ===================================================
source("p800/analyses/CV1_multiMLP_part1.R"); reset()
source("p800/analyses/CV1_multiMLP_part2.R"); reset()
source("p800/misc/CV1_multiMLP_merge.R"); reset()
source("p800/analyses/CV2_multiMLP.R"); reset()

# Step X: phenomics analysis ===================================================
source("p800/analyses/phenomic.R"); reset()

# Step XI: Result Plotting =====================================================
source("p800/plot_p800_results.R"); reset()


