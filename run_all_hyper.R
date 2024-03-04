# This is an overarching script that sources all .R files required to reproduce
# the hyper results and plot. If running in WSL2 using the oneMKL BLAS and LAPACK
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
if (!("gfBLUpold" %in% installed.packages())) {
  install.packages("gfBLUPold_1.3.1.tar.gz", type = "source", repos = NULL)
}

# Step I: data generation ======================================================
source("hyper/data_generation/generate_hyper_datasets.R"); reset()