# This is one of the 7 masterscripts in the repository.
#
# This particular masterscript sources all R scripts to reproduce the timing
# results using simulated data with random residual structure and different,
# numbers of secondary features.
#
# Note that timing results are obviously not reproducible. Timing results also
# depend on hardware.
#
# This masterscript also produces figure S8 of the supplementary material.

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
if (!("gfBLUPold" %in% installed.packages())) {
  install.packages("gfBLUPold_1.3.1.tar.gz", type = "source", repos = NULL)
}



# Data generation, timing, and plotting ========================================
## Step I: data generation =====================================================
# This script generates simulated datasets with random residual structure with
# different numbers of secondary features.
source("timing/generate_timing_data.R"); reset()

## Step II: train/test division of the data ====================================
# This script randomly samples a training and test set for each of the generated
# datasets.
source("timing/traintest_timing_data.R"); reset()

## Step III: timing the analyses ===============================================
# This script performs the timings of glfBLUP using the generated datasets.
source("timing/timing.R"); reset()

## Step IV: plotting the results ===============================================
# This script uses the timing results to produce figure S8 of the supplementary
# material.
source("timing/plotting.R"); reset()


