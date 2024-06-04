# This is an overarching script that sources all .R files required to reproduce
# the timing results and plot. Run in Windows 11 Pro 23H2 (i5-13600KF).

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



# Step I: data generation ======================================================
source("timing/generate_timing_data.R"); reset()

# Step II: train/test division of the data =====================================
source("timing/traintest_timing_data.R"); reset()

# Step III: timing the analyses ================================================
source("timing/timing.R"); reset()

# Step IV: plotting the results ================================================
source("timing/plotting.R"); reset()


