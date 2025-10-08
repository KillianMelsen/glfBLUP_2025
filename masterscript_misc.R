# This is one of the 7 masterscripts in the repository.
#
# This particular masterscript sources two R scripts. The first one produces
# figure S9 of the supplementary material while the second script produces
# figure S10 of the supplementary material.
#
# For more information, please refer to the readme, the short comments before each
# step in this file, or the larger comments at the top of each of the scripts sourced
# below.

# Getting the project working dir:
wd <- getwd()

# Function to reset the environment:
reset <- function() {
  rm(list = setdiff(ls(".GlobalEnv"), c("wd", "reset")), pos = ".GlobalEnv")
  tictoc::tic.clearlog()
  setwd(wd)
}



# Step I: plot feature filtering for B5IR ======================================
# This script produces figure S9 of the supplementary material
source("hyper_B5IR/plot_RF.R"); reset()

# Step II: plot feature filtering for HEAT =====================================
# This script produces figure S10 of the supplementary material
source("hyper_HEAT/plot_RF.R"); reset()


