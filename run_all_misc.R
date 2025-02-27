# This is an overarching script that sources all .R files required to reproduce
# some miscellaneous results and plots.

# Getting the project working dir:
wd <- getwd()

# Function to reset the environment:
reset <- function() {
  rm(list = setdiff(ls(".GlobalEnv"), c("wd", "reset")), pos = ".GlobalEnv")
  tictoc::tic.clearlog()
  setwd(wd)
}



# Step I: plot feature filtering for B5IR ======================================
source("hyper_B5IR/plot_RF.R"); reset()

# Step II: plot feature filtering for HEAT =====================================
source("hyper_HEAT/plot_RF.R"); reset()


