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

# Step II: univariate analysis =================================================
source("hyper/analyses/univariate_hyper.R"); reset()

# Step III: lsBLUP analysis ====================================================
source("hyper/analyses/lsBLUP_hyper.R"); reset()

# Step IV: siBLUP analysis =====================================================
source("hyper/analyses/siBLUP_hyper.R"); reset()

# Step V: gfBLUP analysis ======================================================
source("hyper/analyses/gfBLUP_hyper.R"); reset()

# Step VI: MegaLMM CV1 analysis ================================================
source("hyper/analyses/MegaLMM_hyper_CV1_RF.R"); reset()

# Step VII: MegaLMM CV2 analysis ===============================================
source("hyper/analyses/MegaLMM_hyper_CV2_RF.R"); reset()

# Step VIII: MultiMLP CV1 analysis =============================================
source("hyper/analyses/multiMLP_hyper_CV1_RF.R"); reset()

# Step IX: MultiMLP CV2 analysis ===============================================
source("hyper/analyses/multiMLP_hyper_CV2_RF.R"); reset()

# Step VIII: Result Plotting ===================================================
source("hyper/plot_hyper_results.R"); reset()


