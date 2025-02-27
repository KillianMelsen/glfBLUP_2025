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
if (!("gfBLUPold" %in% installed.packages())) {
  install.packages("gfBLUPold_1.3.1.tar.gz",
                   type = "source", repos = NULL, dependencies = TRUE)
}

# Checking if the glfBLUP package is installed:
if (!("glfBLUP" %in% installed.packages())) {
  install.packages("glfBLUP_1.0.0.tar.gz",
                   type = "source", repos = NULL, dependencies = TRUE)
}



# Hyper b5IR ===================================================================
## Step I: data pre-processing for B5IR ========================================
source("hyper_1415B5IR/data_generation/hyper_preprocessing.R"); reset()

## Step II: data generation for B5IR ===========================================
source("hyper_1415B5IR/data_generation/generate_hyper_datasets.R"); reset()

## Step III: univariate analysis ===============================================
source("hyper_1415B5IR/analyses/univariate_hyper.R"); reset()

## Step IV: lsBLUP analysis ====================================================
source("hyper_1415B5IR/analyses/lsBLUP_hyper.R"); reset()
source("hyper_1415B5IR/analyses/lsBLUP_hyper_CV2VEG.R"); reset()

## Step V: siBLUP analysis =====================================================
source("hyper_1415B5IR/analyses/siBLUP_hyper.R"); reset()
source("hyper_1415B5IR/analyses/siBLUP_hyper_CV2VEG.R"); reset()

## Step VI: glfBLUP analysis ===================================================
source("hyper_1415B5IR/analyses/glfBLUP_hyper.R"); reset()
source("hyper_1415B5IR/analyses/glfBLUP_hyper_CV2VEG.R"); reset()

## Step VII: MegaLMM CV1 analysis ==============================================
source("hyper_1415B5IR/analyses/MegaLMM_hyper_CV1.R"); reset()

## Step VIII: MegaLMM CV2 analysis =============================================
source("hyper_1415B5IR/analyses/MegaLMM_hyper_CV2.R"); reset()
source("hyper_1415B5IR/analyses/MegaLMM_hyper_CV2VEG.R"); reset()

## Step IX: MultiMLP CV1 analysis ==============================================
source("hyper_1415B5IR/analyses/multiMLP_hyper_CV1.R"); reset()

## Step X: MultiMLP CV2 analysis ===============================================
source("hyper_1415B5IR/analyses/multiMLP_hyper_CV2.R"); reset()
source("hyper_1415B5IR/analyses/multiMLP_hyper_CV2_VEG.R"); reset()



# Hyper HEAT ===================================================================
## Step I: data pre-processing for HEAT ========================================
source("hyper_1415HEAT/data_generation/hyper_preprocessing.R"); reset()

## Step II: data generation for HEAT ===========================================
source("hyper_1415HEAT/data_generation/generate_hyper_datasets.R"); reset()

## Step III: univariate analysis ===============================================
source("hyper_1415HEAT/analyses/univariate_hyper.R"); reset()

## Step IV: lsBLUP analysis ====================================================
source("hyper_1415HEAT/analyses/lsBLUP_hyper.R"); reset()
source("hyper_1415HEAT/analyses/lsBLUP_hyper_CV2VEG.R"); reset()

## Step V: siBLUP analysis =====================================================
source("hyper_1415HEAT/analyses/siBLUP_hyper.R"); reset()
source("hyper_1415HEAT/analyses/siBLUP_hyper_CV2VEG.R"); reset()

## Step VI: glfBLUP analysis ===================================================
source("hyper_1415HEAT/analyses/glfBLUP_hyper.R"); reset()
source("hyper_1415HEAT/analyses/glfBLUP_hyper_CV2VEG.R"); reset()

## Step VII: MegaLMM CV1 analysis ==============================================
source("hyper_1415HEAT/analyses/MegaLMM_hyper_CV1.R"); reset()

## Step VIII: MegaLMM CV2 analysis =============================================
source("hyper_1415HEAT/analyses/MegaLMM_hyper_CV2.R"); reset()
source("hyper_1415HEAT/analyses/MegaLMM_hyper_CV2VEG.R"); reset()

## Step IX: MultiMLP CV1 analysis ==============================================
source("hyper_1415HEAT/analyses/multiMLP_hyper_CV1.R"); reset()

## Step X: MultiMLP CV2 analysis ===============================================
source("hyper_1415HEAT/analyses/multiMLP_hyper_CV2.R"); reset()
source("hyper_1415HEAT/analyses/multiMLP_hyper_CV2_VEG.R"); reset()



# Hyper B5IR and HEAT plotting =================================================
## Step I: Result Plotting =====================================================
source("plot_hyper.R"); reset()



# Single date analyses and plotting ============================================
## Step I: lsBLUP single date analysis and plotting ============================
source("hyper_B5IR/analyses_single_date/lsBLUP_hyper_single_date.R"); reset()

## Step II: siBLUP single date analysis and plotting ===========================
source("hyper_B5IR/analyses_single_date/siBLUP_hyper_single_date.R"); reset()

## Step III: glfBLUP single date analysis and plotting =========================
source("hyper_B5IR/analyses_single_date/glfBLUP_hyper_single_date.R"); reset()

## Step IV: MegaLMM single date analysis and plotting ==========================
source("hyper_B5IR/analyses_single_date/MegaLMM_hyper_single_date.R"); reset()


