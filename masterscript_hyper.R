# This is one of the 7 masterscripts in the repository.
#
# This particular masterscript sources all R scripts to reproduce the hyperspectral
# B5IR and HEAT intermediate results, as well as all the figures related to the
# hyperspectral data:
# - Main text figure 3
# - Main text figure 4
# - Supplementary material figure S1A
# - Supplementary material figure S1B
# - Supplementary material figure S2A
# - Supplementary material figure S2B
# - Supplementary material figure S2C
# - Supplementary material figure S3
# - Supplementary material figure S4
# - Supplementary material figure S5
# - Supplementary material figure S6
# - Supplementary material figure S7A
# - Supplementary material figure S7B
# - Supplementary material figure S7C
# - Supplementary material figure S9
# - Supplementary material figure S10
#
# There are 5 sections that should ideally be run in that order if you are interested
# in reproducing all intermediate results. For reproducibility spot checks, make sure
# that the repository either includes all of the datafiles in the folders
# `hyper_1415B5IR/datasets/splines/`, `hyper_1415B5IR/datasets/VEGsplines/`,
# `hyper_1415HEAT/datasets/nosplines/`, and `hyper_1415HEAT/datasets/splines/`, or
# run the following scripts below:
# - Hyper B5IR, Step II: data generation for B5IR
# - Hyper HEAT, Step II: data generation for HEAT
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
if (!("gfBLUPold" %in% installed.packages())) {
  install.packages("gfBLUPold_1.3.1.tar.gz",
                   type = "source", repos = NULL, dependencies = TRUE)
}

# Checking if the glfBLUP package is installed:
if (!("glfBLUP" %in% installed.packages())) {
  install.packages("glfBLUP_1.0.0.tar.gz",
                   type = "source", repos = NULL, dependencies = TRUE)
}



# Hyper B5IR ===================================================================
## Step I: data pre-processing for B5IR ========================================
# This script preprocesses the raw hyperspectral B5IR data. The output of this
# script should already be present in the repository, so running this script is
# not strictly required (and not possible for the public repository due to the lack
# of raw data).
if (!(length(list.files("hyper_datafiles/")) == 1)) {
  source("hyper_1415B5IR/data_generation/hyper_preprocessing.R"); reset()
}

## Step II: data generation for B5IR ===========================================
# This script uses the output of the previous script to create all 250 datafiles
# for the hyperspectral B5IR dataset.
source("hyper_1415B5IR/data_generation/generate_hyper_datasets.R"); reset()

## Step III: univariate analysis ===============================================
# This script runs the 250 replications for the univariate analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/univariate_hyper.R"); reset()

## Step IV: lsBLUP analysis ====================================================
# This script runs the 250 replications for the CV1 and CV2 lsBLUP analyses on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/lsBLUP_hyper.R"); reset()
# This script runs the 250 replications for the CV2VEG lsBLUP analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/lsBLUP_hyper_CV2VEG.R"); reset()

## Step V: siBLUP analysis =====================================================
# This script runs the 250 replications for the CV1 and CV2 siBLUP analyses on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/siBLUP_hyper.R"); reset()
# This script runs the 250 replications for the CV2VEG siBLUP analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/siBLUP_hyper_CV2VEG.R"); reset()

## Step VI: glfBLUP analysis ===================================================
# This script runs the 250 replications for the CV1 and CV2 glfBLUP analyses on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/glfBLUP_hyper.R"); reset()
# This script runs the 250 replications for the CV2VEG glfBLUP analyses on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/glfBLUP_hyper_CV2VEG.R"); reset()

## Step VII: MegaLMM CV1 analysis ==============================================
# This script runs the 250 replications for the CV1 MegaLMM analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/MegaLMM_hyper_CV1.R"); reset()

## Step VIII: MegaLMM CV2 analysis =============================================
# This script runs the 250 replications for the CV2 MegaLMM analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/MegaLMM_hyper_CV2.R"); reset()
# This script runs the 250 replications for the CV2VEG MegaLMM analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/MegaLMM_hyper_CV2VEG.R"); reset()

## Step IX: MultiMLP CV1 analysis ==============================================
# This script runs the 250 replications for the CV1 multiMLP analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/multiMLP_hyper_CV1.R"); reset()

## Step X: MultiMLP CV2 analysis ===============================================
# This script runs the 250 replications for the CV2 multiMLP analysis on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/multiMLP_hyper_CV2.R"); reset()
# This script runs the 250 replications for the CV2VEG multiMLP analyses on the
# hyperspectral B5IR dataset.
source("hyper_1415B5IR/analyses/multiMLP_hyper_CV2_VEG.R"); reset()



# Hyper HEAT ===================================================================
## Step I: data pre-processing for HEAT ========================================
# This script preprocesses the raw hyperspectral HEAT data. The output of this
# script should already be present in the repository, so running this script is
# not strictly required (and not possible for the public repository due to the lack
# of raw data).
if (!(length(list.files("hyper_datafiles/")) == 1)) {
  source("hyper_1415HEAT/data_generation/hyper_preprocessing.R"); reset()
}

## Step II: data generation for HEAT ===========================================
# This script uses the output of the previous script to create all 250 datafiles
# for the hyperspectral HEAT dataset.
source("hyper_1415HEAT/data_generation/generate_hyper_datasets.R"); reset()

## Step III: univariate analysis ===============================================
# This script runs the 250 replications for the univariate analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/univariate_hyper.R"); reset()

## Step IV: lsBLUP analysis ====================================================
# This script runs the 250 replications for the CV1 and CV2 lsBLUP analyses on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/lsBLUP_hyper.R"); reset()
# This script runs the 250 replications for the CV2VEG lsBLUP analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/lsBLUP_hyper_CV2VEG.R"); reset()

## Step V: siBLUP analysis =====================================================
# This script runs the 250 replications for the CV1 and CV2 siBLUP analyses on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/siBLUP_hyper.R"); reset()
# This script runs the 250 replications for the CV2VEG siBLUP analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/siBLUP_hyper_CV2VEG.R"); reset()

## Step VI: glfBLUP analysis ===================================================
# This script runs the 250 replications for the CV1 and CV2 glfBLUP analyses on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/glfBLUP_hyper.R"); reset()
# This script runs the 250 replications for the CV2VEG glfBLUP analyses on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/glfBLUP_hyper_CV2VEG.R"); reset()

## Step VII: MegaLMM CV1 analysis ==============================================
# This script runs the 250 replications for the CV1 MegaLMM analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/MegaLMM_hyper_CV1.R"); reset()

## Step VIII: MegaLMM CV2 analysis =============================================
# This script runs the 250 replications for the CV2 MegaLMM analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/MegaLMM_hyper_CV2.R"); reset()
# This script runs the 250 replications for the CV2VEG MegaLMM analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/MegaLMM_hyper_CV2VEG.R"); reset()

## Step IX: MultiMLP CV1 analysis ==============================================
# This script runs the 250 replications for the CV1 multiMLP analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/multiMLP_hyper_CV1.R"); reset()

## Step X: MultiMLP CV2 analysis ===============================================
# This script runs the 250 replications for the CV2 multiMLP analysis on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/multiMLP_hyper_CV2.R"); reset()
# This script runs the 250 replications for the CV2VEG multiMLP analyses on the
# hyperspectral HEAT dataset.
source("hyper_1415HEAT/analyses/multiMLP_hyper_CV2_VEG.R"); reset()



# Hyper B5IR and HEAT plotting =================================================
## Step I: Result Plotting =====================================================
# This script produces figure 3 of the main text.
source("plot_hyper.R"); reset()



# Single date analyses and plotting ============================================
## Step I: lsBLUP single date analysis and plotting ============================
# This script runs lsBLUP on the hyperspectral B5IR data of a single date, and then
# produces figure S1A of the supplementary material.
source("hyper_B5IR/analyses_single_date/lsBLUP_hyper_single_date.R"); reset()

## Step II: siBLUP single date analysis and plotting ===========================
# This script runs siBLUP on the hyperspectral B5IR data of a single date, and then
# produces figure S1B of the supplementary material.
source("hyper_B5IR/analyses_single_date/siBLUP_hyper_single_date.R"); reset()

## Step III: glfBLUP single date analysis and plotting =========================
# This script runs glfBLUP on the hyperspectral B5IR data of a single date, and then
# produces figure 4 of the main text.
source("hyper_B5IR/analyses_single_date/glfBLUP_hyper_single_date.R"); reset()

## Step IV: MegaLMM single date analyses =======================================
# This script runs MegaLMM on the hyperspectral B5IR data of a single date.
source("hyper_B5IR/analyses_single_date/MegaLMM_hyper_single_date.R"); reset()

## Step V: MegaLMM single date plotting ========================================
# This script uses the output of the previous script to produce the following figures
# of the supplementary material:
# - Figure S2B
# - Figure S2C
# - Figure S2A
# - Figure S7B
# - Figure S4
# - Figure S7C
# - Figure S5
# - Figure S6
# - Figure S7A
# - Figure S3
source("hyper_B5IR/analyses_single_date/MegaLMM_hyper_single_date_plotting.R"); reset()


