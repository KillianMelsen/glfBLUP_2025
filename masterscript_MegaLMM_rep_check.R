# This is one of the 7 masterscripts in the repository.
#
# This particular masterscript sources all R scripts required to run the
# MegaLMM hyper B5IR/HEAT and p800/p800_lowrank reproducibility checks.
#
# This is necessary because MegaLMM uses MCMC sampling which is unfortunately not
# fixed by setting the seed in R due to monte carlo errors. Note that this is a
# limitation of the MegaLMM software and outside of our control.
#
# Running a second set of independent analyses on all datasets allows us to compare
# the mean prediction accuracies reported in the manuscript versus the mean accuracies
# obtained from the second set of independent analyses. The differences between them
# should be very small.
#
# IMPORTANT: note that all datafiles for the hyperspectral B5IR/HEAT and p800/p800_lowrank
# datasets must be available. If these are not present in the repository, please first
# generate them using the appropriate scripts (see the readme).
#
# The preliminaries section checks whether some packages are installed and defines
# a function used to reset the environment between sourcing scripts.
#
# The analyses section then sources the scripts producing all MegaLMM analyses
# for the different datasets and merges some of the intermediate result files.
#
# Finally, the reproducibility checks section shows that the differences between
# mean accuracies obtained from the two independent sets of analyses (i.e., those
# reported in the manuscript and those done from this masterscript) are negligible.
#
# Note that this repository already includes all results from both sets of analyses,
# so the differences between means can immediately be checked using the reproducibility
# checks section.
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



# Analyses =====================================================================
## Hyper B5IR ==================================================================
### Step I: CV1 ================================================================
# This script performs the second, independent set of CV1 MegaLMM analyses for the
# hyperspectral B5IR data.
source("MegaLMM_rep_check_scripts/hyper_1415B5IR/analyses/MegaLMM_hyper_CV1.R"); reset()

### Step II: CV2 ===============================================================
# This script performs the second, independent set of CV2 MegaLMM analyses for the
# hyperspectral B5IR data.
source("MegaLMM_rep_check_scripts/hyper_1415B5IR/analyses/MegaLMM_hyper_CV2.R"); reset()

### Step III: CV2VEG ===========================================================
# This script performs the second, independent set of CV2VEG MegaLMM analyses for the
# hyperspectral B5IR data.
source("MegaLMM_rep_check_scripts/hyper_1415B5IR/analyses/MegaLMM_hyper_CV2VEG.R"); reset()



## Hyper HEAT ==================================================================
### Step I: CV1 ================================================================
# This script performs the second, independent set of CV1 MegaLMM analyses for the
# hyperspectral HEAT data.
source("MegaLMM_rep_check_scripts/hyper_1415HEAT/analyses/MegaLMM_hyper_CV1.R"); reset()

### Step II: CV2 ===============================================================
# This script performs the second, independent set of CV2 MegaLMM analyses for the
# hyperspectral HEAT data.
source("MegaLMM_rep_check_scripts/hyper_1415HEAT/analyses/MegaLMM_hyper_CV2.R"); reset()

### Step III: CV2VEG ===========================================================
# This script performs the second, independent set of CV2VEG MegaLMM analyses for the
# hyperspectral HEAT data.
source("MegaLMM_rep_check_scripts/hyper_1415HEAT/analyses/MegaLMM_hyper_CV2VEG.R"); reset()



## p800 ========================================================================
### Step I: Run part 1 through 5 ===============================================
# The 5 scripts below perform the second, independent set of MegaLMM analyses for the
# simulated data with random residual structure (p800).
source("MegaLMM_rep_check_scripts/p800/analyses/MegaLMM_part1.R"); reset()
source("MegaLMM_rep_check_scripts/p800/analyses/MegaLMM_part2.R"); reset()
source("MegaLMM_rep_check_scripts/p800/analyses/MegaLMM_part3.R"); reset()
source("MegaLMM_rep_check_scripts/p800/analyses/MegaLMM_part4.R"); reset()
source("MegaLMM_rep_check_scripts/p800/analyses/MegaLMM_part5.R"); reset()
# This script simply merges the results produced by the 5 scripts above.
source("MegaLMM_rep_check_scripts/p800/merge_results.R"); reset()


## p800_lowrank ================================================================
### Step I: Run part 1 through 5 ===============================================
# The 5 scripts below perform the second, independent set of MegaLMM analyses for the
# simulated data with low-rank residual structure (p800_lowrank).
source("MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part1.R"); reset()
source("MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part2.R"); reset()
source("MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part3.R"); reset()
source("MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part4.R"); reset()
source("MegaLMM_rep_check_scripts/p800_lowrank/analyses/MegaLMM_part5.R"); reset()
# This script simply merges the results produced by the 5 scripts above.
source("MegaLMM_rep_check_scripts/p800_lowrank/merge_results.R"); reset()



# Reproducibility checks =======================================================
# The code below check the mean accuracies for MegaLMM reported in the manuscript
# against those produced by the second, independent set of analyses. All of these
# lines of code should return numbers that are very small.

## Hyper HEAT dataset ==========================================================
mean(read.csv("hyper_1415HEAT/results/splines/12a_hyper_results_megalmm_CV1.csv")$acc) -
  mean(read.csv("MegaLMM_rep_check_scripts/hyper_1415HEAT/results/splines/12a_hyper_results_megalmm_CV1.csv")$acc)

mean(read.csv("hyper_1415HEAT/results/splines/12b_hyper_results_megalmm_CV2.csv")$acc) -
  mean(read.csv("MegaLMM_rep_check_scripts/hyper_1415HEAT/results/splines/12b_hyper_results_megalmm_CV2.csv")$acc)

mean(read.csv("hyper_1415HEAT/results/nosplines/12b_hyper_results_megalmm_CV2VEG.csv")$acc) -
  mean(read.csv("MegaLMM_rep_check_scripts/hyper_1415HEAT/results/nosplines/12b_hyper_results_megalmm_CV2VEG.csv")$acc)



## Hyper B5IR dataset ==========================================================
mean(read.csv("hyper_1415B5IR/results/splines/12a_hyper_results_megalmm_CV1.csv")$acc) -
  mean(read.csv("MegaLMM_rep_check_scripts/hyper_1415B5IR/results/splines/12a_hyper_results_megalmm_CV1.csv")$acc)

mean(read.csv("hyper_1415B5IR/results/splines/12b_hyper_results_megalmm_CV2.csv")$acc) -
  mean(read.csv("MegaLMM_rep_check_scripts/hyper_1415B5IR/results/splines/12b_hyper_results_megalmm_CV2.csv")$acc)

mean(read.csv("hyper_1415B5IR/results/VEGsplines/12b_hyper_results_megalmm_CV2VEG.csv")$acc) -
  mean(read.csv("MegaLMM_rep_check_scripts/hyper_1415B5IR/results/VEGsplines/12b_hyper_results_megalmm_CV2VEG.csv")$acc)



## p800 ========================================================================
for (h2s in c("05", "07", "09")) {
  files <- list.files(sprintf("MegaLMM_rep_check_scripts/p800/results/h2s%s/", h2s))
  for (f in files) {
    r1 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800/results/h2s%s/%s", h2s, f))$acc
    r2 <- read.csv(sprintf("p800/results/h2s%s/%s", h2s, f))$acc
    cat(sprintf("%f\n", mean(r1) - mean(r2)))
  }
}



## p800_lowrank ================================================================
for (h2s in c("05", "07", "09")) {
  files <- list.files(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/", h2s))
  for (f in files) {
    r1 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/%s", h2s, f))$acc
    r2 <- read.csv(sprintf("p800_lowrank/results/h2s%s/%s", h2s, f))$acc
    cat(sprintf("%f\n", mean(r1) - mean(r2)))
  }
}

