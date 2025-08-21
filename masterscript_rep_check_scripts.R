# This is one of the 7 masterscripts in the repository.
#
# This particular masterscript sources all R scripts to run some example reproducibility
# spot checks for the different datasets. All of the lines below should print `TRUE` at
# the end.
#
# The annotated scripts themselves (in the `rep_check_scripts/` folder), along with
# the comments in all other analysis scripts, also act as a guide for interested people
# to perform other reproducibility spot checks.
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
if (!("gfBLUPold" %in% installed.packages()[,1])) {
  install.packages("gfBLUPold_1.3.1.tar.gz", type = "source", repos = NULL)
}



# Example reproducibility spot checks ==========================================
# This checks the reproducibility of a specific hyperspectral B5IR datafile.
source("rep_check_scripts/rep_check_generate_hyper_datasets_B5IR.R"); reset()

# This checks the reproducibility of a specific random residual structure simulation
# (p800) datafile.
source("rep_check_scripts/rep_check_generate_p800_datasets.R"); reset()

# This checks the reproducibility of a specific hyperspectral B5IR glfBLUP CV1/CV2
# intermediate result.
source("rep_check_scripts/rep_check_glfBLUP_B5IR.R"); reset()

# This checks the reproducibility of a specific hyperspectral HEAT glfBLUP CV2VEG
# intermediate result.
source("rep_check_scripts/rep_check_glfBLUP_CV2VEG_HEAT.R"); reset()

# This checks the reproducibility of a specific random residual structure simulation
# (p800) glfBLUP CV1/CV2 intermediate result.
source("rep_check_scripts/rep_check_glfBLUP_p800.R"); reset()

# This checks the reproducibility of a specific hyperspectral B5IR lsBLUP CV2
# intermediate result.
source("rep_check_scripts/rep_check_lsBLUP_B5IR.R"); reset()

# This checks the reproducibility of a specific low-rank residual structure simulation
# (p800) lsBLUP CV1 intermediate result.
source("rep_check_scripts/rep_check_lsBLUP_p800_lowrank.R"); reset()

# This checks the reproducibility of a specific hyperspectral HEAT siBLUP CV1
# intermediate result.
source("rep_check_scripts/rep_check_siBLUP_HEAT.R"); reset()

# This checks the reproducibility of a specific low-rank residual structure simulation
# (p800) siBLUP CV1 intermediate result.
source("rep_check_scripts/rep_check_siBLUP_part2_p800_lowrank.R"); reset()

# This checks the reproducibility of a specific random residual structure simulation
# (p800) training/test set division.
source("rep_check_scripts/rep_check_traintest_p800_datasets.R"); reset()




