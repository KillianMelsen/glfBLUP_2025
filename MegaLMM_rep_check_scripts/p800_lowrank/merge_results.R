# This script simply merges the 5 parts of the second, independent set of
# intermediate results for the MegaLMM p800_lowrank (simulated data with low-rank
# residual structure) analyses. It then deletes the non-merged csv files.

# Loading libraries:
library(rlist)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Simulated genetic parameters:
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")
h2.foc <- c("01", "03", "05", "07", "09")

# Scenarios:
CVs <- c("CV1", "CV2")

# For printing progress:
progress <- 1

for (CV in CVs) {
  for (h2s in h2.sec) {
    for (comm in comms) {
      for (h2y in h2.foc) {
        
        print(progress)
        
        # Making correct CV label:
        if (CV == "CV1") {
          lab <- "a"
        } else if (CV == "CV2") {
          lab <- "b"
        }
        
        # Load part 1, 2, 3, 4, and 5:
        part_1 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_1to4.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        part_2 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_5to8.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        part_3 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_9to12.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        part_4 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_13to16.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        part_5 <- read.csv(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_17to20.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        # Merging:
        total <- rbind(part_1, part_2, part_3, part_4, part_5)
        
        # Saving the full results:
        write.csv(total, sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s.csv",
                                 h2s, lab, CV, h2y, comm, h2s))
        
        # Deleting part 1 - 5:
        unlink(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_1to4.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_5to8.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_9to12.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_13to16.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("MegaLMM_rep_check_scripts/p800_lowrank/results/h2s%s/12%s_p800_lowrank_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_17to20.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        
        progress <- progress + 1
      }
    }
  }
}


