# This script simply merges the intermediate results produced by the part 1 and
# part 2 scripts for the p800_lowrank siBLUP analyses. It then deletes the non-merged
# csv files.

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
        
        # Load part 1 and part 2:
        part_1 <- read.csv(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.csv",
                                    h2s, lab, CV, h2y, comm, h2s))
        part_2 <- read.csv(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.csv",
                                   h2s, lab, CV, h2y, comm, h2s))

        # Merging:
        total <- rbind(part_1, part_2)

        # Saving the full results:
        write.csv(total, sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s.csv",
                                 h2s, lab, CV, h2y, comm, h2s))

        # Delete part 1 and 2:
        unlink(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.csv",
                       h2s, lab, CV, h2y, comm, h2s))
        
        # Doing the same for the additional results:
        part_1 <- list.load(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.RData",
                                    h2s, lab, CV, h2y, comm, h2s))
        part_2 <- list.load(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.RData",
                                    h2s, lab, CV, h2y, comm, h2s))
        total <- c(part_1, part_2)
        list.save(total, sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s.RData",
                                 h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.RData",
                       h2s, lab, CV, h2y, comm, h2s))
        unlink(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.RData",
                       h2s, lab, CV, h2y, comm, h2s))
        
        progress <- progress + 1
      }
    }
  }
}


