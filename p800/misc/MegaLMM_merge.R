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
        
        # Load part 1, 2, 3, and 4:
        part_1 <- read.csv(sprintf("p800/results/h2s%s/12%s_p800_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_1to6.csv",
                                    h2s, lab, CV, h2y, comm, h2s))
        
        part_2 <- read.csv(sprintf("p800/results/h2s%s/12%s_p800_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_7to12.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        part_3 <- read.csv(sprintf("p800/results/h2s%s/12%s_p800_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_13to16.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        part_4 <- read.csv(sprintf("p800/results/h2s%s/12%s_p800_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_17to20.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        # Merging:
        total <- rbind(part_1, part_2, part_3, part_4)
        
        # Saving the full results:
        write.csv(total, sprintf("p800/results/h2s%s/12%s_p800_results_MegaLMM_%s_h2y%s_comm%s_h2s%s.csv",
                                 h2s, lab, CV, h2y, comm, h2s))
        
        progress <- progress + 1
      }
    }
  }
}


