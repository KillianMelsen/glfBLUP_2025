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
CVs <- c("CV1")

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
        part_1 <- read.csv(sprintf("p800/results/h2s%s/8%s_p800_results_multiMLP_%s_h2y%s_comm%s_h2s%s_1to10.csv",
                                    h2s, lab, CV, h2y, comm, h2s))
        
        part_2 <- read.csv(sprintf("p800/results/h2s%s/8%s_p800_results_multiMLP_%s_h2y%s_comm%s_h2s%s_11to20.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
        # Merging:
        total <- rbind(part_1, part_2)
        
        # Saving the full results:
        write.csv(total, sprintf("p800/results/h2s%s/8%s_p800_results_multiMLP_%s_h2y%s_comm%s_h2s%s.csv",
                                 h2s, lab, CV, h2y, comm, h2s))
        
        progress <- progress + 1
      }
    }
  }
}


