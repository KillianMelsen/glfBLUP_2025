# This script is a copy of the script at `p800_lowrank/analyses/siBLUP_part2.R`.
# The purpose of this script is to provide an annotated example of how to perform
# a reproducibility spot check. All code that is not needed for the spot check has
# been commented out.
# Please carefully read the instructions below, and make sure you manually define
# the variables as explained below.

#!/usr/bin/env Rscript
# This script produces intermediate results 51 to 100 for the siBLUP CV1 and CV2
# analyses on the simulated datasets with low-rank residual structure (p800_lowrank).
#
# The script is structured as follows:
# - There are two outer loops: one for the CVs (CV1 and CV2), and a second one
#   for the focal trait heritabilities ("01", "03", "05", "07", "09").
# - Then there is a third loop over the 3 x 3 = 9 combinations of the
#   communalities ("02", "05", "08") and secondary feature heritabilities ("05", "07", "09").
# - The combination of communality and secondary feature heritability produced by each iteration
#   of the third are defined on lines 64 and 65.
# - A seed is then set inside each iteration of the third loop and the replicates 51 to 100 for the specific
#   combination of CV, secondary feature heritability, communality, and focal trait
#   heritability are produced in a final inner loop starting on line 80.
#
# Suppose you want to reproduce the result of replication 59 for `CV = "CV1"`, `lab = "a"`,
# `h2s = "07"`, `h2y = "01"`, and `comm = "02"`. The first step is to set these values
# manually instead of relying on the two outer loops (CV and h2y), and assignments on
# lines 64 (comm) and 65 (h2s). Then we set `last <- 59` instead of setting it to
# 100 on line 69. Then running the inner loop will produce the correct accuracy for
# replication 59 as the 9th value in `acc`.
# This can be compared to the 59th accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[59],
#           acc[9])
#
# [1] TRUE

# Define these manually:
CV = "CV1"
lab = "a"
h2s = "07"
h2y = "01"
comm = "02"
last <- 59

# Loading libraries:
library(gfBLUPold)
library(rlist)
library(tictoc)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading kinship and marker data:
load("genotypes/K_sim.RData")

# Simulated genetic parameters:
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")
h2.foc <- c("01", "03", "05", "07", "09")
CVs <- c("CV1", "CV2")

combis <- length(h2.foc) * length(CVs) * length(h2.sec) * length(comms)
combi <- 1

par.combis <- data.frame(comm = rep(comms, each = length(h2.sec)),
                         h2s = rep(h2.sec, length(comms)))

# All simulated data for both CV1 and CV2:
# tic("siBLUP")
# for (CV in CVs) {
  # for (h2y in h2.foc) {
    
    # for (i in 1:nrow(par.combis)) {
      
      # comm <- par.combis[i, "comm"]
      # h2s <- par.combis[i, "h2s"]
      
      # Number of simulated datasets to load:
      first <- 51
      # last <- 100
      n.sim <- length(first:last)
      
      # Setting up result storage:
      acc <- numeric(n.sim)
      pen <- numeric(n.sim)
      extra <- vector("list", n.sim)
      
      # Running (SET SEED IN EACH PARALLEL WORKER!):
      set.seed(1997)
      j <- 1
      for (sim in first:last) {
        
        cat(sprintf("Running siBLUP on %s p800_lowrank dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                    CV, j, n.sim, h2s, comm, h2y, combi, combis))
        
        # Loading simulated dataset:
        datalist <- list.load(file = sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
        
        # Storing data and prediction target:
        d <- datalist$data.real
        pred.target <- datalist$pred.target
        
        ### Model ############################################################
        tic(j)
        RESULT <- gfBLUPold::siBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, t.RF = 1, verbose = FALSE)
        toc(log = TRUE)
        ######################################################################
        
        RESULT$preds <- RESULT$preds[match(names(pred.target), names(RESULT$preds))]
        
        acc[j] <- cor(RESULT$preds, pred.target)
        pen[j] <- RESULT$regPen
        extra[[j]] <- list(gamma = RESULT$gamma,
                           Vg = RESULT$Vg,
                           Ve = RESULT$Ve,
                           H2s = RESULT$H2s,
                           penSI = RESULT$penSI,
                           AM.direct = RESULT$AM.direct,
                           AM.indirect = RESULT$AM.indirect)
        
        j <- j + 1
      }
     
cat(paste0(all.equal(read.csv(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[59],
                     acc[9]), "\n"))
      
# We don't need the code below for the reproducibility spot check.
#       # Retrieve computational times:
#       tictoc.logs <- tic.log(format = FALSE)
#       tic.clearlog()
#       comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
#       
#       # Collect results:
#       results <- data.frame(acc = acc,
#                             comptimes = comptimes,
#                             pen = pen)
#       
#       # Making correct CV label:
#       if (CV == "CV1") {
#         lab <- "a"
#       } else if (CV == "CV2") {
#         lab <- "b"
#       }
#       
#       # Export results:
#       write.csv(results, sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.csv",
#                                  h2s, lab, CV, h2y, comm, h2s))
#       
#       list.save(extra, file = sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.RData",
#                                       h2s, lab, CV, h2y, comm, h2s))
#       
#       combi <- combi + 1
#     }
#   }
# }
# toc()

