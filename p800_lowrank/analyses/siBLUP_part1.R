#!/usr/bin/env Rscript
# This script produces intermediate results 1 to 50 for the siBLUP CV1 and CV2
# analyses on the simulated datasets with low-rank residual structure (p800_lowrank).
#
# The script is structured as follows:
# - There are two outer loops: one for the CVs (CV1 and CV2), and a second one
#   for the focal trait heritabilities ("01", "03", "05", "07", "09").
# - Then there is a third loop over the 3 x 3 = 9 combinations of the
#   communalities ("02", "05", "08") and secondary feature heritabilities ("05", "07", "09").
# - The combination of communality and secondary feature heritability produced by each iteration
#   of the third are defined on lines 64 and 65.
# - A seed is then set inside each iteration of the third loop and the replicates 1 to 50 for the specific
#   combination of CV, secondary feature heritability, communality, and focal trait
#   heritability are produced in a final inner loop starting on line 80.
#
# Suppose you want to reproduce the result of replication 3 for `CV = "CV2"`, `lab = "b"`,
# `h2s = "09"`, `h2y = "09"`, and `comm = "05"`. The first step is to set these values
# manually instead of relying on the two outer loops (CV and h2y), and assignments on
# lines 64 (comm) and 65 (h2s). Then we set `last <- 3` instead of setting it to
# 50 on line 69. Then running the inner loop will produce the correct accuracy for
# replication 3 as the 3rd value in `acc`.
# This can be compared to the 3rd accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[3],
#           acc[3])
#
# [1] TRUE

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
tic("siBLUP")
for (CV in CVs) {
  for (h2y in h2.foc) {
    
    for (i in 1:nrow(par.combis)) {
      
      comm <- par.combis[i, "comm"]
      h2s <- par.combis[i, "h2s"]
        
      # Number of simulated datasets to load:
      first <- 1
      last <- 50
      n.sim <- length(first:last)
      
      # Setting up result storage:
      acc <- numeric(n.sim)
      pen <- numeric(n.sim)
      extra <- vector("list", n.sim)
      
      # Running (SET SEED IN EACH LOOP ITERATION!):
      set.seed(1997)
      
      for (sim in first:last) {
        
        cat(sprintf("Running siBLUP on %s p800_lowrank dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                    CV, sim, n.sim, h2s, comm, h2y, combi, combis))
        
        # Loading simulated dataset:
        datalist <- list.load(file = sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
        
        # Storing data and prediction target:
        d <- datalist$data.real
        pred.target <- datalist$pred.target
        
        ### Model ############################################################
        tic(sim)
        RESULT <- gfBLUPold::siBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, t.RF = 1, verbose = FALSE)
        toc(log = TRUE)
        ######################################################################
        
        RESULT$preds <- RESULT$preds[match(names(pred.target), names(RESULT$preds))]
        
        acc[sim] <- cor(RESULT$preds, pred.target)
        pen[sim] <- RESULT$regPen
        extra[[sim]] <- list(gamma = RESULT$gamma,
                             Vg = RESULT$Vg,
                             Ve = RESULT$Ve,
                             H2s = RESULT$H2s,
                             penSI = RESULT$penSI,
                             AM.direct = RESULT$AM.direct,
                             AM.indirect = RESULT$AM.indirect)
        
        
      }
      
      # Retrieve computational times:
      tictoc.logs <- tic.log(format = FALSE)
      tic.clearlog()
      comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
      
      # Collect results:
      results <- data.frame(acc = acc,
                            comptimes = comptimes,
                            pen = pen)
      
      # Making correct CV label:
      if (CV == "CV1") {
        lab <- "a"
      } else if (CV == "CV2") {
        lab <- "b"
      }
      
      # Export results:
      write.csv(results, sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.csv",
                                 h2s, lab, CV, h2y, comm, h2s))

      list.save(extra, file = sprintf("p800_lowrank/results/h2s%s/5%s_p800_lowrank_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.RData",
                                      h2s, lab, CV, h2y, comm, h2s))
      combi <- combi + 1
    }
  }
}
toc()

