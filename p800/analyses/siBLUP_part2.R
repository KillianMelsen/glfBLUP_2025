# This script produces intermediate results 51 to 100 for the siBLUP CV1 and CV2
# analyses on the simulated datasets with random residual structure (p800).
#
# The script is structured as follows:
# - There are two outer loops: one for the CVs (CV1 and CV2), and a second one
#   for the focal trait heritabilities ("01", "03", "05", "07", "09").
# - Then there is a parallel loop dividing the 3 x 3 = 9 combinations of the
#   communalities ("02", "05", "08") and secondary feature heritabilities ("05", "07", "09")
#   over 9 parallel workers.
# - The combination of communality and secondary feature heritability produced by each worker
#   are defined on lines 71 and 72.
# - A seed is then set inside each worker and the replicates 1 to 50 for the specific
#   combination of CV, secondary feature heritability, communality, and focal trait
#   heritability are produced in a final inner loop starting on line 87.
#
# Suppose you want to reproduce the result of replication 52 for `CV = "CV2"`, `lab = "b"`,
# `h2s = "09"`, `h2y = "05"`, and `comm = "08"`. The first step is to set these values
# manually instead of relying on the two outer loops (CV and h2y), and assignments on
# lines 71 (comm) and 72 (h2s). Then we set `last <- 52` instead of setting it to
# 100 on line 76. Then running the inner loop will produce the correct accuracy for
# replication 52 as the 2nd value in `acc`.
# This can be compared to the 52nd accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("p800/results/h2s%s/5%s_p800_results_siblup_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[52],
#           acc[2])
#
# [1] TRUE

# Loading libraries:
library(gfBLUPold)
library(rlist)
library(tictoc)
library(doParallel)

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

combis <- length(h2.foc) * length(CVs)
combi <- 1

par.combis <- data.frame(comm = rep(comms, each = length(h2.sec)),
                         h2s = rep(h2.sec, length(comms)))

# All simulated data for both CV1 and CV2:
tic("siBLUP")
for (CV in CVs) {
  for (h2y in h2.foc) {
    
    tic(sprintf("siBLUP combi %d / %d, (CV = %s, h2y = %s)", combi, combis, CV, h2y))
    
    cl <- parallel::makeCluster(9, outfile = sprintf("logs/siBLUP_sim_p800_h2y%s_%s.txt", h2y, CV))
    registerDoParallel(cl)
    
    invisible(
    foreach::foreach(i = 1:nrow(par.combis), .packages = c("rlist", "tictoc")) %dopar% {
      
      comm <- par.combis[i, "comm"]
      h2s <- par.combis[i, "h2s"]
        
      # Number of simulated datasets to load:
      first <- 51
      last <- 100
      n.sim <- length(first:last)
      
      # Setting up result storage:
      acc <- numeric(n.sim)
      pen <- numeric(n.sim)
      extra <- vector("list", n.sim)
      
      # Running (SET SEED IN EACH PARALLEL WORKER!):
      set.seed(1997)
      i <- 1
      for (sim in first:last) {
        
        cat(sprintf("Running siBLUP on %s p800 dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                    CV, i, n.sim, h2s, comm, h2y, combi, combis))
        
        # Loading simulated dataset:
        datalist <- list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
        
        # Storing data and prediction target:
        d <- datalist$data.real
        pred.target <- datalist$pred.target
        
        ### Model ############################################################
        tic(sim)
        RESULT <- gfBLUPold::siBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, t.RF = 1, verbose = FALSE)
        toc(log = TRUE)
        ######################################################################
        
        RESULT$preds <- RESULT$preds[match(names(pred.target), names(RESULT$preds))]
        
        acc[i] <- cor(RESULT$preds, pred.target)
        pen[i] <- RESULT$regPen
        extra[[i]] <- list(gamma = RESULT$gamma,
                           Vg = RESULT$Vg,
                           Ve = RESULT$Ve,
                           H2s = RESULT$H2s,
                           penSI = RESULT$penSI,
                           AM.direct = RESULT$AM.direct,
                           AM.indirect = RESULT$AM.indirect)
        
        i <- i + 1
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
      write.csv(results, sprintf("p800/results/h2s%s/5%s_p800_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.csv",
                                 h2s, lab, CV, h2y, comm, h2s))
      
      list.save(extra, file = sprintf("p800/results/h2s%s/5%s_p800_extra_results_siblup_%s_h2y%s_comm%s_h2s%s_51to100.RData",
                                      h2s, lab, CV, h2y, comm, h2s))
        
    })
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    toc()
    combi <- combi + 1
  }
}
toc()

